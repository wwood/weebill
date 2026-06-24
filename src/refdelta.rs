//! Reference-based delta compression of read (sample) sketches.
//!
//! A read sketch from an organism present in a reference database shares most of
//! its FracMinHash hashes with that organism's genome sketch. Instead of storing
//! those hashes as ~56-bit values, we store *which* reference hashes are present.
//! Only hashes explained by no reference (sequencing errors, strain divergence,
//! unknown organisms) are stored explicitly.
//!
//! ## K-mer dereplicated reference database
//!
//! To make "which reference" unambiguous, every k-mer in the database is
//! assigned to at most one owner, using a two-level (species rep / strain)
//! scheme that is robust to strains contaminated with sequence from other
//! species:
//!   * Build order places strains of a species *contiguously*, species
//!     representatives first.
//!   * A k-mer present in exactly one species **representative** is owned by that
//!     rep (representatives are higher quality, so they win over strains — a
//!     contaminant k-mer that really belongs to another species' rep is
//!     attributed there, not to the strain carrying the contamination).
//!   * By default, a k-mer in three or more reps is ambiguous → the shared
//!     **pool**. `ref-build --pool-min-genomes` can change that threshold; with
//!     3, exactly two same-tier genomes keep the k-mer assigned to the first
//!     owner.
//!   * A k-mer in no rep but exactly one strain is owned by that strain.
//!   * A k-mer in no rep and enough strains to meet the threshold → the shared
//!     pool.
//!
//! ## Two-stage reference database layout (`.sylref`)
//!
//! Loading a whole reference into a single hash→location table is expensive for
//! large databases. The `.sylref` file is therefore seekable and read in two
//! stages:
//!   * **Stage 1 (sparse, loaded fully):** a 1/N subsample of every genome's
//!     *distinctive* k-mers, stored uncompressed, mapping each sparse hash to its
//!     owning genome. Cheap to load; querying a sample against it yields the small
//!     set of genomes the sample actually contains ("hit genomes").
//!   * **Stage 2 (dense, loaded on demand):** each genome's distinctive set,
//!     minus the hashes already held in stage 1, is an independently Golomb–Rice-
//!     coded block at a known offset (so each k-mer is stored only once). Only the
//!     hit genomes' blocks are decoded and merged back with their stage-1 hashes
//!     to reconstruct the full distinctive array (cached across samples). The
//!     shared **pool** is the one exception — it is conserved across samples and
//!     is loaded once when the index is opened.
//!
//! Missing a genome in stage 1 only costs compression ratio (its hashes fall back
//! to "novel" full-price coding), never correctness: decompression reads the exact
//! genome ids recorded in the sample and loads precisely those dense blocks.
//!
//! ## Encoding a read sketch
//!
//! Read hashes are partitioned into (a) per-genome distinctive hits, (b) pool
//! hits, (c) novel hashes. Each matched set is encoded as positions into the
//! relevant sorted reference array using whichever of {bitmask, present-Rice,
//! absent-Rice} is smallest (a well-covered genome costs ~0 bytes via an empty
//! absent set). The hit genome ids are themselves delta-coded: because strains of
//! a species are contiguous, a sample's hits cluster into small id ranges.
//! Novel hashes use Golomb–Rice (as for normal sketches); counts are varints.
//! The whole payload is zstd-framed.

use crate::cmdline::{RefBuildArgs, RefCompressArgs};
use crate::compress::{
    self, read_hashes, read_seq_sketch_compressed, read_string, read_uvarint, write_hashes,
    write_string, write_uvarint,
};
use crate::constants::*;
use crate::seeding::{mm_hash64, rev_mm_hash64};
use crate::types::*;
use boomphf::Mphf;
use fxhash::{FxHashMap, FxHashSet};
use log::*;
use rayon::prelude::*;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};
use std::time::Instant;

const REFDB_MAGIC: &[u8; 4] = b"SYLR";
/// v5 adds an optional per-rep-genome 2-bit nucleotide section (`--store-genomes`)
/// used by error-k-mer encoding. v4 files (no genomes) are still readable.
const REFDB_VERSION: u8 = 5;
const REFDB_MIN_VERSION: u8 = 4;
const SKETCH_MAGIC: &[u8; 4] = b"SYLD"; // reference-Delta sample
/// v4 adds the single-substitution error-k-mer section. v3 samples still decode.
const SKETCH_VERSION: u8 = 4;

const SCHEME_BITMASK: u8 = 0;
const SCHEME_PRESENT_RICE: u8 = 1;
const SCHEME_ABSENT_RICE: u8 = 2;

/// zstd level for the read-sketch payload (matches the normal sketch format).
const ZSTD_LEVEL: i32 = 3;

/// Fixed size of the `.sylref` header: magic (4) + version (1) + footer offset (8).
const HEADER_LEN: u64 = 13;

/// A genome is treated as a "hit" (its dense block is loaded) once the sample
/// shares at least this many of the genome's stage-1 sparse k-mers. Distinctive
/// k-mers are genome-specific, so even 1 is a meaningful signal; over-loading a
/// genome only costs time, under-loading only costs compression ratio.
const SPARSE_MIN_HITS: u32 = 1;
const SPARSE_MPHF_GAMMA: f64 = 2.0;
const POOL_MPHF_GAMMA: f64 = 2.0;

/// Metadata for one reference genome, in the dereplicated build order.
#[derive(Clone, Debug, PartialEq)]
pub struct RefGenome {
    pub file_name: String,
    pub species: String,
    pub is_rep: bool,
}

/// A k-mer dereplicated reference database: per-genome distinctive hash sets plus
/// one shared pool, with genomes ordered so a species' strains are contiguous.
/// This is the in-memory form produced by `build_refdb` and serialized by
/// `write_refdb`; querying/compressing uses the seekable `RefIndex` instead.
#[derive(Clone, Debug, PartialEq)]
pub struct RefDb {
    pub c: usize,
    pub k: usize,
    pub genomes: Vec<RefGenome>,
    /// `distinctive[g]` = sorted, deduplicated hashes owned uniquely by genome `g`.
    pub distinctive: Vec<Vec<u64>>,
    /// Shared-pool hashes (after rep preference), sorted and deduplicated.
    pub pool: Vec<u64>,
    /// Optional nucleotide-sequence source per genome (only set for species
    /// representatives, used by error-k-mer encoding). `rep_seqs[g] == None` for
    /// a genome with no stored sequence. An empty vec means "no genomes stored".
    /// At write time each source is read/packed and streamed to disk one genome
    /// at a time, so building a reference over hundreds of thousands of genomes
    /// never holds every sequence in RAM.
    pub rep_seqs: Vec<Option<GenomeSource>>,
}

/// Where a genome's nucleotide sequence comes from when writing the reference.
#[derive(Clone, Debug, PartialEq)]
pub enum GenomeSource {
    /// Already-decoded sequence (used by tests and small inputs).
    InMemory(GenomeSeq),
    /// A FASTA path read and 2-bit packed on demand at write time (streaming).
    Path(String),
}

/// A genome's nucleotide sequence, one `Vec` of 2-bit base codes (0=A,1=C,2=G,
/// 3=T; any non-ACGT base maps to 0, exactly as the sketcher does) per contig.
/// K-mers never span contigs.
#[derive(Clone, Debug, PartialEq)]
pub struct GenomeSeq {
    pub contigs: Vec<Vec<u8>>,
}

#[derive(Clone, Debug)]
pub struct RefCompressTelemetry {
    pub sample_name: String,
    pub ref_db_name: String,
    pub ref_screen_ani: f64,
    pub hit_genomes_total: usize,
    pub genome_id: u32,
    pub genome_file: String,
    pub species: String,
    pub sparse_hits: u32,
    pub sparse_total: u32,
    pub sparse_ani: f64,
    pub assigned_kmers: usize,
    pub error_kmers: usize,
}

impl GenomeSeq {
    /// Cumulative base offset of each contig (so a global base index uniquely maps
    /// back to a single contig and local position). `cum[i]` = sum of contig
    /// lengths < `i`; the final entry is the total base count.
    fn cumulative(&self) -> Vec<u64> {
        let mut cum = Vec::with_capacity(self.contigs.len() + 1);
        let mut acc = 0u64;
        cum.push(0);
        for c in &self.contigs {
            acc += c.len() as u64;
            cum.push(acc);
        }
        cum
    }

    /// Map a global base index (a k-mer start) to `(contig, local_start)`.
    fn locate(&self, cum: &[u64], global: u64) -> (usize, usize) {
        // contigs are few; a linear scan with the precomputed prefix sums is fine.
        let mut ci = 0;
        while ci + 1 < cum.len() && cum[ci + 1] <= global {
            ci += 1;
        }
        (ci, (global - cum[ci]) as usize)
    }
}

/// Forward and reverse-complement 2-bit packings of the k-mer at `codes[start..start+k]`,
/// matching `seeding::fmh_seeds` exactly: the forward k-mer has base 0 in the high
/// bits, the reverse complement is built from the per-base complements.
#[inline]
fn window_fr(codes: &[u8], start: usize, k: usize) -> (u64, u64) {
    let mut f = 0u64;
    let mut r = 0u64;
    for j in 0..k {
        let nuc_f = codes[start + j] as u64;
        let nuc_r = 3 - nuc_f;
        f = (f << 2) | nuc_f;
        r |= nuc_r << (2 * j);
    }
    (f, r)
}

/// The canonical (min of forward / reverse-complement) FracMinHash hash of the
/// k-mer at `start` with the base at offset `off` replaced by `base`. Mirrors the
/// read sketcher so a reconstructed error k-mer hashes identically.
#[inline]
fn substituted_hash(f: u64, r: u64, k: usize, off: usize, base: u8) -> u64 {
    let shift_f = 2 * (k - 1 - off);
    let shift_r = 2 * off;
    let f2 = (f & !(3u64 << shift_f)) | ((base as u64) << shift_f);
    let r2 = (r & !(3u64 << shift_r)) | (((3 - base) as u64) << shift_r);
    mm_hash64(f2.min(r2))
}

// --- genome sequence (de)serialization (2-bit packed) -----------------------

/// Read a FASTA/FASTQ into a `GenomeSeq` of 2-bit base codes (one code per byte),
/// matching the sketcher's `BYTE_TO_SEQ` mapping. Returns `None` if unreadable.
fn read_genome_seq_from_fasta(path: &str) -> Option<GenomeSeq> {
    let reader = needletail::parse_fastx_file(path);
    let mut reader = match reader {
        Ok(r) => r,
        Err(_) => return None,
    };
    let mut contigs = Vec::new();
    while let Some(record) = reader.next() {
        let record = match record {
            Ok(r) => r,
            Err(_) => return None,
        };
        let seq = record.seq();
        let codes: Vec<u8> = seq.iter().map(|&b| BYTE_TO_SEQ[b as usize]).collect();
        contigs.push(codes);
    }
    Some(GenomeSeq { contigs })
}

/// Serialize a genome's contigs as: contig count, then per contig its base length
/// followed by the 2-bit-packed bases (4 bases/byte, contig-aligned).
fn write_genome_seq_block<W: Write>(w: &mut W, seq: &GenomeSeq) -> io::Result<()> {
    write_uvarint(w, seq.contigs.len() as u64)?;
    for contig in &seq.contigs {
        write_uvarint(w, contig.len() as u64)?;
        let mut byte = 0u8;
        let mut nfilled = 0u8;
        for &code in contig {
            byte |= (code & 3) << (2 * nfilled);
            nfilled += 1;
            if nfilled == 4 {
                w.write_all(&[byte])?;
                byte = 0;
                nfilled = 0;
            }
        }
        if nfilled > 0 {
            w.write_all(&[byte])?;
        }
    }
    Ok(())
}

/// Materialize a genome source and serialize it into a self-contained block, or
/// `None` if a `Path` source cannot be read (that genome then simply gets no
/// stored sequence and falls back to full-price coding).
fn genome_block_bytes(src: &GenomeSource) -> io::Result<Option<Vec<u8>>> {
    let seq = match src {
        GenomeSource::InMemory(s) => Some(s.clone()),
        GenomeSource::Path(p) => read_genome_seq_from_fasta(p),
    };
    match seq {
        Some(s) => {
            let mut v = Vec::new();
            write_genome_seq_block(&mut v, &s)?;
            Ok(Some(v))
        }
        None => Ok(None),
    }
}

fn read_genome_seq_block<R: Read>(r: &mut R) -> io::Result<GenomeSeq> {
    let n_contigs = read_uvarint(r)? as usize;
    let mut contigs = Vec::with_capacity(n_contigs);
    for _ in 0..n_contigs {
        let len = read_uvarint(r)? as usize;
        let nbytes = (len + 3) / 4;
        let mut packed = vec![0u8; nbytes];
        r.read_exact(&mut packed)?;
        let mut codes = Vec::with_capacity(len);
        for i in 0..len {
            let byte = packed[i / 4];
            codes.push((byte >> (2 * (i % 4))) & 3);
        }
        contigs.push(codes);
    }
    Ok(GenomeSeq { contigs })
}

/// A stable fingerprint so a compressed sample can be matched to its reference DB.
/// This mixes only array lengths and boundary hashes (not full contents): a
/// deliberately different DB with matching shape/boundaries could in principle
/// collide, but the chance of that in practice is negligible and it keeps the
/// digest cheap.
fn fingerprint(db: &RefDb) -> u64 {
    let mut h: u64 = 1469598103934665603; // FNV-1a offset
    let mut mix = |x: u64| {
        h ^= x;
        h = h.wrapping_mul(1099511628211);
    };
    mix(db.c as u64);
    mix(db.k as u64);
    mix(db.genomes.len() as u64);
    mix(db.pool.len() as u64);
    for (i, d) in db.distinctive.iter().enumerate() {
        mix(i as u64);
        mix(d.len() as u64);
        if let Some(&first) = d.first() {
            mix(first);
        }
        if let Some(&last) = d.last() {
            mix(last);
        }
    }
    if let Some(&p) = db.pool.first() {
        mix(p);
    }
    if let Some(&p) = db.pool.last() {
        mix(p);
    }
    h
}

fn owner_for_accum(a: &OwnAccum, pool_min_genomes: u32) -> u32 {
    if a.rep_count > 0 {
        if a.rep_count >= pool_min_genomes {
            POOL
        } else {
            a.rep_id
        }
    } else if a.strain_count >= pool_min_genomes {
        POOL
    } else {
        a.strain_id
    }
}

fn build_pool_mphf(pool: &[u64]) -> io::Result<(Mphf<u64>, Vec<u64>)> {
    let mphf = Mphf::new_parallel(POOL_MPHF_GAMMA, pool, Some(0));
    let mut ordered = vec![0u64; pool.len()];
    let mut filled = vec![false; pool.len()];
    for &h in pool {
        let slot = mphf.hash(&h) as usize;
        if slot >= pool.len() || filled[slot] {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "MPHF construction produced an invalid pool slot",
            ));
        }
        ordered[slot] = h;
        filled[slot] = true;
    }
    Ok((mphf, ordered))
}

fn write_pool_block<W: Write>(w: &mut W, pool: &[u64]) -> io::Result<()> {
    let (mphf, ordered) = build_pool_mphf(pool)?;
    let mphf_bytes = bincode::serialize(&mphf)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e.to_string()))?;
    write_uvarint(w, mphf_bytes.len() as u64)?;
    w.write_all(&mphf_bytes)?;
    for &h in &ordered {
        w.write_all(&h.to_le_bytes())?;
    }
    Ok(())
}

fn read_pool_block<R: Read>(r: &mut R, pool_len: usize) -> io::Result<(Mphf<u64>, Vec<u64>)> {
    let mphf_len = read_uvarint(r)? as usize;
    let mut mphf_bytes = vec![0u8; mphf_len];
    r.read_exact(&mut mphf_bytes)?;
    let mphf: Mphf<u64> = bincode::deserialize(&mphf_bytes)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e.to_string()))?;
    let mut pool = Vec::with_capacity(pool_len);
    for _ in 0..pool_len {
        let mut buf = [0u8; 8];
        r.read_exact(&mut buf)?;
        pool.push(u64::from_le_bytes(buf));
    }
    Ok((mphf, pool))
}

#[inline]
fn sparse_fingerprint(h: u64) -> u32 {
    (h ^ (h >> 32)) as u32
}

fn write_sparse_index_block<W: Write>(w: &mut W, keys: &[u64], owners: &[u32]) -> io::Result<()> {
    if keys.len() != owners.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "sparse index keys and owners have different lengths",
        ));
    }
    let mphf = Mphf::new_parallel(SPARSE_MPHF_GAMMA, keys, Some(0));
    let mut fingerprints = vec![0u32; keys.len()];
    let mut ordered_owners = vec![0u32; keys.len()];
    let mut filled = vec![false; keys.len()];
    for (&h, &owner) in keys.iter().zip(owners.iter()) {
        let slot = mphf.hash(&h) as usize;
        if slot >= keys.len() || filled[slot] {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "MPHF construction produced an invalid sparse slot",
            ));
        }
        fingerprints[slot] = sparse_fingerprint(h);
        ordered_owners[slot] = owner;
        filled[slot] = true;
    }

    let mphf_bytes = bincode::serialize(&mphf)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e.to_string()))?;
    write_uvarint(w, mphf_bytes.len() as u64)?;
    w.write_all(&mphf_bytes)?;
    for &fp in &fingerprints {
        w.write_all(&fp.to_le_bytes())?;
    }
    for &owner in &ordered_owners {
        w.write_all(&owner.to_le_bytes())?;
    }
    Ok(())
}

#[inline]
fn sparse_threshold(sparse_c: usize) -> u64 {
    u64::MAX / sparse_c.max(1) as u64
}

#[inline]
fn keep_sparse_hash(h: u64, sparse_c: usize) -> bool {
    h < sparse_threshold(sparse_c)
}

#[inline]
fn sparse_naive_ani(matches: u32, domain: u32, k: usize) -> f64 {
    if matches == 0 || domain == 0 || k == 0 {
        return 0.0;
    }
    f64::powf(matches as f64 / domain as f64, 1.0 / k as f64) * 100.0
}

fn read_sparse_index_block<R: Read>(
    r: &mut R,
    sparse_len: usize,
) -> io::Result<(Mphf<u64>, Vec<u32>, Vec<u32>)> {
    let mphf_len = read_uvarint(r)? as usize;
    let mut mphf_bytes = vec![0u8; mphf_len];
    r.read_exact(&mut mphf_bytes)?;
    let mphf: Mphf<u64> = bincode::deserialize(&mphf_bytes)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e.to_string()))?;

    let mut fingerprints = Vec::with_capacity(sparse_len);
    for _ in 0..sparse_len {
        let mut buf = [0u8; 4];
        r.read_exact(&mut buf)?;
        fingerprints.push(u32::from_le_bytes(buf));
    }

    let mut owners = Vec::with_capacity(sparse_len);
    for _ in 0..sparse_len {
        let mut buf = [0u8; 4];
        r.read_exact(&mut buf)?;
        owners.push(u32::from_le_bytes(buf));
    }
    Ok((mphf, fingerprints, owners))
}

// --- building the reference DB ---------------------------------------------

struct OwnAccum {
    rep_count: u32,
    rep_id: u32,
    strain_count: u32,
    strain_id: u32,
}

const POOL: u32 = u32::MAX;

/// Build a dereplicated reference DB from genome sketches and an optional
/// `file_name -> (species, is_rep)` taxonomy. Genomes absent from the taxonomy
/// are treated as their own single-genome species representative.
pub fn build_refdb(
    sketches: &[GenomeSketch],
    taxonomy: &FxHashMap<String, (String, bool)>,
) -> RefDb {
    build_refdb_with_pool_min_genomes(sketches, taxonomy, 3)
}

pub fn build_refdb_with_pool_min_genomes(
    sketches: &[GenomeSketch],
    taxonomy: &FxHashMap<String, (String, bool)>,
    pool_min_genomes: u32,
) -> RefDb {
    let pool_min_genomes = pool_min_genomes.max(2);
    let c = sketches.first().map(|s| s.c).unwrap_or(0);
    let k = sketches.first().map(|s| s.k).unwrap_or(0);

    // resolve species/is_rep per input genome (match on full path or basename)
    let resolve = |file_name: &str| -> (String, bool) {
        if let Some(v) = taxonomy.get(file_name) {
            return v.clone();
        }
        let base = std::path::Path::new(file_name)
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or(file_name);
        if let Some(v) = taxonomy.get(base) {
            return v.clone();
        }
        (file_name.to_string(), true)
    };

    // order: species, then reps before strains, then file name -> strains contiguous
    let mut order: Vec<usize> = (0..sketches.len()).collect();
    let meta: Vec<(String, bool)> = sketches.iter().map(|s| resolve(&s.file_name)).collect();
    order.sort_by(|&a, &b| {
        meta[a]
            .0
            .cmp(&meta[b].0)
            .then(meta[b].1.cmp(&meta[a].1)) // is_rep true first
            .then(sketches[a].file_name.cmp(&sketches[b].file_name))
    });

    let genomes: Vec<RefGenome> = order
        .iter()
        .map(|&i| RefGenome {
            file_name: sketches[i].file_name.clone(),
            species: meta[i].0.clone(),
            is_rep: meta[i].1,
        })
        .collect();

    // accumulate ownership over all k-mers
    let mut acc: FxHashMap<u64, OwnAccum> = FxHashMap::default();
    for (new_id, &orig) in order.iter().enumerate() {
        let is_rep = meta[orig].1;
        let new_id = new_id as u32;
        // dedup within a genome first
        let mut seen: FxHashMap<u64, ()> = FxHashMap::default();
        for &h in &sketches[orig].genome_kmers {
            if seen.insert(h, ()).is_some() {
                continue;
            }
            let e = acc.entry(h).or_insert(OwnAccum {
                rep_count: 0,
                rep_id: 0,
                strain_count: 0,
                strain_id: 0,
            });
            if is_rep {
                if e.rep_count == 0 {
                    e.rep_id = new_id;
                }
                e.rep_count += 1;
            } else {
                if e.strain_count == 0 {
                    e.strain_id = new_id;
                }
                e.strain_count += 1;
            }
        }
    }

    // assign owners
    let mut distinctive: Vec<Vec<u64>> = vec![Vec::new(); genomes.len()];
    let mut pool: Vec<u64> = Vec::new();
    for (&h, a) in acc.iter() {
        let owner = owner_for_accum(a, pool_min_genomes);
        if owner == POOL {
            pool.push(h);
        } else {
            distinctive[owner as usize].push(h);
        }
    }
    for d in distinctive.iter_mut() {
        d.sort_unstable();
        d.dedup();
    }
    pool.sort_unstable();
    pool.dedup();

    RefDb {
        c,
        k,
        genomes,
        distinctive,
        pool,
        rep_seqs: Vec::new(),
    }
}

// --- reference DB serialization (seekable two-stage container) --------------

/// Per-genome species id: genomes are already grouped contiguously by species, so
/// a running counter assigns each species a dense integer id.
fn species_ids(genomes: &[RefGenome]) -> Vec<u32> {
    let mut ids = vec![0u32; genomes.len()];
    let mut cur = 0u32;
    for g in 0..genomes.len() {
        if g > 0 && genomes[g].species != genomes[g - 1].species {
            cur += 1;
        }
        ids[g] = cur;
    }
    ids
}

/// Write a `.sylref` in the seekable two-stage layout. `sparse_c` controls the
/// stage-1 sparse FracMinHash rate; hashes below `u64::MAX / sparse_c` are kept.
pub fn write_refdb<W: Write + Seek>(mut w: W, db: &RefDb, sparse_c: usize) -> io::Result<()> {
    let ng = db.genomes.len();
    let sp_ids = species_ids(&db.genomes);
    let sparse_c = sparse_c.max(db.c).max(1);

    // Everything is streamed straight to the writer so neither the serialized body
    // (the sparse/index/pool/dense sections) nor the genome sequences are ever held
    // in RAM in full — important for references over hundreds of thousands of
    // genomes, where the body alone is tens of GB. `pos` tracks the absolute file
    // offset; a small `buf` is reused to serialize one section at a time. The footer
    // offset is only known at the end, so the header keeps a placeholder that is
    // patched once everything is written (needs Seek).
    w.write_all(REFDB_MAGIC)?;
    w.write_all(&[REFDB_VERSION])?;
    w.write_all(&0u64.to_le_bytes())?; // placeholder footer offset
    let mut pos = HEADER_LEN;
    let mut buf: Vec<u8> = Vec::new();

    // Sparse section: FracMinHash subset of each genome's distinctive k-mers,
    // uncompressed. The MPHF below keeps the stage-1 screen RAM-light.
    let sparse_off = pos;
    let mut sparse_count = vec![0usize; ng];
    let mut sparse_keys: Vec<u64> = Vec::new();
    let mut sparse_owners: Vec<u32> = Vec::new();
    for g in 0..ng {
        buf.clear();
        let mut cnt = 0usize;
        for &h in &db.distinctive[g] {
            if keep_sparse_hash(h, sparse_c) {
                buf.extend_from_slice(&h.to_le_bytes());
                sparse_keys.push(h);
                sparse_owners.push(g as u32);
                cnt += 1;
            }
        }
        w.write_all(&buf)?;
        pos += buf.len() as u64;
        sparse_count[g] = cnt;
    }

    let sparse_index_off = pos;
    buf.clear();
    write_sparse_index_block(&mut buf, &sparse_keys, &sparse_owners)?;
    drop(sparse_keys);
    drop(sparse_owners);
    w.write_all(&buf)?;
    pos += buf.len() as u64;

    let pool_off = pos;
    buf.clear();
    write_pool_block(&mut buf, &db.pool)?;
    w.write_all(&buf)?;
    pos += buf.len() as u64;

    // Dense blocks store only the *complement* of the stage-1 sparse subset, so a
    // k-mer kept in stage 1 is not also stored here; load merges the two back.
    let mut dense_off = vec![0u64; ng];
    for g in 0..ng {
        dense_off[g] = pos;
        let complement: Vec<u64> = db.distinctive[g]
            .iter()
            .copied()
            .filter(|&h| !keep_sparse_hash(h, sparse_c))
            .collect();
        buf.clear();
        write_hashes(&mut buf, &complement)?;
        w.write_all(&buf)?;
        pos += buf.len() as u64;
    }

    // Optional per-rep-genome 2-bit nucleotide blocks (loaded on demand for
    // error-k-mer encoding). `seq_off[g] == 0` means no sequence is stored. Blocks
    // are decoded/packed in parallel a chunk at a time and written in order.
    let mut seq_off = vec![0u64; ng];
    if db.rep_seqs.iter().any(|s| s.is_some()) {
        const SEQ_CHUNK: usize = 512;
        let mut g0 = 0usize;
        while g0 < ng {
            let g1 = (g0 + SEQ_CHUNK).min(ng);
            let blocks: Vec<(usize, Option<Vec<u8>>)> = (g0..g1)
                .into_par_iter()
                .map(|g| {
                    let blk = match db.rep_seqs.get(g).and_then(|s| s.as_ref()) {
                        Some(src) => genome_block_bytes(src).unwrap_or(None),
                        None => None,
                    };
                    (g, blk)
                })
                .collect();
            for (g, blk) in blocks {
                if let Some(b) = blk {
                    seq_off[g] = pos;
                    w.write_all(&b)?;
                    pos += b.len() as u64;
                }
            }
            g0 = g1;
        }
    }
    let has_seqs = seq_off.iter().any(|&o| o != 0);
    let footer_offset = pos;

    // Footer (zstd-compressed): metadata + absolute offsets into the file.
    let mut footer: Vec<u8> = Vec::new();
    write_uvarint(&mut footer, db.c as u64)?;
    write_uvarint(&mut footer, db.k as u64)?;
    write_uvarint(&mut footer, sparse_c as u64)?;
    footer.extend_from_slice(&fingerprint(db).to_le_bytes());
    write_uvarint(&mut footer, ng as u64)?;
    write_uvarint(&mut footer, sparse_off)?;
    for g in 0..ng {
        write_uvarint(&mut footer, sparse_count[g] as u64)?;
    }
    write_uvarint(&mut footer, sparse_index_off)?;
    write_uvarint(&mut footer, pool_off)?;
    write_uvarint(&mut footer, db.pool.len() as u64)?;
    for g in 0..ng {
        write_string(&mut footer, &db.genomes[g].file_name)?;
        write_string(&mut footer, &db.genomes[g].species)?;
        footer.push(db.genomes[g].is_rep as u8);
        write_uvarint(&mut footer, sp_ids[g] as u64)?;
        write_uvarint(&mut footer, dense_off[g])?;
        write_uvarint(&mut footer, db.distinctive[g].len() as u64)?;
    }
    // v5: optional genome-sequence directory (a per-genome offset; 0 = none).
    footer.push(has_seqs as u8);
    if has_seqs {
        for g in 0..ng {
            write_uvarint(&mut footer, seq_off[g])?;
        }
    }
    let footer_comp = zstd::stream::encode_all(&footer[..], 9)?;
    w.write_all(&footer_comp)?;

    // Patch the header's footer offset now that the genome section length is known.
    w.flush()?;
    w.seek(SeekFrom::Start(5))?;
    w.write_all(&footer_offset.to_le_bytes())?;
    w.flush()?;
    Ok(())
}

/// Per-genome metadata held by an opened `RefIndex` (everything but the dense
/// distinctive block, which is loaded on demand).
#[derive(Clone, Debug)]
pub struct RefGenomeMeta {
    pub file_name: String,
    pub species: String,
    pub is_rep: bool,
    pub species_id: u32,
    sparse_offset: u64,
    sparse_count: u32,
    dense_offset: u64,
    dense_domain: u32,
    /// Absolute offset of this genome's 2-bit sequence block, or 0 if not stored.
    seq_offset: u64,
}

/// Any source we can both stream and seek (a file or an in-memory cursor).
pub trait ReadSeek: Read + Seek + Send {}
impl<T: Read + Seek + Send> ReadSeek for T {}

/// An opened two-stage reference. Constructing it loads only the stage-1 sparse
/// index and the shared pool; per-genome dense blocks are decoded lazily (and
/// cached) as samples need them.
pub struct RefIndex {
    pub c: usize,
    pub k: usize,
    sparse_c: usize,
    fingerprint: u64,
    pub genomes: Vec<RefGenomeMeta>,
    /// stage-1: sparse distinctive hash -> owning genome id via MPHF slot arrays.
    sparse_mphf: Mphf<u64>,
    sparse_fingerprints: Vec<u32>,
    sparse_owners: Vec<u32>,
    /// shared pool, loaded once.
    pool: Vec<u64>,
    pool_mphf: Mphf<u64>,
    reader: Mutex<Box<dyn ReadSeek>>,
    cache: Mutex<FxHashMap<u32, Arc<Vec<u64>>>>,
    /// decoded 2-bit genome sequences, loaded on demand and cached.
    seq_cache: Mutex<FxHashMap<u32, Arc<GenomeSeq>>>,
    compression_only: bool,
}

impl RefIndex {
    /// Whether this reference stores genome sequences (i.e. supports error-k-mer
    /// encoding). True if any genome has a non-zero sequence offset.
    pub fn has_genome_seqs(&self) -> bool {
        self.genomes.iter().any(|g| g.seq_offset != 0)
    }
}

enum RefOpenMode {
    Full,
    CompressionOnly,
}

/// Open a `.sylref`, loading stage 1 (sparse index) and the shared pool.
pub fn open_ref_index<R: Read + Seek + Send + 'static>(r: R) -> io::Result<RefIndex> {
    open_ref_index_with_mode(r, RefOpenMode::Full)
}

pub fn open_ref_index_for_compress<R: Read + Seek + Send + 'static>(r: R) -> io::Result<RefIndex> {
    open_ref_index_with_mode(r, RefOpenMode::CompressionOnly)
}

fn open_ref_index_with_mode<R: Read + Seek + Send + 'static>(
    mut r: R,
    mode: RefOpenMode,
) -> io::Result<RefIndex> {
    let mut hdr = [0u8; HEADER_LEN as usize];
    r.read_exact(&mut hdr)?;
    if &hdr[0..4] != REFDB_MAGIC {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "not a sylph reference DB",
        ));
    }
    let version = hdr[4];
    if !(REFDB_MIN_VERSION..=REFDB_VERSION).contains(&version) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "unsupported reference DB version",
        ));
    }
    let footer_offset = u64::from_le_bytes(hdr[5..13].try_into().unwrap());

    // footer
    r.seek(SeekFrom::Start(footer_offset))?;
    let mut comp = Vec::new();
    r.read_to_end(&mut comp)?;
    let fbytes = zstd::stream::decode_all(&comp[..])?;
    let mut f = &fbytes[..];
    let c = read_uvarint(&mut f)? as usize;
    let k = read_uvarint(&mut f)? as usize;
    let sparse_c = read_uvarint(&mut f)? as usize;
    let mut fpb = [0u8; 8];
    f.read_exact(&mut fpb)?;
    let fingerprint = u64::from_le_bytes(fpb);
    let ng = read_uvarint(&mut f)? as usize;
    let sparse_offset = read_uvarint(&mut f)?;
    let mut sparse_count = Vec::with_capacity(ng);
    for _ in 0..ng {
        sparse_count.push(read_uvarint(&mut f)? as usize);
    }
    let sparse_index_offset = read_uvarint(&mut f)?;
    let pool_offset = read_uvarint(&mut f)?;
    let pool_domain = read_uvarint(&mut f)? as usize;
    let mut per_genome_sparse_offsets = Vec::with_capacity(ng);
    let mut sparse_pos = sparse_offset;
    for &cnt in &sparse_count {
        per_genome_sparse_offsets.push(sparse_pos);
        sparse_pos += cnt as u64 * 8;
    }
    let mut genomes = Vec::with_capacity(ng);
    for g in 0..ng {
        let file_name = read_string(&mut f)?;
        let species = read_string(&mut f)?;
        let mut rep = [0u8; 1];
        f.read_exact(&mut rep)?;
        let species_id = read_uvarint(&mut f)? as u32;
        let dense_offset = read_uvarint(&mut f)?;
        let dense_domain = read_uvarint(&mut f)? as u32;
        genomes.push(RefGenomeMeta {
            file_name,
            species,
            is_rep: rep[0] != 0,
            species_id,
            sparse_offset: per_genome_sparse_offsets[g],
            sparse_count: sparse_count[g] as u32,
            dense_offset,
            dense_domain,
            seq_offset: 0,
        });
    }

    // v5: optional per-genome sequence-block offset directory.
    if version >= 5 {
        let mut has_seqs = [0u8; 1];
        f.read_exact(&mut has_seqs)?;
        if has_seqs[0] != 0 {
            for g in 0..ng {
                genomes[g].seq_offset = read_uvarint(&mut f)?;
            }
        }
    }

    let total_sparse: usize = sparse_count.iter().sum();
    r.seek(SeekFrom::Start(sparse_index_offset))?;
    let (sparse_mphf, sparse_fingerprints, sparse_owners) =
        read_sparse_index_block(&mut r, total_sparse)?;

    // shared pool (loaded once). Hashes are stored in MPHF slot order, so
    // compression can compute a candidate slot and verify it without a HashMap.
    r.seek(SeekFrom::Start(pool_offset))?;
    let (pool_mphf, pool) = read_pool_block(&mut r, pool_domain)?;

    Ok(RefIndex {
        c,
        k,
        sparse_c,
        fingerprint,
        genomes,
        sparse_mphf,
        sparse_fingerprints,
        sparse_owners,
        pool,
        pool_mphf,
        reader: Mutex::new(Box::new(r)),
        cache: Mutex::new(FxHashMap::default()),
        seq_cache: Mutex::new(FxHashMap::default()),
        compression_only: matches!(mode, RefOpenMode::CompressionOnly),
    })
}

impl RefIndex {
    /// Stage-1 query: sparse hit counts for genomes passing the sparse FracMinHash ANI screen.
    pub fn hit_genome_counts(&self, sketch: &SequencesSketch, min_ani: f64) -> Vec<(u32, u32)> {
        let mut counts: FxHashMap<u32, u32> = FxHashMap::default();
        for &h in sketch.kmer_counts.keys() {
            if !keep_sparse_hash(h, self.sparse_c) {
                continue;
            }
            if let Some(slot) = self.sparse_mphf.try_hash(&h) {
                let slot = slot as usize;
                if slot < self.sparse_fingerprints.len()
                    && self.sparse_fingerprints[slot] == sparse_fingerprint(h)
                {
                    let g = self.sparse_owners[slot];
                    *counts.entry(g).or_insert(0) += 1;
                }
            }
        }
        counts
            .into_iter()
            .filter(|&(g, c)| {
                c >= SPARSE_MIN_HITS
                    && sparse_naive_ani(c, self.genomes[g as usize].sparse_count, self.k) >= min_ani
            })
            .collect()
    }

    /// Stage-1 query: genomes passing the sparse FracMinHash ANI screen.
    pub fn hit_genomes(&self, sketch: &SequencesSketch, min_ani: f64) -> Vec<u32> {
        self.hit_genome_counts(sketch, min_ani)
            .into_iter()
            .map(|(g, _)| g)
            .collect()
    }

    /// Stage-2: the genome's full distinctive array, decoded on demand and cached.
    /// The dense block holds only the non-sparse complement; it is merged with the
    /// genome's resident stage-1 hashes to reconstruct the complete sorted array.
    fn load_genome(&self, g: u32) -> io::Result<Arc<Vec<u64>>> {
        if let Some(a) = self.cache.lock().unwrap().get(&g) {
            return Ok(a.clone());
        }
        let meta = &self.genomes[g as usize];
        let (complement, sparse) = {
            let mut rd = self.reader.lock().unwrap();
            rd.seek(SeekFrom::Start(meta.dense_offset))?;
            let mut r: &mut dyn ReadSeek = &mut **rd;
            let complement = read_hashes(&mut r)?;

            rd.seek(SeekFrom::Start(meta.sparse_offset))?;
            let mut sparse = Vec::with_capacity(meta.sparse_count as usize);
            for _ in 0..meta.sparse_count {
                let mut buf = [0u8; 8];
                rd.read_exact(&mut buf)?;
                sparse.push(u64::from_le_bytes(buf));
            }
            (complement, sparse)
        };
        // merge two sorted, disjoint ascending runs into the full distinctive array
        let mut arr = Vec::with_capacity(complement.len() + sparse.len());
        let (mut i, mut j) = (0usize, 0usize);
        while i < complement.len() && j < sparse.len() {
            if complement[i] < sparse[j] {
                arr.push(complement[i]);
                i += 1;
            } else {
                arr.push(sparse[j]);
                j += 1;
            }
        }
        arr.extend_from_slice(&complement[i..]);
        arr.extend_from_slice(&sparse[j..]);
        let arc = Arc::new(arr);
        self.cache.lock().unwrap().insert(g, arc.clone());
        Ok(arc)
    }

    /// The 2-bit nucleotide sequence of genome `g`, decoded on demand and cached.
    /// Returns `None` if no sequence is stored for that genome.
    fn load_genome_seq(&self, g: u32) -> io::Result<Option<Arc<GenomeSeq>>> {
        let off = self.genomes[g as usize].seq_offset;
        if off == 0 {
            return Ok(None);
        }
        if let Some(a) = self.seq_cache.lock().unwrap().get(&g) {
            return Ok(Some(a.clone()));
        }
        let seq = {
            let mut rd = self.reader.lock().unwrap();
            rd.seek(SeekFrom::Start(off))?;
            let mut r: &mut dyn ReadSeek = &mut **rd;
            read_genome_seq_block(&mut r)?
        };
        let arc = Arc::new(seq);
        self.seq_cache.lock().unwrap().insert(g, arc.clone());
        Ok(Some(arc))
    }

    #[inline]
    fn pool_index(&self, h: u64) -> Option<u32> {
        let slot = self.pool_mphf.try_hash(&h)? as usize;
        if slot < self.pool.len() && self.pool[slot] == h {
            Some(slot as u32)
        } else {
            None
        }
    }

    fn ensure_can_decompress(&self) -> io::Result<()> {
        if self.compression_only {
            Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "reference DB was opened in compression-only mode",
            ))
        } else {
            Ok(())
        }
    }
}

// --- adaptive present/absent subset coding ----------------------------------

/// Encode the sorted `present` indices into `domain` using the smallest of a
/// bitmask, present-Rice, or absent(complement)-Rice. Self-delimiting given the
/// `domain`, which the decoder knows from the reference DB.
fn encode_subset(out: &mut Vec<u8>, present: &[u64], domain: u64) -> io::Result<()> {
    // present-Rice (cheap: O(present))
    let mut p_rice = Vec::new();
    write_hashes(&mut p_rice, present)?;
    let bm_len = ((domain + 7) / 8) as usize;

    let mut best = (SCHEME_BITMASK, bm_len);
    if p_rice.len() < best.1 {
        best = (SCHEME_PRESENT_RICE, p_rice.len());
    }

    // absent-Rice can only beat present-Rice when there are fewer absent than
    // present indices, i.e. the subset is more than half full. For sparse subsets
    // present-Rice already dominates, so skip materializing the (up to `domain`-
    // sized) complement — otherwise a sparse sample against a huge pool/genome
    // would allocate the whole complement just to discard it.
    let mut a_rice = Vec::new();
    if present.len().saturating_mul(2) >= domain as usize {
        let mut absent = Vec::with_capacity((domain as usize).saturating_sub(present.len()));
        let mut it = present.iter().copied().peekable();
        for i in 0..domain {
            match it.peek() {
                Some(&p) if p == i => {
                    it.next();
                }
                _ => absent.push(i),
            }
        }
        write_hashes(&mut a_rice, &absent)?;
        if a_rice.len() < best.1 {
            best = (SCHEME_ABSENT_RICE, a_rice.len());
        }
    }

    let (scheme, _len) = best;
    out.push(scheme);
    match scheme {
        SCHEME_PRESENT_RICE => out.extend_from_slice(&p_rice),
        SCHEME_ABSENT_RICE => out.extend_from_slice(&a_rice),
        _ => {
            let mut bm = vec![0u8; bm_len];
            for &i in present {
                bm[(i / 8) as usize] |= 1 << (i % 8);
            }
            out.extend_from_slice(&bm);
        }
    }
    Ok(())
}

fn decode_subset<R: Read>(r: &mut R, domain: u64) -> io::Result<Vec<u64>> {
    let mut scheme = [0u8; 1];
    r.read_exact(&mut scheme)?;
    match scheme[0] {
        SCHEME_PRESENT_RICE => read_hashes(r),
        SCHEME_ABSENT_RICE => {
            let absent = read_hashes(r)?;
            let mut present = Vec::with_capacity((domain as usize).saturating_sub(absent.len()));
            let mut it = absent.iter().copied().peekable();
            for i in 0..domain {
                match it.peek() {
                    Some(&a) if a == i => {
                        it.next();
                    }
                    _ => present.push(i),
                }
            }
            Ok(present)
        }
        SCHEME_BITMASK => {
            let bm_len = ((domain + 7) / 8) as usize;
            let mut bm = vec![0u8; bm_len];
            r.read_exact(&mut bm)?;
            let mut present = Vec::new();
            for i in 0..domain {
                if bm[(i / 8) as usize] & (1 << (i % 8)) != 0 {
                    present.push(i);
                }
            }
            Ok(present)
        }
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "bad subset scheme",
        )),
    }
}

// --- compressing / decompressing a read sketch ------------------------------

/// One encoded single-substitution error k-mer: the k-mer starting at global base
/// index `pos` in a genome, with the base at offset `off` replaced by `base`.
struct ErrorEntry {
    pos: u64,
    off: u8,
    base: u8,
}

/// 2-bit reverse complement of the low `k` bases of `kmer`.
#[inline]
fn revcomp_kmer(kmer: u64, k: usize) -> u64 {
    let mut x = kmer;
    let mut rc = 0u64;
    for _ in 0..k {
        rc = (rc << 2) | (3 - (x & 3));
        x >>= 2;
    }
    rc
}

#[inline]
fn low_bits_mask(bits: u32) -> u64 {
    if bits >= 64 {
        u64::MAX
    } else {
        (1u64 << bits) - 1
    }
}

/// If `a` and `b` (each `k` 2-bit bases) differ at exactly one base, return that
/// base offset (0 = highest base); otherwise `None`.
#[inline]
fn single_base_diff(a: u64, b: u64, k: usize) -> Option<usize> {
    let d = a ^ b;
    // collapse each 2-bit group to one bit that is set iff the group is non-zero
    let nz = (d | (d >> 1)) & 0x5555_5555_5555_5555;
    if nz.count_ones() == 1 {
        let j = (nz.trailing_zeros() / 2) as usize; // group index from the low end
        Some(k - 1 - j)
    } else {
        None
    }
}

/// Try to express each `novel` hash as a single-base substitution of a hit
/// **representative** genome's k-mer (sequencing errors produce exactly such
/// k-mers). Reconstruction is deterministic and verified by hash equality, so a
/// match is never a false positive. Returns the matched entries grouped by genome
/// (already sorted by genomic position) and the set of consumed novel hashes.
///
/// Rather than enumerating all 3k substitutions of every genome k-mer, we exploit
/// two facts: (1) `mm_hash64` is a bijection, so the canonical k-mer behind a novel
/// hash is recovered exactly by `rev_mm_hash64`; (2) a single substitution leaves
/// one half of the k-mer untouched, so a genome k-mer at Hamming distance 1 shares
/// either its first or its second half exactly.
///
/// We index the (small) novel set once: each novel hash's k-mer, in both strand
/// orientations, keyed by each half (a hash map, since it is probed millions of
/// times — once per genome position — where O(1) lookups beat binary search). Then
/// each hit representative genome is scanned a single time, probing each genome
/// k-mer's two halves against the index. This amortizes one index build across all
/// hit genomes (instead of rebuilding a big per-genome index) and is O(genome)
/// scanning instead of O(genome × 3k) hashing. The genome sequence is still the
/// perfect-hash "fingerprint": each candidate is verified by recomputing the hash.
fn find_error_kmers(
    idx: &RefIndex,
    hits: &[u32],
    novel: &[u64],
    _c: usize,
    k: usize,
) -> io::Result<(Vec<(u32, Vec<ErrorEntry>)>, FxHashSet<u64>)> {
    let mut by_genome: Vec<(u32, Vec<ErrorEntry>)> = Vec::new();
    let mut consumed: FxHashSet<u64> = FxHashSet::default();
    if novel.is_empty() || k == 0 {
        return Ok((by_genome, consumed));
    }
    // each half is at most 16 bases, so a half-key fits in u32.
    let k2 = k / 2;
    let shift_hi = 2 * k2 as u32;
    let mask_full = low_bits_mask(2 * k as u32);
    let mask_low = low_bits_mask(shift_hi);

    // Index the novel k-mers (both strands) by each half: half-key -> list of
    // (novel hash, oriented k-mer). A genome k-mer that is one substitution from an
    // oriented novel k-mer shares one half exactly, so it will hit this index.
    let mut index: FxHashMap<u32, Vec<(u64, u64)>> = FxHashMap::default();
    index.reserve(novel.len() * 4);
    for &h in novel {
        let kmer = rev_mm_hash64(h) & mask_full;
        for eo in [kmer, revcomp_kmer(kmer, k)] {
            index
                .entry((eo >> shift_hi) as u32)
                .or_default()
                .push((h, eo));
            index
                .entry((eo & mask_low) as u32)
                .or_default()
                .push((h, eo));
        }
    }

    // Blocked Bloom filter: both bits for a key live in the SAME u64 word, so each
    // probe costs exactly one cache miss instead of two. The block (u64 word) is
    // chosen by the upper bits of hash1; the two bit positions within the block come
    // from the lower bits of hash1 and hash2. Capped at 64 MB so the filter always
    // fits in L3 regardless of sketch density; above the cap FPR rises but the
    // filter remains fast (L3 hits beat RAM misses from an uncapped filter).
    let n_blocks = (index.len() / 8 + 1)
        .next_power_of_two()
        .min(64 * 1024 * 1024 / 8);
    let n_block_mask = n_blocks - 1;
    let mut bloom = vec![0u64; n_blocks];
    for &key in index.keys() {
        let h1 = (key as usize).wrapping_mul(0x9e37_79b9);
        let block = (h1 >> 6) & n_block_mask;
        let bit1 = h1 & 63;
        let bit2 = (key as usize).wrapping_mul(0x517c_c1b7) & 63;
        bloom[block] |= (1u64 << bit1) | (1u64 << bit2);
    }

    let mut sorted_hits = hits.to_vec();
    sorted_hits.sort_unstable();
    let eligible: Vec<u32> = sorted_hits
        .into_iter()
        .filter(|&g| idx.genomes[g as usize].is_rep)
        .collect();

    // Phase 1 (I/O): load all eligible genome sequences into memory with the
    // producer-consumer approach so the parallel scan below has no mutex contention.
    let mut genome_seqs: Vec<(u32, Arc<GenomeSeq>)> = Vec::with_capacity(eligible.len());
    let eligible_slice: &[u32] = &eligible;
    std::thread::scope(|s| -> io::Result<()> {
        let (tx, rx) = std::sync::mpsc::sync_channel::<io::Result<(u32, Arc<GenomeSeq>)>>(8);
        s.spawn(move || {
            for &g in eligible_slice {
                match idx.load_genome_seq(g) {
                    Ok(Some(seq)) => {
                        if tx.send(Ok((g, seq))).is_err() {
                            return;
                        }
                    }
                    Ok(None) => {}
                    Err(e) => {
                        let _ = tx.send(Err(e));
                        return;
                    }
                }
            }
        });
        while let Ok(item) = rx.recv() {
            genome_seqs.push(item?);
        }
        Ok(())
    })?;

    // Phase 2 (CPU): parallel k-mer scan over all loaded genomes. `index` and `bloom`
    // are read-only so they are safe to share across Rayon threads. Each task keeps a
    // per-genome `local_consumed` set so the same novel hash is not attributed to
    // multiple positions within the same genome. Cross-genome dedup is handled by the
    // sequential pass below.
    //
    // par_iter() on a Vec is an IndexedParallelIterator, so collect() preserves
    // genome order (same as the sorted eligible order), which the dedup pass requires.
    let per_genome: Vec<(u32, Vec<(u64, ErrorEntry)>)> = genome_seqs
        .par_iter()
        .map(|(g, seq)| {
            let mut matches: Vec<(u64, ErrorEntry)> = Vec::new();
            let mut local_consumed: FxHashSet<u64> = FxHashSet::default();
            let mut base_global = 0u64;
            for contig in &seq.contigs {
                let clen = contig.len();
                if clen >= k {
                    let (mut f, _) = window_fr(contig, 0, k);
                    for start in 0..=(clen - k) {
                        if start > 0 {
                            f = ((f << 2) | contig[start + k - 1] as u64) & mask_full;
                        }
                        for key in [(f >> shift_hi) as u32, (f & mask_low) as u32] {
                            // Blocked Bloom pre-filter: one cache miss loads the block,
                            // then a pure bitmask check rejects ~97% of probes.
                            let h1 = (key as usize).wrapping_mul(0x9e37_79b9);
                            let block = (h1 >> 6) & n_block_mask;
                            let bit1 = h1 & 63;
                            let bit2 = (key as usize).wrapping_mul(0x517c_c1b7) & 63;
                            let bmask = (1u64 << bit1) | (1u64 << bit2);
                            if bloom[block] & bmask != bmask {
                                continue;
                            }
                            let Some(cands) = index.get(&key) else {
                                continue;
                            };
                            for &(h, eo) in cands {
                                if local_consumed.contains(&h) {
                                    continue;
                                }
                                if let Some(off) = single_base_diff(f, eo, k) {
                                    let j = k - 1 - off;
                                    let base = ((eo >> (2 * j)) & 3) as u8;
                                    let r = revcomp_kmer(f, k);
                                    if substituted_hash(f, r, k, off, base) == h {
                                        let global = base_global + start as u64;
                                        matches.push((
                                            h,
                                            ErrorEntry {
                                                pos: global.min(u32::MAX as u64),
                                                off: off as u8,
                                                base,
                                            },
                                        ));
                                        local_consumed.insert(h);
                                    }
                                }
                            }
                        }
                    }
                }
                base_global += clen as u64;
            }
            (*g, matches)
        })
        .collect();

    // Phase 3 (sequential): cross-genome dedup. Genomes are in sorted order; the
    // first genome to claim a novel hash wins; later genomes' entries for the same
    // hash are discarded.
    for (g, matches) in per_genome {
        let mut entries: Vec<ErrorEntry> = matches
            .into_iter()
            .filter_map(|(h, entry)| consumed.insert(h).then_some(entry))
            .collect();
        if !entries.is_empty() {
            entries.sort_unstable_by_key(|e| (e.pos, e.off, e.base));
            by_genome.push((g, entries));
        }
    }

    Ok((by_genome, consumed))
}

/// Encode the per-genome error entries as three grouped arrays (positions,
/// offsets, bases) so the outer zstd frame sees long runs of similar bytes.
/// Returns the encoded section and the total number of entries.
fn encode_error_section(by_genome: &[(u32, Vec<ErrorEntry>)]) -> io::Result<(Vec<u8>, usize)> {
    let mut out = Vec::new();
    let total: usize = by_genome.iter().map(|(_, v)| v.len()).sum();
    write_uvarint(&mut out, by_genome.len() as u64)?;
    // directory: delta-coded genome id + entry count (already ascending by id)
    let mut prev = 0u64;
    for (g, v) in by_genome {
        write_uvarint(&mut out, *g as u64 - prev)?;
        prev = *g as u64;
        write_uvarint(&mut out, v.len() as u64)?;
    }
    // array 1: per-genome Golomb-Rice position blocks (entries are pre-sorted)
    for (_, v) in by_genome {
        let positions: Vec<u64> = v.iter().map(|e| e.pos).collect();
        write_hashes(&mut out, &positions)?;
    }
    // array 2: one offset byte per entry
    for (_, v) in by_genome {
        for e in v {
            out.push(e.off);
        }
    }
    // array 3: 2-bit-packed replacement bases
    let mut byte = 0u8;
    let mut n = 0u8;
    for (_, v) in by_genome {
        for e in v {
            byte |= (e.base & 3) << (2 * n);
            n += 1;
            if n == 4 {
                out.push(byte);
                byte = 0;
                n = 0;
            }
        }
    }
    if n > 0 {
        out.push(byte);
    }
    Ok((out, total))
}

/// Compress a read sketch against the reference index. Only the sample's hit
/// genomes' dense blocks are loaded; the pool is already resident in `idx`.
pub fn compress_seq<W: Write>(
    inner: W,
    sketch: &SequencesSketch,
    idx: &RefIndex,
    ref_db_name: &str,
) -> io::Result<()> {
    compress_seq_with_screen_ani(
        inner,
        sketch,
        idx,
        ref_db_name,
        ReadSketchMeta::default(),
        REF_SCREEN_ANI_DEFAULT,
        MIN_DENSE_KMERS_FOR_ERROR_DEFAULT,
    )
}

pub fn compress_seq_with_meta<W: Write>(
    inner: W,
    sketch: &SequencesSketch,
    idx: &RefIndex,
    ref_db_name: &str,
    meta: ReadSketchMeta,
) -> io::Result<()> {
    compress_seq_with_screen_ani(
        inner,
        sketch,
        idx,
        ref_db_name,
        meta,
        REF_SCREEN_ANI_DEFAULT,
        MIN_DENSE_KMERS_FOR_ERROR_DEFAULT,
    )
}

pub fn compress_seq_with_screen_ani<W: Write>(
    inner: W,
    sketch: &SequencesSketch,
    idx: &RefIndex,
    ref_db_name: &str,
    meta: ReadSketchMeta,
    ref_screen_ani: f64,
    min_dense_kmers_for_error: usize,
) -> io::Result<()> {
    compress_seq_with_screen_ani_and_telemetry(
        inner,
        sketch,
        idx,
        ref_db_name,
        meta,
        ref_screen_ani,
        min_dense_kmers_for_error,
    )
    .map(|_| ())
}

pub fn compress_seq_with_screen_ani_and_telemetry<W: Write>(
    inner: W,
    sketch: &SequencesSketch,
    idx: &RefIndex,
    ref_db_name: &str,
    meta: ReadSketchMeta,
    ref_screen_ani: f64,
    min_dense_kmers_for_error: usize,
) -> io::Result<Vec<RefCompressTelemetry>> {
    let total_start = Instant::now();
    // stage 1 -> stage 2: load the distinctive blocks of the hit genomes and
    // build a per-sample hash -> (genome, index) lookup over just those.
    let stage1_start = Instant::now();
    let hit_counts = idx.hit_genome_counts(sketch, ref_screen_ani);
    let hits: Vec<u32> = hit_counts.iter().map(|&(g, _)| g).collect();
    info!(
        "ref-compress stage1 sparse screen: {} hit genomes in {:.3}s",
        hits.len(),
        stage1_start.elapsed().as_secs_f64()
    );

    let dense_load_start = Instant::now();
    let mut map: FxHashMap<u64, (u32, u32)> = FxHashMap::default();
    for &g in &hits {
        let arr = idx.load_genome(g)?;
        for (i, &h) in arr.iter().enumerate() {
            map.insert(h, (g, i as u32));
        }
    }
    info!(
        "ref-compress stage2 dense load + dense map build: {} genomes, {} dense hashes in {:.3}s",
        hits.len(),
        map.len(),
        dense_load_start.elapsed().as_secs_f64()
    );

    let assign_start = Instant::now();
    let mut per_genome: FxHashMap<u32, Vec<u64>> = FxHashMap::default();
    let mut pool_hits: Vec<u64> = Vec::new();
    let mut novel: Vec<u64> = Vec::new();
    for &h in sketch.kmer_counts.keys() {
        if let Some(pidx) = idx.pool_index(h) {
            pool_hits.push(pidx as u64);
        } else if let Some(&(g, i)) = map.get(&h) {
            per_genome.entry(g).or_default().push(i as u64);
        } else {
            novel.push(h);
        }
    }
    info!(
        "ref-compress stage3 sample partition: {} assigned kmers, {} pool hits, {} novel kmers in {:.3}s",
        per_genome.values().map(|v| v.len()).sum::<usize>(),
        pool_hits.len(),
        novel.len(),
        assign_start.elapsed().as_secs_f64()
    );

    // Reclassify novel hashes that are single-substitution variants of a hit
    // representative genome's k-mers into compact (position, offset, base) entries.
    let error_start = Instant::now();
    let error_by_genome = if idx.has_genome_seqs() && !novel.is_empty() {
        let error_hits: Vec<u32> = hits
            .iter()
            .copied()
            .filter(|g| {
                per_genome.get(g).map(|v| v.len()).unwrap_or(0) >= min_dense_kmers_for_error
            })
            .collect();
        let (by_genome, consumed) = find_error_kmers(idx, &error_hits, &novel, sketch.c, sketch.k)?;
        if !consumed.is_empty() {
            novel.retain(|h| !consumed.contains(h));
        }
        by_genome
    } else {
        Vec::new()
    };
    info!(
        "ref-compress stage4 error-kmer scan: {} eligible genomes, {} error genomes, {} error kmers in {:.3}s",
        hits
            .iter()
            .filter(|g| per_genome.get(g).map(|v| v.len()).unwrap_or(0) >= min_dense_kmers_for_error)
            .count(),
        error_by_genome.len(),
        error_by_genome.iter().map(|(_, v)| v.len()).sum::<usize>(),
        error_start.elapsed().as_secs_f64()
    );
    let encode_start = Instant::now();
    let (error_section, error_count) = encode_error_section(&error_by_genome)?;
    let error_counts: FxHashMap<u32, usize> =
        error_by_genome.iter().map(|(g, v)| (*g, v.len())).collect();

    // hit genomes: sorted global ids, delta-coded (strains of a species are
    // contiguous, so a sample's hits cluster into small gaps)
    let mut hit_section = Vec::new();
    let mut hit_ids: Vec<u32> = per_genome.keys().copied().collect();
    hit_ids.sort_unstable();
    write_uvarint(&mut hit_section, hit_ids.len() as u64)?;
    let mut prev = 0u64;
    let mut assigned_hits = 0usize;
    for &g in &hit_ids {
        write_uvarint(&mut hit_section, g as u64 - prev)?;
        prev = g as u64;
        let v = per_genome.get_mut(&g).unwrap();
        v.sort_unstable();
        assigned_hits += v.len();
        encode_subset(
            &mut hit_section,
            v,
            idx.genomes[g as usize].dense_domain as u64,
        )?;
    }

    // pool
    let mut pool_section = Vec::new();
    pool_hits.sort_unstable();
    write_uvarint(&mut pool_section, pool_hits.len() as u64)?;
    if !pool_hits.is_empty() {
        encode_subset(&mut pool_section, &pool_hits, idx.pool.len() as u64)?;
    }

    // novel hashes (Rice)
    let mut novel_section = Vec::new();
    write_hashes(&mut novel_section, &novel)?;

    // counts, in ascending-hash order (reproducible on decode)
    let mut count_section = Vec::new();
    let mut keys: Vec<u64> = sketch.kmer_counts.keys().copied().collect();
    keys.sort_unstable();
    for h in &keys {
        write_uvarint(&mut count_section, sketch.kmer_counts[h] as u64)?;
    }

    let mut payload = Vec::new();
    payload.extend_from_slice(SKETCH_MAGIC);
    payload.push(SKETCH_VERSION);
    payload.extend_from_slice(&idx.fingerprint.to_le_bytes());
    write_string(&mut payload, ref_db_name)?;
    write_uvarint(&mut payload, sketch.c as u64)?;
    write_uvarint(&mut payload, sketch.k as u64)?;
    write_string(&mut payload, &sketch.file_name)?;
    match &sketch.sample_name {
        Some(name) => {
            payload.push(1);
            write_string(&mut payload, name)?;
        }
        None => payload.push(0),
    }
    payload.push(sketch.paired as u8);
    payload.extend_from_slice(&sketch.mean_read_length.to_le_bytes());
    write_uvarint(&mut payload, meta.num_reads)?;
    write_uvarint(&mut payload, hit_ids.len() as u64)?;
    write_uvarint(&mut payload, assigned_hits as u64)?;
    write_uvarint(&mut payload, pool_hits.len() as u64)?;
    write_uvarint(&mut payload, novel.len() as u64)?;
    write_uvarint(&mut payload, hit_section.len() as u64)?;
    write_uvarint(&mut payload, pool_section.len() as u64)?;
    write_uvarint(&mut payload, novel_section.len() as u64)?;
    write_uvarint(&mut payload, count_section.len() as u64)?;
    // v4: single-substitution error-k-mer section
    write_uvarint(&mut payload, error_count as u64)?;
    write_uvarint(&mut payload, error_section.len() as u64)?;
    payload.extend_from_slice(&hit_section);
    payload.extend_from_slice(&pool_section);
    payload.extend_from_slice(&novel_section);
    payload.extend_from_slice(&error_section);
    payload.extend_from_slice(&count_section);

    let mut enc = zstd::stream::write::Encoder::new(inner, ZSTD_LEVEL)?;
    enc.write_all(&payload)?;
    enc.finish()?;
    info!(
        "ref-compress stage5 payload encode + zstd: payload {} bytes in {:.3}s",
        payload.len(),
        encode_start.elapsed().as_secs_f64()
    );

    let sample_name = sketch
        .sample_name
        .clone()
        .unwrap_or_else(|| sketch.file_name.clone());
    let hit_genomes_total = hits.len();
    let telemetry = hit_counts
        .into_iter()
        .map(|(g, sparse_hits)| {
            let meta = &idx.genomes[g as usize];
            RefCompressTelemetry {
                sample_name: sample_name.clone(),
                ref_db_name: ref_db_name.to_string(),
                ref_screen_ani,
                hit_genomes_total,
                genome_id: g,
                genome_file: meta.file_name.clone(),
                species: meta.species.clone(),
                sparse_hits,
                sparse_total: meta.sparse_count,
                sparse_ani: sparse_naive_ani(sparse_hits, meta.sparse_count, idx.k),
                assigned_kmers: per_genome.get(&g).map(|v| v.len()).unwrap_or(0),
                error_kmers: *error_counts.get(&g).unwrap_or(&0),
            }
        })
        .collect();
    info!(
        "ref-compress total for {}: {:.3}s",
        sample_name,
        total_start.elapsed().as_secs_f64()
    );
    Ok(telemetry)
}

/// Decompress a reference-delta read sketch. Only the genomes referenced by the
/// sample (plus the resident pool) are loaded from the index.
pub fn decompress_seq<R: Read>(inner: R, idx: &RefIndex) -> io::Result<SequencesSketch> {
    decompress_seq_with_meta(inner, idx).map(|(sketch, _)| sketch)
}

pub fn decompress_seq_with_meta<R: Read>(
    inner: R,
    idx: &RefIndex,
) -> io::Result<(SequencesSketch, ReadSketchMeta)> {
    idx.ensure_can_decompress()?;
    let raw = zstd::stream::decode_all(inner)?;
    let mut r = &raw[..];
    let mut magic = [0u8; 4];
    r.read_exact(&mut magic)?;
    if &magic != SKETCH_MAGIC {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "not a reference-delta sketch",
        ));
    }
    let mut ver = [0u8; 1];
    r.read_exact(&mut ver)?;
    if ver[0] == 0 || ver[0] > SKETCH_VERSION {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "unsupported reference-delta version",
        ));
    }
    let mut fp = [0u8; 8];
    r.read_exact(&mut fp)?;
    if u64::from_le_bytes(fp) != idx.fingerprint {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "reference DB does not match the one used to compress this sample",
        ));
    }
    if ver[0] >= 2 {
        let _ref_db_name = read_string(&mut r)?;
    }
    let c = read_uvarint(&mut r)? as usize;
    let k = read_uvarint(&mut r)? as usize;
    let file_name = read_string(&mut r)?;
    let mut tag = [0u8; 1];
    r.read_exact(&mut tag)?;
    let sample_name = if tag[0] != 0 {
        Some(read_string(&mut r)?)
    } else {
        None
    };
    let mut paired = [0u8; 1];
    r.read_exact(&mut paired)?;
    let mut mrl = [0u8; 8];
    r.read_exact(&mut mrl)?;
    let mean_read_length = f64::from_le_bytes(mrl);
    let meta = if ver[0] >= 3 {
        ReadSketchMeta {
            num_reads: read_uvarint(&mut r)?,
        }
    } else {
        ReadSketchMeta::default()
    };

    if ver[0] >= 2 {
        for _ in 0..8 {
            let _ = read_uvarint(&mut r)?;
        }
    }
    if ver[0] >= 4 {
        // error-entry count and error-section length (header metadata only)
        for _ in 0..2 {
            let _ = read_uvarint(&mut r)?;
        }
    }

    let mut hashes: Vec<u64> = Vec::new();
    let nhit = read_uvarint(&mut r)? as usize;
    let mut g = 0u64;
    for _ in 0..nhit {
        g += read_uvarint(&mut r)?;
        let domain = idx.genomes[g as usize].dense_domain as u64;
        let indices = decode_subset(&mut r, domain)?;
        let arr = idx.load_genome(g as u32)?;
        for i in indices {
            hashes.push(arr[i as usize]);
        }
    }
    let npool = read_uvarint(&mut r)? as usize;
    if npool > 0 {
        let indices = decode_subset(&mut r, idx.pool.len() as u64)?;
        for i in indices {
            hashes.push(idx.pool[i as usize]);
        }
    }
    let novel = read_hashes(&mut r)?;
    hashes.extend_from_slice(&novel);

    // error k-mers: reconstruct each (genome, position, offset, base) entry back
    // to its canonical FracMinHash hash using the stored genome sequence.
    if ver[0] >= 4 {
        let n_err_genomes = read_uvarint(&mut r)? as usize;
        let mut gids = Vec::with_capacity(n_err_genomes);
        let mut counts_e = Vec::with_capacity(n_err_genomes);
        let mut gg = 0u64;
        for _ in 0..n_err_genomes {
            gg += read_uvarint(&mut r)?;
            let cnt = read_uvarint(&mut r)? as usize;
            gids.push(gg as u32);
            counts_e.push(cnt);
        }
        // array 1: per-genome Golomb-Rice position blocks
        let mut positions: Vec<Vec<u64>> = Vec::with_capacity(n_err_genomes);
        for _ in 0..n_err_genomes {
            positions.push(read_hashes(&mut r)?);
        }
        let total: usize = counts_e.iter().sum();
        // array 2: one offset byte per entry
        let mut offsets = vec![0u8; total];
        r.read_exact(&mut offsets)?;
        // array 3: 2-bit-packed replacement bases
        let mut packed = vec![0u8; (total + 3) / 4];
        r.read_exact(&mut packed)?;

        let mut e = 0usize;
        for gi in 0..n_err_genomes {
            let seq = idx.load_genome_seq(gids[gi])?.ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "error k-mers reference a genome with no stored sequence",
                )
            })?;
            let cum = seq.cumulative();
            for &pos in &positions[gi] {
                let off = offsets[e] as usize;
                let base = (packed[e / 4] >> (2 * (e % 4))) & 3;
                let (ci, local) = seq.locate(&cum, pos);
                let (f, rr) = window_fr(&seq.contigs[ci], local, k);
                hashes.push(substituted_hash(f, rr, k, off, base));
                e += 1;
            }
        }
    }

    hashes.sort_unstable();
    let mut kmer_counts = FxHashMap::with_capacity_and_hasher(hashes.len(), Default::default());
    for &h in &hashes {
        let count = read_uvarint(&mut r)? as u32;
        kmer_counts.insert(h, count);
    }

    Ok((
        SequencesSketch {
            kmer_counts,
            c,
            k,
            file_name,
            sample_name,
            paired: paired[0] != 0,
            mean_read_length,
        },
        meta,
    ))
}

// --- CLI handlers -----------------------------------------------------------

fn init_logger(trace: bool) {
    let level = if trace {
        log::LevelFilter::Trace
    } else {
        log::LevelFilter::Info
    };
    let _ = simple_logger::SimpleLogger::new().with_level(level).init();
}

fn load_seq_sketch(path: &str) -> SequencesSketch {
    let file = File::open(path).unwrap_or_else(|_| panic!("Could not open {}", path));
    let mut reader = BufReader::with_capacity(10_000_000, file);
    if compress::peek_is_compressed(&mut reader).unwrap_or(false) {
        read_seq_sketch_compressed(&mut reader)
            .unwrap_or_else(|_| panic!("{} is not a valid sample sketch", path))
    } else {
        bincode::deserialize_from(&mut reader)
            .unwrap_or_else(|_| panic!("{} is not a valid sample sketch", path))
    }
}

fn compare_seq_sketches(
    original: &SequencesSketch,
    decoded: &SequencesSketch,
) -> Result<(), String> {
    if original.c != decoded.c {
        return Err(format!(
            "c differs: original={} decoded={}",
            original.c, decoded.c
        ));
    }
    if original.k != decoded.k {
        return Err(format!(
            "k differs: original={} decoded={}",
            original.k, decoded.k
        ));
    }
    if original.file_name != decoded.file_name {
        return Err(format!(
            "file_name differs: original={:?} decoded={:?}",
            original.file_name, decoded.file_name
        ));
    }
    if original.sample_name != decoded.sample_name {
        return Err(format!(
            "sample_name differs: original={:?} decoded={:?}",
            original.sample_name, decoded.sample_name
        ));
    }
    if original.paired != decoded.paired {
        return Err(format!(
            "paired differs: original={} decoded={}",
            original.paired, decoded.paired
        ));
    }
    if original.mean_read_length != decoded.mean_read_length {
        return Err(format!(
            "mean_read_length differs: original={} decoded={}",
            original.mean_read_length, decoded.mean_read_length
        ));
    }
    if original.kmer_counts.len() != decoded.kmer_counts.len() {
        return Err(format!(
            "hash count differs: original={} decoded={}",
            original.kmer_counts.len(),
            decoded.kmer_counts.len()
        ));
    }
    for (&hash, &count) in &original.kmer_counts {
        match decoded.kmer_counts.get(&hash) {
            Some(&decoded_count) if decoded_count == count => {}
            Some(&decoded_count) => {
                return Err(format!(
                    "count differs for hash {}: original={} decoded={}",
                    hash, count, decoded_count
                ));
            }
            None => return Err(format!("decoded sketch is missing hash {}", hash)),
        }
    }
    if let Some((&hash, _)) = decoded
        .kmer_counts
        .iter()
        .find(|(hash, _)| !original.kmer_counts.contains_key(hash))
    {
        return Err(format!("decoded sketch has extra hash {}", hash));
    }
    Ok(())
}

fn verify_ref_sketch(path: &Path, original: &SequencesSketch, idx: &RefIndex) -> io::Result<()> {
    let r = BufReader::with_capacity(
        10_000_000,
        File::open(path).unwrap_or_else(|_| panic!("Could not open {:?}", path)),
    );
    let decoded = decompress_seq(r, idx)?;
    compare_seq_sketches(original, &decoded)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn original_path_candidates(ref_path: &Path, sample_file: &str) -> Vec<PathBuf> {
    let mut candidates = Vec::new();

    if let Some(stem) = ref_path
        .file_name()
        .and_then(|name| name.to_str())
        .and_then(|name| name.strip_suffix(REF_SAMPLE_SUFFIX))
    {
        if let Some(parent) = ref_path.parent() {
            candidates.push(parent.join(format!("{}{}", stem, SAMPLE_FILE_SUFFIX)));
            candidates.push(parent.join(format!("{}{}", stem, SAMPLE_COMP_FILE_SUFFIX)));
        }
    }

    if sample_file.ends_with(SAMPLE_FILE_SUFFIX) || sample_file.ends_with(SAMPLE_COMP_FILE_SUFFIX) {
        candidates.push(PathBuf::from(sample_file));
        if let Some(name) = Path::new(sample_file).file_name() {
            if let Some(parent) = ref_path.parent() {
                candidates.push(parent.join(name));
            }
        }
    }

    candidates
}

fn find_original_sketch_path(ref_path: &Path, sample_file: &str) -> io::Result<PathBuf> {
    for candidate in original_path_candidates(ref_path, sample_file) {
        if candidate.exists() {
            return Ok(candidate);
        }
    }
    Err(io::Error::new(
        io::ErrorKind::NotFound,
        format!(
            "could not find original sketch for {:?}; tried stored sample path {:?} and same-directory .sylsp/.sylspc fallbacks",
            ref_path, sample_file
        ),
    ))
}

#[derive(Debug)]
struct RefSketchInspect {
    path: String,
    compressed_bytes: u64,
    payload_bytes: usize,
    version: u8,
    reference_fingerprint: u64,
    reference_db: String,
    c: usize,
    k: usize,
    sample_file: String,
    sample_name: String,
    paired: bool,
    mean_read_length: f64,
    num_reads: Option<u64>,
    header_metadata_bytes: usize,
    hit_genomes: Option<u64>,
    assigned_to_genomes: Option<u64>,
    shared_pool: Option<u64>,
    novel: Option<u64>,
    hit_section_bytes: Option<u64>,
    pool_section_bytes: Option<u64>,
    novel_section_bytes: Option<u64>,
    count_section_bytes: Option<u64>,
    error_kmers: Option<u64>,
    error_section_bytes: Option<u64>,
}

fn inspect_ref_sketch(path: &str) -> io::Result<RefSketchInspect> {
    let compressed_bytes = std::fs::metadata(path)?.len();
    let raw = zstd::stream::decode_all(BufReader::new(File::open(path)?))?;
    let mut r = &raw[..];
    let mut magic = [0u8; 4];
    r.read_exact(&mut magic)?;
    if &magic != SKETCH_MAGIC {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "not a reference-delta sketch",
        ));
    }
    let mut ver = [0u8; 1];
    r.read_exact(&mut ver)?;
    if ver[0] == 0 || ver[0] > SKETCH_VERSION {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "unsupported reference-delta version",
        ));
    }
    let mut fp = [0u8; 8];
    r.read_exact(&mut fp)?;
    let reference_fingerprint = u64::from_le_bytes(fp);
    let reference_db = if ver[0] >= 2 {
        read_string(&mut r)?
    } else {
        "unknown".to_string()
    };
    let c = read_uvarint(&mut r)? as usize;
    let k = read_uvarint(&mut r)? as usize;
    let sample_file = read_string(&mut r)?;
    let mut tag = [0u8; 1];
    r.read_exact(&mut tag)?;
    let sample_name = if tag[0] != 0 {
        read_string(&mut r)?
    } else {
        String::new()
    };
    let mut paired = [0u8; 1];
    r.read_exact(&mut paired)?;
    let mut mrl = [0u8; 8];
    r.read_exact(&mut mrl)?;
    let mean_read_length = f64::from_le_bytes(mrl);
    let num_reads = if ver[0] >= 3 {
        Some(read_uvarint(&mut r)?)
    } else {
        None
    };

    let (
        hit_genomes,
        assigned_to_genomes,
        shared_pool,
        novel,
        hit_section_bytes,
        pool_section_bytes,
        novel_section_bytes,
        count_section_bytes,
    ) = if ver[0] >= 2 {
        (
            Some(read_uvarint(&mut r)?),
            Some(read_uvarint(&mut r)?),
            Some(read_uvarint(&mut r)?),
            Some(read_uvarint(&mut r)?),
            Some(read_uvarint(&mut r)?),
            Some(read_uvarint(&mut r)?),
            Some(read_uvarint(&mut r)?),
            Some(read_uvarint(&mut r)?),
        )
    } else {
        (None, None, None, None, None, None, None, None)
    };
    let (error_kmers, error_section_bytes) = if ver[0] >= 4 {
        (Some(read_uvarint(&mut r)?), Some(read_uvarint(&mut r)?))
    } else {
        (None, None)
    };
    let header_metadata_bytes = raw.len() - r.len();

    Ok(RefSketchInspect {
        path: path.to_string(),
        compressed_bytes,
        payload_bytes: raw.len(),
        version: ver[0],
        reference_fingerprint,
        reference_db,
        c,
        k,
        sample_file,
        sample_name,
        paired: paired[0] != 0,
        mean_read_length,
        num_reads,
        header_metadata_bytes,
        hit_genomes,
        assigned_to_genomes,
        shared_pool,
        novel,
        hit_section_bytes,
        pool_section_bytes,
        novel_section_bytes,
        count_section_bytes,
        error_kmers,
        error_section_bytes,
    })
}

fn opt_u64(x: Option<u64>) -> String {
    x.map(|v| v.to_string()).unwrap_or_else(|| "-".to_string())
}

fn run_ref_inspect(files: &[String]) {
    println!(
        "file\tversion\tcompressed_bytes\tpayload_bytes\treference_fingerprint\treference_db\tsample_file\tsample_name\tc\tk\tpaired\tmean_read_length\tnum_reads\theader_metadata_payload_bytes\thit_genomes\tassigned_to_genomes\tshared_pool\tnovel\terror_kmers\ttotal_hashes\thit_section_payload_bytes\tpool_section_payload_bytes\tnovel_section_payload_bytes\terror_section_payload_bytes\tcount_section_payload_bytes"
    );
    for f in files {
        let x = inspect_ref_sketch(f).unwrap_or_else(|e| panic!("Failed to inspect {}: {}", f, e));
        let err = x.error_kmers.unwrap_or(0);
        let total_hashes = match (x.assigned_to_genomes, x.shared_pool, x.novel) {
            (Some(a), Some(p), Some(n)) => (a + p + n + err).to_string(),
            _ => "-".to_string(),
        };
        println!(
            "{}\t{}\t{}\t{}\t{:016x}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            x.path,
            x.version,
            x.compressed_bytes,
            x.payload_bytes,
            x.reference_fingerprint,
            x.reference_db,
            x.sample_file,
            x.sample_name,
            x.c,
            x.k,
            x.paired,
            x.mean_read_length,
            opt_u64(x.num_reads),
            x.header_metadata_bytes,
            opt_u64(x.hit_genomes),
            opt_u64(x.assigned_to_genomes),
            opt_u64(x.shared_pool),
            opt_u64(x.novel),
            opt_u64(x.error_kmers),
            total_hashes,
            opt_u64(x.hit_section_bytes),
            opt_u64(x.pool_section_bytes),
            opt_u64(x.novel_section_bytes),
            opt_u64(x.error_section_bytes),
            opt_u64(x.count_section_bytes),
        );
    }
}

fn run_ref_verify(files: &[String], idx: &RefIndex) {
    for f in files {
        let ref_path = Path::new(f);
        let inspect =
            inspect_ref_sketch(f).unwrap_or_else(|e| panic!("Failed to inspect {}: {}", f, e));
        let original_path = find_original_sketch_path(ref_path, &inspect.sample_file)
            .unwrap_or_else(|e| panic!("Failed to locate original sketch for {}: {}", f, e));
        let original = load_seq_sketch(original_path.to_str().unwrap_or_else(|| {
            panic!(
                "Original sketch path is not valid UTF-8: {:?}",
                original_path
            )
        }));
        verify_ref_sketch(ref_path, &original, idx).unwrap_or_else(|e| {
            panic!(
                "Verification failed for {} against {:?}: {}",
                f, original_path, e
            )
        });
        info!("Verified {} against {:?}", f, original_path);
    }
}

fn parse_taxonomy(path: &str) -> FxHashMap<String, (String, bool)> {
    let text = std::fs::read_to_string(path)
        .unwrap_or_else(|_| panic!("Could not read taxonomy file {}", path));
    let mut map = FxHashMap::default();
    for (lineno, line) in text.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 3 {
            error!(
                "Taxonomy line {} has fewer than 3 tab-separated columns; exiting",
                lineno + 1
            );
            std::process::exit(1);
        }
        let is_rep = match cols[2].trim().to_ascii_lowercase().as_str() {
            "rep" | "representative" | "1" => true,
            "strain" | "2" => false,
            other => {
                error!(
                    "Taxonomy line {}: level must be rep/strain (got '{}'); exiting",
                    lineno + 1,
                    other
                );
                std::process::exit(1);
            }
        };
        map.insert(
            cols[0].trim().to_string(),
            (cols[1].trim().to_string(), is_rep),
        );
    }
    map
}

// --- streaming, hash-partitioned, parallel build ---------------------------

/// Map a hash to one of `p` partitions by its low bits. FracMinHash hashes all
/// sit below the MinHash threshold (~2^64/c), so high-bit partitioning would pile
/// everything into the lowest partitions; the low bits are uniform. This is not
/// order-preserving, so each genome's hashes are sorted once after merging.
#[inline]
fn partition_of(h: u64, p: u64) -> usize {
    (h % p) as usize
}

/// Read a LEB128 uvarint, returning `None` at a clean end of stream (used to
/// iterate the variable number of genome blocks in a shard file).
fn read_uvarint_opt<R: Read>(r: &mut R) -> io::Result<Option<u64>> {
    let mut first = [0u8; 1];
    match r.read(&mut first)? {
        0 => return Ok(None),
        _ => {}
    }
    let mut x = (first[0] & 0x7f) as u64;
    if first[0] & 0x80 == 0 {
        return Ok(Some(x));
    }
    let mut shift = 7u32;
    loop {
        let mut b = [0u8; 1];
        r.read_exact(&mut b)?;
        x |= ((b[0] & 0x7f) as u64) << shift;
        if b[0] & 0x80 == 0 {
            break;
        }
        shift += 7;
    }
    Ok(Some(x))
}

/// Stream every genome sketch from the input databases, one at a time. Compressed
/// databases are streamed without materializing all sketches; legacy bincode
/// databases are loaded per-file (they cannot be streamed).
fn for_each_genome<F: FnMut(GenomeSketch) -> io::Result<()>>(files: &[String], mut f: F) {
    for path in files {
        info!("Streaming genome sketches from {}", path);
        let file = File::open(path).unwrap_or_else(|_| panic!("Could not open {}", path));
        let mut reader = BufReader::with_capacity(10_000_000, file);
        if compress::peek_is_compressed(&mut reader).unwrap_or(false) {
            compress::stream_genome_sketches_compressed(reader, &mut f)
                .unwrap_or_else(|e| panic!("{} is not a valid database sketch: {}", path, e));
        } else {
            warn!(
                "{} is a legacy (uncompressed) database; it must be loaded in full. Re-sketch with --compressed-output for low-RAM builds.",
                path
            );
            let sketches: Vec<GenomeSketch> = bincode::deserialize_from(&mut reader)
                .unwrap_or_else(|_| panic!("{} is not a valid database sketch", path));
            for s in sketches {
                f(s).unwrap();
            }
        }
    }
}

#[derive(Default)]
struct ShardOut {
    /// distinctive hashes by genome id (unsorted; the merge step sorts them).
    dist: Vec<(u32, Vec<u64>)>,
    /// pool hashes owned by this shard (unsorted; the merge step sorts them).
    pool: Vec<u64>,
}

/// Build ownership for one shard (one hash class): read its genome blocks, tally
/// rep/strain occurrences per hash, and assign each hash to a single owner genome
/// or the shared pool, exactly as `build_refdb` does globally.
fn process_shard(path: &Path, remap: &[u32], pool_min_genomes: u32) -> io::Result<ShardOut> {
    let mut r = BufReader::with_capacity(1 << 20, File::open(path)?);
    let mut acc: FxHashMap<u64, OwnAccum> = FxHashMap::default();
    while let Some(fid) = read_uvarint_opt(&mut r)? {
        let mut rep = [0u8; 1];
        r.read_exact(&mut rep)?;
        let is_rep = rep[0] != 0;
        let gid = remap[fid as usize];
        let hashes = read_hashes(&mut r)?;
        for h in hashes {
            let e = acc.entry(h).or_insert(OwnAccum {
                rep_count: 0,
                rep_id: 0,
                strain_count: 0,
                strain_id: 0,
            });
            if is_rep {
                if e.rep_count == 0 {
                    e.rep_id = gid;
                }
                e.rep_count += 1;
            } else {
                if e.strain_count == 0 {
                    e.strain_id = gid;
                }
                e.strain_count += 1;
            }
        }
    }

    let mut by_gid: FxHashMap<u32, Vec<u64>> = FxHashMap::default();
    let mut pool: Vec<u64> = Vec::new();
    for (h, a) in acc {
        let owner = owner_for_accum(&a, pool_min_genomes);
        if owner == POOL {
            pool.push(h);
        } else {
            by_gid.entry(owner).or_default().push(h);
        }
    }
    // Not sorted here: modular partitioning isn't order-preserving, so the merge
    // step sorts each genome's hashes and the pool once after concatenation.
    let dist: Vec<(u32, Vec<u64>)> = by_gid.into_iter().collect();
    Ok(ShardOut { dist, pool })
}

/// Mark every species-representative genome for sequence storage, pointing at the
/// source FASTA the sketch recorded. Strains are skipped (errors are only encoded
/// against representatives). The FASTAs are not read here — `write_refdb` streams
/// and 2-bit packs them one genome at a time, so even a 200k-genome reference does
/// not hold every sequence in RAM. A FASTA that turns out to be unreadable simply
/// yields no stored sequence (that genome falls back to full-price coding).
fn set_rep_seq_sources(db: &mut RefDb) {
    let sources: Vec<Option<GenomeSource>> = db
        .genomes
        .iter()
        .map(|g| {
            if g.is_rep {
                Some(GenomeSource::Path(g.file_name.clone()))
            } else {
                None
            }
        })
        .collect();
    let reps = sources.iter().filter(|s| s.is_some()).count();
    info!(
        "Will store nucleotide sequences for {} representative genome(s) (streamed + 2-bit packed at write time)",
        reps
    );
    db.rep_seqs = sources;
}

pub fn run_ref_build(args: RefBuildArgs) {
    init_logger(args.trace);
    let threads = args.threads.max(1);
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .ok();
    if args.files.is_empty() {
        error!("No genome database sketches (*.syldb) supplied; exiting");
        std::process::exit(1);
    }
    let taxonomy = match &args.taxonomy {
        Some(p) => parse_taxonomy(p),
        None => FxHashMap::default(),
    };
    if args.sparse_div_compat.is_some() {
        warn!("--sparse-subsample is deprecated and ignored; use --sparse-c to size the sparse MPHF index");
    }
    let sparse_c = args.sparse_c.max(1);
    let pool_min_genomes = args.pool_min_genomes.max(2);

    let resolve = |file_name: &str| -> (String, bool) {
        if let Some(v) = taxonomy.get(file_name) {
            return v.clone();
        }
        let base = Path::new(file_name)
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or(file_name);
        if let Some(v) = taxonomy.get(base) {
            return v.clone();
        }
        (file_name.to_string(), true)
    };

    // Choose the partition count from the RAM target (soft) and the input size.
    // Distinct k-mers are bounded by k-mer instances, which we estimate from the
    // compressed input size; per-shard `acc` peak is ~ (distinct / p) * entry size.
    let total_input_bytes: u64 = args
        .files
        .iter()
        .filter_map(|f| std::fs::metadata(f).ok())
        .map(|m| m.len())
        .sum();
    const EST_BYTES_PER_KMER: f64 = 2.0;
    const ACC_BYTES_PER_ENTRY: f64 = 48.0;
    let est_kmers = (total_input_bytes as f64 / EST_BYTES_PER_KMER).max(1.0);
    let p = match args.max_ram {
        Some(gb) if gb > 0 => {
            // reserve ~half the budget for `acc`, the rest for output + buffers
            let budget = (gb as f64) * 1e9 * 0.5;
            let needed = est_kmers * ACC_BYTES_PER_ENTRY * threads as f64 / budget;
            (needed.ceil() as usize).max(threads).min(512).max(1)
        }
        _ => 256usize.min(512).max(threads).max(1),
    };
    info!(
        "Building reference with {} hash partitions ({} threads), shared-pool threshold >= {} genome(s){}",
        p,
        threads,
        pool_min_genomes,
        match args.max_ram {
            Some(gb) => format!(", ~{} GB RAM target", gb),
            None => String::new(),
        }
    );

    let out = format!("{}{}", args.output, REF_DB_SUFFIX);
    if let Some(parent) = Path::new(&out).parent() {
        std::fs::create_dir_all(parent).ok();
    }
    let tmp_dir = match &args.tmp_dir {
        Some(d) => PathBuf::from(d),
        None => Path::new(&out)
            .parent()
            .unwrap_or_else(|| Path::new("."))
            .to_path_buf(),
    }
    .join(format!(".sylref_build_{}", std::process::id()));
    std::fs::create_dir_all(&tmp_dir)
        .unwrap_or_else(|e| panic!("Could not create scratch dir {:?}: {}", tmp_dir, e));
    let shard_path = |i: usize| tmp_dir.join(format!("shard_{}.bin", i));

    // --- pass 1: stream genomes, route k-mers into hash-partitioned shards ---
    let mut writers: Vec<BufWriter<File>> = (0..p)
        .map(|i| {
            BufWriter::with_capacity(
                1 << 16,
                File::create(shard_path(i))
                    .unwrap_or_else(|e| panic!("Could not create scratch shard: {}", e)),
            )
        })
        .collect();
    let mut meta: Vec<(String, String, bool)> = Vec::new();
    let mut c = 0usize;
    let mut k = 0usize;
    let mut total_instances: u64 = 0;
    let pu = p as u64;
    {
        let mut fid = 0u32;
        // reused per genome: one bucket of k-mers per partition
        let mut buckets: Vec<Vec<u64>> = vec![Vec::new(); p];
        for_each_genome(&args.files, |s| {
            if meta.is_empty() {
                c = s.c;
                k = s.k;
            }
            let (species, is_rep) = resolve(&s.file_name);
            meta.push((s.file_name.clone(), species, is_rep));
            let rep_byte = is_rep as u8;
            let mut kmers = s.genome_kmers;
            kmers.sort_unstable();
            kmers.dedup();
            total_instances += kmers.len() as u64;
            // distribute this genome's (sorted) k-mers into partition buckets; each
            // bucket stays sorted because we scan in order.
            for b in buckets.iter_mut() {
                b.clear();
            }
            for &h in &kmers {
                buckets[partition_of(h, pu)].push(h);
            }
            for (part, b) in buckets.iter().enumerate() {
                if b.is_empty() {
                    continue;
                }
                let w = &mut writers[part];
                write_uvarint(w, fid as u64)?;
                w.write_all(&[rep_byte])?;
                write_hashes(w, b)?;
            }
            fid += 1;
            Ok(())
        });
    }
    for mut w in writers {
        w.flush()
            .unwrap_or_else(|e| panic!("Failed to flush scratch shard: {}", e));
    }

    let ng = meta.len();
    if ng == 0 {
        std::fs::remove_dir_all(&tmp_dir).ok();
        error!("No genome sketches found; exiting");
        std::process::exit(1);
    }
    info!(
        "Routed {} genomes ({} k-mers) into {} partitions; assigning owners...",
        ng, total_instances, p
    );

    // build order (species, reps first, then file name) and file-id -> genome-id remap
    let mut order: Vec<usize> = (0..ng).collect();
    order.sort_by(|&a, &b| {
        meta[a]
            .1
            .cmp(&meta[b].1)
            .then(meta[b].2.cmp(&meta[a].2))
            .then(meta[a].0.cmp(&meta[b].0))
    });
    let mut remap = vec![0u32; ng];
    let mut genomes: Vec<RefGenome> = Vec::with_capacity(ng);
    for (gid, &fid) in order.iter().enumerate() {
        remap[fid] = gid as u32;
        genomes.push(RefGenome {
            file_name: meta[fid].0.clone(),
            species: meta[fid].1.clone(),
            is_rep: meta[fid].2,
        });
    }

    // --- pass 2: build ownership per shard in parallel ----------------------
    let shard_outs: Vec<ShardOut> = (0..p)
        .into_par_iter()
        .map(|pi| {
            process_shard(&shard_path(pi), &remap, pool_min_genomes)
                .unwrap_or_else(|e| panic!("Failed to process scratch shard {}: {}", pi, e))
        })
        .collect();
    std::fs::remove_dir_all(&tmp_dir).ok();

    // merge shards (any order), then sort each genome's hashes and the pool, since
    // modular partitioning does not preserve hash order.
    let mut distinctive: Vec<Vec<u64>> = vec![Vec::new(); ng];
    let mut pool: Vec<u64> = Vec::new();
    for mut shard in shard_outs {
        for (gid, mut hashes) in shard.dist.drain(..) {
            distinctive[gid as usize].append(&mut hashes);
        }
        pool.append(&mut shard.pool);
    }
    distinctive.par_iter_mut().for_each(|d| d.sort_unstable());
    pool.par_sort_unstable();

    let mut db = RefDb {
        c,
        k,
        genomes,
        distinctive,
        pool,
        rep_seqs: Vec::new(),
    };
    if args.store_genomes {
        set_rep_seq_sources(&mut db);
    }
    let n_dist: usize = db.distinctive.iter().map(|d| d.len()).sum();
    let effective_sparse_c = sparse_c.max(db.c).max(1);
    let n_sparse: usize = db
        .distinctive
        .iter()
        .map(|d| {
            d.iter()
                .filter(|&&h| keep_sparse_hash(h, effective_sparse_c))
                .count()
        })
        .sum();
    info!(
        "Reference: {} genomes, {} distinctive k-mers, {} shared-pool k-mers, {} stage-1 sparse k-mers (sparse c={})",
        db.genomes.len(),
        n_dist,
        db.pool.len(),
        n_sparse,
        effective_sparse_c
    );

    let w =
        BufWriter::new(File::create(&out).unwrap_or_else(|_| panic!("Could not create {}", out)));
    write_refdb(w, &db, sparse_c).unwrap_or_else(|e| panic!("Failed to write {}: {}", out, e));
    info!("Wrote reference database to {}", out);
}

/// Open a `.sylref` for querying/compression.
fn open_refdb_file(path: &str) -> RefIndex {
    open_refdb_file_with_mode(path, false)
}

fn open_refdb_file_for_compress(path: &str) -> RefIndex {
    open_refdb_file_with_mode(path, true)
}

fn open_refdb_file_with_mode(path: &str, compression_only: bool) -> RefIndex {
    let r = BufReader::with_capacity(
        10_000_000,
        File::open(path).unwrap_or_else(|_| panic!("Could not open {}", path)),
    );
    if compression_only {
        open_ref_index_for_compress(r)
            .unwrap_or_else(|e| panic!("{} is not a valid reference DB: {}", path, e))
    } else {
        open_ref_index(r).unwrap_or_else(|e| panic!("{} is not a valid reference DB: {}", path, e))
    }
}

fn telemetry_field(s: &str) -> String {
    s.replace(['\t', '\n', '\r'], " ")
}

fn write_telemetry_header<W: Write>(w: &mut W) -> io::Result<()> {
    writeln!(
        w,
        "sample\tref_db\tref_screen_ani\thit_genomes_total\tgenome_id\tgenome_file\tspecies\tsparse_hits\tsparse_total\tsparse_ani\tassigned_kmers\terror_kmers"
    )
}

fn write_telemetry_rows<W: Write>(w: &mut W, rows: &[RefCompressTelemetry]) -> io::Result<()> {
    for row in rows {
        writeln!(
            w,
            "{}\t{}\t{:.3}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{}\t{}",
            telemetry_field(&row.sample_name),
            telemetry_field(&row.ref_db_name),
            row.ref_screen_ani,
            row.hit_genomes_total,
            row.genome_id,
            telemetry_field(&row.genome_file),
            telemetry_field(&row.species),
            row.sparse_hits,
            row.sparse_total,
            row.sparse_ani,
            row.assigned_kmers,
            row.error_kmers
        )?;
    }
    Ok(())
}

pub fn run_ref_compress(args: RefCompressArgs) {
    init_logger(args.trace);
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .ok();
    if args.files.is_empty() {
        error!("No sample sketches supplied; exiting");
        std::process::exit(1);
    }
    if args.inspect {
        run_ref_inspect(&args.files);
        return;
    }
    let ref_db = match &args.ref_db {
        Some(r) => r,
        None => {
            error!("--reference is required unless --inspect is used; exiting");
            std::process::exit(1);
        }
    };
    if args.verify && args.decompress {
        error!("--verify expects existing *.sylspr inputs and cannot be combined with --decompress; exiting");
        std::process::exit(1);
    }
    if args.telemetry.is_some() && (args.decompress || args.verify) {
        error!("--telemetry is only supported while compressing; exiting");
        std::process::exit(1);
    }
    let idx = if args.decompress || args.verify {
        info!("Loading reference database {} (full mode)", ref_db);
        open_refdb_file(ref_db)
    } else {
        info!(
            "Loading reference database {} (compression-only mode: stage-1 sparse index + MPHF pool)",
            ref_db
        );
        open_refdb_file_for_compress(ref_db)
    };

    if args.verify {
        run_ref_verify(&args.files, &idx);
        info!("Done ({} sample(s)).", args.files.len());
        return;
    }

    std::fs::create_dir_all(&args.output_dir).expect("Could not create output directory; exiting");

    let outdir = Path::new(&args.output_dir);
    let counter = Mutex::new(0usize);
    let telemetry_writer = args.telemetry.as_ref().map(|path| {
        let mut w = BufWriter::new(
            File::create(path)
                .unwrap_or_else(|_| panic!("Could not create telemetry file {}", path)),
        );
        write_telemetry_header(&mut w)
            .unwrap_or_else(|e| panic!("Could not write telemetry header to {}: {}", path, e));
        Mutex::new(w)
    });
    if args.decompress {
        args.files.par_iter().for_each(|f| {
            let r = BufReader::with_capacity(
                10_000_000,
                File::open(f).unwrap_or_else(|_| panic!("Could not open {}", f)),
            );
            let sketch = decompress_seq(r, &idx).unwrap_or_else(|e| {
                error!("Failed to decompress {}: {}", f, e);
                std::process::exit(1);
            });
            let base = Path::new(f).file_name().unwrap().to_str().unwrap();
            let stem = base.strip_suffix(REF_SAMPLE_SUFFIX).unwrap_or(base);
            let out = outdir.join(format!("{}{}", stem, SAMPLE_FILE_SUFFIX));
            let mut w = BufWriter::new(
                File::create(&out).unwrap_or_else(|_| panic!("Could not create {:?}", out)),
            );
            bincode::serialize_into(&mut w, &sketch).unwrap();
            let mut c = counter.lock().unwrap();
            *c += 1;
            info!("Decompressed {} -> {:?}", f, out);
        });
    } else {
        args.files.par_iter().for_each(|f| {
            let sketch = load_seq_sketch(f);
            let base = Path::new(f).file_name().unwrap().to_str().unwrap();
            let stem = base
                .strip_suffix(SAMPLE_FILE_SUFFIX)
                .or_else(|| base.strip_suffix(SAMPLE_COMP_FILE_SUFFIX))
                .unwrap_or(base);
            let out = outdir.join(format!("{}{}", stem, REF_SAMPLE_SUFFIX));
            let w = BufWriter::new(
                File::create(&out).unwrap_or_else(|_| panic!("Could not create {:?}", out)),
            );
            let telemetry = compress_seq_with_screen_ani_and_telemetry(
                w,
                &sketch,
                &idx,
                ref_db,
                ReadSketchMeta::default(),
                args.ref_screen_ani,
                args.min_dense_kmers_for_error,
            )
            .unwrap_or_else(|e| panic!("Failed to compress {}: {}", f, e));
            if let Some(writer) = &telemetry_writer {
                let mut writer = writer.lock().unwrap();
                write_telemetry_rows(&mut *writer, &telemetry)
                    .unwrap_or_else(|e| panic!("Failed to write telemetry for {}: {}", f, e));
            }
            let mut c = counter.lock().unwrap();
            *c += 1;
            info!("Compressed {} -> {:?}", f, out);
        });
    }
    info!("Done ({} sample(s)).", *counter.lock().unwrap());
}
