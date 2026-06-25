//! Building, writing, and reading the seekable `.sylref` reference database.

use crate::cmdline::RefBuildArgs;
use crate::compress::{
    self, read_hashes, read_string, read_uvarint, write_hashes, write_string, write_uvarint,
};
use crate::constants::*;
use crate::types::*;
use boomphf::Mphf;
use fxhash::FxHashMap;
use log::*;
use rayon::prelude::*;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

const REFDB_MAGIC: &[u8; 4] = b"SYLR";
/// v5 adds an optional per-rep-genome 2-bit nucleotide section (`--store-genomes`)
/// used by error-k-mer encoding. v4 files (no genomes) are still readable.
const REFDB_VERSION: u8 = 5;
const REFDB_MIN_VERSION: u8 = 4;

/// Fixed size of the `.sylref` header: magic (4) + version (1) + footer offset (8).
const HEADER_LEN: u64 = 13;

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

/// A genome's nucleotide sequence in 2-bit-packed form (4 bases per byte,
/// same layout as the on-disk format). Each contig is `(length_in_bases,
/// packed_data)`. K-mers never span contigs.
#[derive(Clone, Debug, PartialEq)]
pub struct GenomeSeq {
    pub contigs: Vec<(usize, Vec<u8>)>,
}

impl GenomeSeq {
    /// Cumulative base offset of each contig (so a global base index uniquely maps
    /// back to a single contig and local position). `cum[i]` = sum of contig
    /// lengths < `i`; the final entry is the total base count.
    pub(crate) fn cumulative(&self) -> Vec<u64> {
        let mut cum = Vec::with_capacity(self.contigs.len() + 1);
        let mut acc = 0u64;
        cum.push(0);
        for (len, _) in &self.contigs {
            acc += *len as u64;
            cum.push(acc);
        }
        cum
    }

    /// Map a global base index (a k-mer start) to `(contig, local_start)`.
    pub(crate) fn locate(&self, cum: &[u64], global: u64) -> (usize, usize) {
        // contigs are few; a linear scan with the precomputed prefix sums is fine.
        let mut ci = 0;
        while ci + 1 < cum.len() && cum[ci + 1] <= global {
            ci += 1;
        }
        (ci, (global - cum[ci]) as usize)
    }
}

// --- genome sequence (de)serialization (2-bit packed) -----------------------

/// Read a FASTA/FASTQ into a packed `GenomeSeq` (4 bases per byte, same layout
/// as the on-disk format), using the sketcher's `BYTE_TO_SEQ` mapping.
/// Returns `None` if unreadable.
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
        let len = seq.len();
        let mut packed = vec![0u8; (len + 3) / 4];
        for (i, &b) in seq.iter().enumerate() {
            let code = BYTE_TO_SEQ[b as usize];
            packed[i / 4] |= (code & 3) << (2 * (i % 4));
        }
        contigs.push((len, packed));
    }
    Some(GenomeSeq { contigs })
}

/// Serialize a genome's contigs as: contig count, then per contig its base length
/// followed by the 2-bit-packed bases (4 bases/byte, contig-aligned).
fn write_genome_seq_block<W: Write>(w: &mut W, seq: &GenomeSeq) -> io::Result<()> {
    write_uvarint(w, seq.contigs.len() as u64)?;
    for (len, data) in &seq.contigs {
        write_uvarint(w, *len as u64)?;
        w.write_all(data)?;
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

pub(crate) fn read_genome_seq_block<R: Read>(r: &mut R) -> io::Result<GenomeSeq> {
    let n_contigs = read_uvarint(r)? as usize;
    let mut contigs = Vec::with_capacity(n_contigs);
    for _ in 0..n_contigs {
        let len = read_uvarint(r)? as usize;
        let nbytes = (len + 3) / 4;
        let mut packed = vec![0u8; nbytes];
        r.read_exact(&mut packed)?;
        contigs.push((len, packed));
    }
    Ok(GenomeSeq { contigs })
}

/// A stable fingerprint so a compressed sample can be matched to its reference DB.
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

struct OwnAccum {
    rep_count: u32,
    rep_id: u32,
    strain_count: u32,
    strain_id: u32,
}

const POOL: u32 = u32::MAX;

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
pub(crate) fn sparse_threshold(sparse_c: usize) -> u64 {
    u64::MAX / sparse_c.max(1) as u64
}

#[inline]
pub(crate) fn keep_sparse_hash(h: u64, sparse_c: usize) -> bool {
    h < sparse_threshold(sparse_c)
}

#[inline]
pub(crate) fn sparse_naive_ani(matches: u32, domain: u32, k: usize) -> f64 {
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
    pub(crate) sparse_count: u32,
    dense_offset: u64,
    pub(crate) dense_domain: u32,
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
    pub(crate) fingerprint: u64,
    pub genomes: Vec<RefGenomeMeta>,
    /// stage-1: sparse distinctive hash -> owning genome id via MPHF slot arrays.
    sparse_mphf: Mphf<u64>,
    sparse_fingerprints: Vec<u32>,
    sparse_owners: Vec<u32>,
    /// shared pool, loaded once.
    pub(crate) pool: Vec<u64>,
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
    pub(crate) fn load_genome(&self, g: u32) -> io::Result<Arc<Vec<u64>>> {
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
    pub(crate) fn load_genome_seq(&self, g: u32) -> io::Result<Option<Arc<GenomeSeq>>> {
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
    pub(crate) fn pool_index(&self, h: u64) -> Option<u32> {
        let slot = self.pool_mphf.try_hash(&h)? as usize;
        if slot < self.pool.len() && self.pool[slot] == h {
            Some(slot as u32)
        } else {
            None
        }
    }

    pub(crate) fn ensure_can_decompress(&self) -> io::Result<()> {
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

/// Open a `.sylref` for querying/compression (CLI helper).
pub(super) fn open_refdb_file(path: &str) -> RefIndex {
    open_refdb_file_with_mode(path, false)
}

pub(super) fn open_refdb_file_for_compress(path: &str) -> RefIndex {
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

pub fn run_ref_build(args: RefBuildArgs) {
    super::init_logger(args.trace);
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
