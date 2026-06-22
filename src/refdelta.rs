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
use crate::seeding::{canonical_kmer, mm_hash64, mutate_kmer, rev_hash_64};
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
const REFDB_VERSION: u8 = 5; // v4 = no all-kmers index; v5 = with all-kmers index
const SKETCH_MAGIC: &[u8; 4] = b"SYLD"; // reference-Delta sample
const SKETCH_VERSION: u8 = 4; // v3 = no error kmers; v4 = with error-kmer section

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
const ALLKMERS_MPHF_GAMMA: f64 = 1.7;

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

/// Build and serialize an MPHF over `all_keys` (sorted, deduplicated), returning the
/// MPHF and a slot-indexed hash array for zero-false-positive verification.
fn build_allkmers_index(all_keys: &[u64]) -> io::Result<(Mphf<u64>, Vec<u64>)> {
    let mphf = Mphf::new_parallel(ALLKMERS_MPHF_GAMMA, all_keys, Some(0));
    let n = all_keys.len();
    let mut by_slot = vec![0u64; n];
    let mut filled = vec![false; n];
    for &h in all_keys {
        let slot = mphf.hash(&h) as usize;
        if slot >= n || filled[slot] {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "all-kmers MPHF construction produced an invalid slot",
            ));
        }
        by_slot[slot] = h;
        filled[slot] = true;
    }
    Ok((mphf, by_slot))
}

fn write_allkmers_block<W: Write>(w: &mut W, mphf: &Mphf<u64>, by_slot: &[u64]) -> io::Result<()> {
    let mphf_bytes = bincode::serialize(mphf)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e.to_string()))?;
    write_uvarint(w, mphf_bytes.len() as u64)?;
    w.write_all(&mphf_bytes)?;
    for &h in by_slot {
        w.write_all(&h.to_le_bytes())?;
    }
    Ok(())
}

fn read_allkmers_block<R: Read>(r: &mut R, count: usize) -> io::Result<(Mphf<u64>, Vec<u64>)> {
    let mphf_len = read_uvarint(r)? as usize;
    let mut mphf_bytes = vec![0u8; mphf_len];
    r.read_exact(&mut mphf_bytes)?;
    let mphf: Mphf<u64> = bincode::deserialize(&mphf_bytes)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e.to_string()))?;
    let mut by_slot = Vec::with_capacity(count);
    for _ in 0..count {
        let mut buf = [0u8; 8];
        r.read_exact(&mut buf)?;
        by_slot.push(u64::from_le_bytes(buf));
    }
    Ok((mphf, by_slot))
}

/// For a novel read k-mer (with canonical raw form `read_raw`), search the all-kmers MPHF for a
/// reference k-mer that is one substitution away. Returns `(ref_slot, pos_in_ref, new_base)` —
/// the mutation to apply to the *reference* k-mer to reconstruct the read k-mer.
fn find_ref_error_neighbor(
    read_raw: u64,
    k: usize,
    mphf: &Mphf<u64>,
    by_slot: &[u64],
) -> Option<(u64, u8, u8)> {
    for pos in 0..k {
        let read_base = ((read_raw >> (2 * pos)) & 3) as u8;
        for new_base in 0..4u8 {
            if new_base == read_base {
                continue;
            }
            let neighbor_raw = mutate_kmer(read_raw, pos, new_base);
            let canonical_nb = canonical_kmer(neighbor_raw, k);
            let nb_hash = mm_hash64(canonical_nb);
            if let Some(slot) = mphf.try_hash(&nb_hash) {
                let slot = slot as usize;
                if slot < by_slot.len() && by_slot[slot] == nb_hash {
                    // Found the reference k-mer. Find which mutation on the reference
                    // recovers the read k-mer (needed for lossless decompression).
                    let ref_raw = rev_hash_64(nb_hash);
                    if let Some((rpos, rbase)) = find_reverse_mutation(ref_raw, read_raw, k) {
                        return Some((slot as u64, rpos, rbase));
                    }
                }
            }
        }
    }
    None
}

/// Find `(pos, new_base)` such that `canonical_kmer(mutate_kmer(ref_raw, pos, new_base), k) == read_raw`.
/// Always succeeds for a valid Hamming-1 pair; returns None only if called incorrectly.
fn find_reverse_mutation(ref_raw: u64, read_raw: u64, k: usize) -> Option<(u8, u8)> {
    for pos in 0..k {
        let ref_base = ((ref_raw >> (2 * pos)) & 3) as u8;
        for new_base in 0..4u8 {
            if new_base == ref_base {
                continue;
            }
            let mutated = mutate_kmer(ref_raw, pos, new_base);
            if canonical_kmer(mutated, k) == read_raw {
                return Some((pos as u8, new_base));
            }
        }
    }
    None
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

/// Write a `.sylref` in the seekable two-stage layout. `sparse_div` controls the
/// stage-1 subsampling rate (keep 1/`sparse_div` of each genome's distinctive
/// k-mers); `sparse_div <= 1` keeps all of them. When `all_kmers_index` is true, an
/// additional MPHF over every reference k-mer is embedded (v5 format); otherwise v4.
pub fn write_refdb<W: Write>(
    mut w: W,
    db: &RefDb,
    sparse_div: u64,
    all_kmers_index: bool,
) -> io::Result<()> {
    let ng = db.genomes.len();
    let sp_ids = species_ids(&db.genomes);

    // Body: [sparse section][MPHF sparse index][MPHF pool block][dense blocks...]
    let mut body: Vec<u8> = Vec::new();

    let sparse_off = body.len() as u64; // 0
    let mut sparse_count = vec![0usize; ng];
    let mut sparse_keys: Vec<u64> = Vec::new();
    let mut sparse_owners: Vec<u32> = Vec::new();
    for g in 0..ng {
        let mut cnt = 0usize;
        for &h in &db.distinctive[g] {
            if sparse_div <= 1 || h % sparse_div == 0 {
                body.extend_from_slice(&h.to_le_bytes());
                sparse_keys.push(h);
                sparse_owners.push(g as u32);
                cnt += 1;
            }
        }
        sparse_count[g] = cnt;
    }

    let sparse_index_off = body.len() as u64;
    write_sparse_index_block(&mut body, &sparse_keys, &sparse_owners)?;
    drop(sparse_keys);
    drop(sparse_owners);

    let pool_off = body.len() as u64;
    write_pool_block(&mut body, &db.pool)?;

    // Dense blocks store only the *complement* of the stage-1 sparse subset, so a
    // k-mer kept in stage 1 is not also stored here; load merges the two back.
    let mut dense_off = vec![0u64; ng];
    for g in 0..ng {
        dense_off[g] = body.len() as u64;
        let complement: Vec<u64> = db.distinctive[g]
            .iter()
            .copied()
            .filter(|&h| !(sparse_div <= 1 || h % sparse_div == 0))
            .collect();
        write_hashes(&mut body, &complement)?;
    }

    // Optional all-kmers MPHF block (v5 only).
    let allkmers_info: Option<(u64, usize)> = if all_kmers_index {
        let mut all_keys: Vec<u64> = db
            .distinctive
            .iter()
            .flatten()
            .copied()
            .chain(db.pool.iter().copied())
            .collect();
        all_keys.sort_unstable();
        all_keys.dedup();
        let n = all_keys.len();
        let (mphf, by_slot) = build_allkmers_index(&all_keys)?;
        let off = body.len() as u64;
        write_allkmers_block(&mut body, &mphf, &by_slot)?;
        Some((off, n))
    } else {
        None
    };

    // Footer (zstd-compressed): metadata + absolute offsets into the file.
    let mut footer: Vec<u8> = Vec::new();
    write_uvarint(&mut footer, db.c as u64)?;
    write_uvarint(&mut footer, db.k as u64)?;
    write_uvarint(&mut footer, sparse_div)?;
    footer.extend_from_slice(&fingerprint(db).to_le_bytes());
    write_uvarint(&mut footer, ng as u64)?;
    write_uvarint(&mut footer, HEADER_LEN + sparse_off)?;
    for g in 0..ng {
        write_uvarint(&mut footer, sparse_count[g] as u64)?;
    }
    write_uvarint(&mut footer, HEADER_LEN + sparse_index_off)?;
    write_uvarint(&mut footer, HEADER_LEN + pool_off)?;
    write_uvarint(&mut footer, db.pool.len() as u64)?;
    for g in 0..ng {
        write_string(&mut footer, &db.genomes[g].file_name)?;
        write_string(&mut footer, &db.genomes[g].species)?;
        footer.push(db.genomes[g].is_rep as u8);
        write_uvarint(&mut footer, sp_ids[g] as u64)?;
        write_uvarint(&mut footer, HEADER_LEN + dense_off[g])?;
        write_uvarint(&mut footer, db.distinctive[g].len() as u64)?;
    }
    // v5 extension: all-kmers index location (0 = not present).
    if let Some((off, count)) = allkmers_info {
        write_uvarint(&mut footer, HEADER_LEN + off)?;
        write_uvarint(&mut footer, count as u64)?;
    } else {
        write_uvarint(&mut footer, 0u64)?; // 0 = no all-kmers index
        write_uvarint(&mut footer, 0u64)?;
    }
    let footer_comp = zstd::stream::encode_all(&footer[..], 9)?;
    let footer_offset = HEADER_LEN + body.len() as u64;

    let version = if all_kmers_index { 5u8 } else { 4u8 };
    w.write_all(REFDB_MAGIC)?;
    w.write_all(&[version])?;
    w.write_all(&footer_offset.to_le_bytes())?;
    w.write_all(&body)?;
    w.write_all(&footer_comp)?;
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
    sparse_div: u64,
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
    compression_only: bool,
    /// All-kmers index: MPHF over every reference k-mer (distinctive + pool), built with
    /// `ref-build --all-kmers-index`. When present, enables 1-error k-mer recovery during
    /// compression. Also required to decompress v4 sketches.
    all_kmers_mphf: Option<Mphf<u64>>,
    /// Hash value at each MPHF slot (slot → hash), used to verify membership and to
    /// reconstruct the reference k-mer sequence for error k-mer decoding.
    all_kmers_hashes: Option<Vec<u64>>,
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
    let db_ver = hdr[4];
    if db_ver < 4 || db_ver > REFDB_VERSION {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("unsupported reference DB version {} (supported: 4-{})", db_ver, REFDB_VERSION),
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
    let sparse_div = read_uvarint(&mut f)?;
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
        });
    }

    // v5 extension: all-kmers index location (written for both v4 and v5 refdbs,
    // but only non-zero when an all-kmers index was built).
    let (allkmers_offset, allkmers_count) = if db_ver >= 5 {
        let off = read_uvarint(&mut f)?;
        let cnt = read_uvarint(&mut f)? as usize;
        (off, cnt)
    } else {
        (0u64, 0usize)
    };

    let total_sparse: usize = sparse_count.iter().sum();
    r.seek(SeekFrom::Start(sparse_index_offset))?;
    let (sparse_mphf, sparse_fingerprints, sparse_owners) =
        read_sparse_index_block(&mut r, total_sparse)?;

    // shared pool (loaded once). Hashes are stored in MPHF slot order, so
    // compression can compute a candidate slot and verify it without a HashMap.
    r.seek(SeekFrom::Start(pool_offset))?;
    let (pool_mphf, pool) = read_pool_block(&mut r, pool_domain)?;

    // All-kmers index (optional): loaded eagerly when present.
    let (all_kmers_mphf, all_kmers_hashes) = if allkmers_offset > 0 && allkmers_count > 0 {
        r.seek(SeekFrom::Start(allkmers_offset))?;
        let (mphf, by_slot) = read_allkmers_block(&mut r, allkmers_count)?;
        (Some(mphf), Some(by_slot))
    } else {
        (None, None)
    };

    Ok(RefIndex {
        c,
        k,
        sparse_div,
        fingerprint,
        genomes,
        sparse_mphf,
        sparse_fingerprints,
        sparse_owners,
        pool,
        pool_mphf,
        reader: Mutex::new(Box::new(r)),
        cache: Mutex::new(FxHashMap::default()),
        compression_only: matches!(mode, RefOpenMode::CompressionOnly),
        all_kmers_mphf,
        all_kmers_hashes,
    })
}

impl RefIndex {
    /// Stage-1 query: the genomes a sample contains, by sparse-hit count.
    pub fn hit_genomes(&self, sketch: &SequencesSketch) -> Vec<u32> {
        let mut counts: FxHashMap<u32, u32> = FxHashMap::default();
        for &h in sketch.kmer_counts.keys() {
            if self.sparse_div > 1 && h % self.sparse_div != 0 {
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
            .filter(|&(_, c)| c >= SPARSE_MIN_HITS)
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

/// Compress a read sketch against the reference index. Only the sample's hit
/// genomes' dense blocks are loaded; the pool is already resident in `idx`.
pub fn compress_seq<W: Write>(
    inner: W,
    sketch: &SequencesSketch,
    idx: &RefIndex,
    ref_db_name: &str,
) -> io::Result<()> {
    compress_seq_with_meta(inner, sketch, idx, ref_db_name, ReadSketchMeta::default())
}

pub fn compress_seq_with_meta<W: Write>(
    inner: W,
    sketch: &SequencesSketch,
    idx: &RefIndex,
    ref_db_name: &str,
    meta: ReadSketchMeta,
) -> io::Result<()> {
    // stage 1 -> stage 2: load the distinctive blocks of the hit genomes and
    // build a per-sample hash -> (genome, index) lookup over just those.
    let hits = idx.hit_genomes(sketch);
    let mut map: FxHashMap<u64, (u32, u32)> = FxHashMap::default();
    for &g in &hits {
        let arr = idx.load_genome(g)?;
        for (i, &h) in arr.iter().enumerate() {
            map.insert(h, (g, i as u32));
        }
    }

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

    // Error k-mers: novel hashes that are 1 substitution away from a reference k-mer.
    // Only possible when the reference DB was built with --all-kmers-index (v5).
    // Encoded as delta-coded MPHF slot + 1-byte packed (pos<<2|base) for the mutation
    // that, applied to the *reference* k-mer, reconstructs the read k-mer.
    let mut error_section: Vec<u8> = Vec::new();
    let truly_novel: Vec<u64>;
    if let (Some(mphf), Some(by_slot)) = (&idx.all_kmers_mphf, &idx.all_kmers_hashes) {
        let mut enc_slots: Vec<(u64, u8, u8)> = Vec::new(); // (slot, pos_in_ref, new_base)
        let mut still_novel: Vec<u64> = Vec::new();
        for &h in &novel {
            let read_raw = rev_hash_64(h);
            if let Some((slot, rpos, rbase)) =
                find_ref_error_neighbor(read_raw, idx.k, mphf, by_slot)
            {
                enc_slots.push((slot, rpos, rbase));
            } else {
                still_novel.push(h);
            }
        }
        enc_slots.sort_unstable_by_key(|&(s, _, _)| s);
        write_uvarint(&mut error_section, enc_slots.len() as u64)?;
        let mut prev_slot = 0u64;
        for (slot, rpos, rbase) in enc_slots {
            write_uvarint(&mut error_section, slot - prev_slot)?;
            prev_slot = slot;
            error_section.push((rpos << 2) | (rbase & 3));
        }
        truly_novel = still_novel;
    } else {
        truly_novel = novel;
    };

    // Truly novel hashes (Rice) — those unexplained even by 1-error matching.
    let mut novel_section = Vec::new();
    write_hashes(&mut novel_section, &truly_novel)?;

    // counts, in ascending-hash order (reproducible on decode)
    let mut count_section = Vec::new();
    let mut keys: Vec<u64> = sketch.kmer_counts.keys().copied().collect();
    keys.sort_unstable();
    for h in &keys {
        write_uvarint(&mut count_section, sketch.kmer_counts[h] as u64)?;
    }

    // Use sketch version 4 when an error section is present (even if empty, for
    // consistency), 3 otherwise.
    let sketch_ver = if idx.all_kmers_mphf.is_some() { 4u8 } else { 3u8 };

    let mut payload = Vec::new();
    payload.extend_from_slice(SKETCH_MAGIC);
    payload.push(sketch_ver);
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
    write_uvarint(&mut payload, truly_novel.len() as u64)?;
    write_uvarint(&mut payload, hit_section.len() as u64)?;
    write_uvarint(&mut payload, pool_section.len() as u64)?;
    write_uvarint(&mut payload, novel_section.len() as u64)?;
    write_uvarint(&mut payload, count_section.len() as u64)?;
    if sketch_ver >= 4 {
        // error_section starts with the count (uvarint); store byte length too for inspect.
        write_uvarint(&mut payload, error_section.len() as u64)?;
    }
    payload.extend_from_slice(&hit_section);
    payload.extend_from_slice(&pool_section);
    payload.extend_from_slice(&novel_section);
    if sketch_ver >= 4 {
        payload.extend_from_slice(&error_section);
    }
    payload.extend_from_slice(&count_section);

    let mut enc = zstd::stream::write::Encoder::new(inner, ZSTD_LEVEL)?;
    enc.write_all(&payload)?;
    enc.finish()?;
    Ok(())
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

    // Header metadata: 8 uvarints for v2-3; 9 for v4 (adds error section byte length).
    if ver[0] >= 2 {
        let n_meta = if ver[0] >= 4 { 9 } else { 8 };
        for _ in 0..n_meta {
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

    // v4: error k-mers section — reconstruct read k-mers from (ref_slot, pos, base).
    if ver[0] >= 4 {
        let by_slot = idx.all_kmers_hashes.as_deref().ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "sketch was compressed with error-kmer matching (v4) but the reference DB has \
                 no all-kmers index; rebuild with `ref-build --all-kmers-index`",
            )
        })?;
        let n_error = read_uvarint(&mut r)? as usize;
        let mut prev_slot = 0u64;
        for _ in 0..n_error {
            let delta = read_uvarint(&mut r)?;
            let slot = prev_slot + delta;
            prev_slot = slot;
            let mut packed = [0u8; 1];
            r.read_exact(&mut packed)?;
            let rpos = (packed[0] >> 2) as usize;
            let rbase = packed[0] & 3;
            let ref_hash = by_slot[slot as usize];
            let ref_raw = rev_hash_64(ref_hash);
            let mutated_raw = mutate_kmer(ref_raw, rpos, rbase);
            let canonical_mutated = canonical_kmer(mutated_raw, k);
            hashes.push(mm_hash64(canonical_mutated));
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
        error_section_bytes,
    ) = if ver[0] >= 2 {
        let hg = Some(read_uvarint(&mut r)?);
        let ag = Some(read_uvarint(&mut r)?);
        let sp = Some(read_uvarint(&mut r)?);
        let nv = Some(read_uvarint(&mut r)?);
        let hs = Some(read_uvarint(&mut r)?);
        let ps = Some(read_uvarint(&mut r)?);
        let ns = Some(read_uvarint(&mut r)?);
        let cs = Some(read_uvarint(&mut r)?);
        let es = if ver[0] >= 4 { Some(read_uvarint(&mut r)?) } else { None };
        (
            hg, ag, sp, nv, hs, ps, ns, cs, es,
        )
    } else {
        (None, None, None, None, None, None, None, None, None)
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
        error_section_bytes,
    })
}

fn opt_u64(x: Option<u64>) -> String {
    x.map(|v| v.to_string()).unwrap_or_else(|| "-".to_string())
}

fn run_ref_inspect(files: &[String]) {
    println!(
        "file\tversion\tcompressed_bytes\tpayload_bytes\treference_fingerprint\treference_db\tsample_file\tsample_name\tc\tk\tpaired\tmean_read_length\tnum_reads\theader_metadata_payload_bytes\thit_genomes\tassigned_to_genomes\tshared_pool\tnovel\ttotal_hashes\thit_section_payload_bytes\tpool_section_payload_bytes\tnovel_section_payload_bytes\tcount_section_payload_bytes\terror_section_payload_bytes"
    );
    for f in files {
        let x = inspect_ref_sketch(f).unwrap_or_else(|e| panic!("Failed to inspect {}: {}", f, e));
        let total_hashes = match (x.assigned_to_genomes, x.shared_pool, x.novel) {
            (Some(a), Some(p), Some(n)) => (a + p + n).to_string(),
            _ => "-".to_string(),
        };
        println!(
            "{}\t{}\t{}\t{}\t{:016x}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
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
            total_hashes,
            opt_u64(x.hit_section_bytes),
            opt_u64(x.pool_section_bytes),
            opt_u64(x.novel_section_bytes),
            opt_u64(x.count_section_bytes),
            opt_u64(x.error_section_bytes),
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
    let sparse_div = args.sparse_div.max(1);
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

    let db = RefDb {
        c,
        k,
        genomes,
        distinctive,
        pool,
    };
    let n_dist: usize = db.distinctive.iter().map(|d| d.len()).sum();
    let n_sparse: usize = if sparse_div <= 1 {
        n_dist
    } else {
        db.distinctive
            .iter()
            .map(|d| d.iter().filter(|&&h| h % sparse_div == 0).count())
            .sum()
    };
    info!(
        "Reference: {} genomes, {} distinctive k-mers, {} shared-pool k-mers, {} stage-1 sparse k-mers (1/{})",
        db.genomes.len(),
        n_dist,
        db.pool.len(),
        n_sparse,
        sparse_div
    );

    let w =
        BufWriter::new(File::create(&out).unwrap_or_else(|_| panic!("Could not create {}", out)));
    write_refdb(w, &db, sparse_div, args.all_kmers_index)
        .unwrap_or_else(|e| panic!("Failed to write {}: {}", out, e));
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
            compress_seq(w, &sketch, &idx, ref_db)
                .unwrap_or_else(|e| panic!("Failed to compress {}: {}", f, e));
            let mut c = counter.lock().unwrap();
            *c += 1;
            info!("Compressed {} -> {:?}", f, out);
        });
    }
    info!("Done ({} sample(s)).", *counter.lock().unwrap());
}
