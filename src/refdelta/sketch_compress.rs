//! Compression of read sketches into the reference-delta `.sylspr` format.

use super::ref_build::{open_refdb_file, open_refdb_file_for_compress, GenomeSeq, RefIndex};
use super::{
    substituted_hash, window_fr, SCHEME_ABSENT_RICE, SCHEME_BITMASK, SCHEME_PRESENT_RICE,
    SKETCH_MAGIC, SKETCH_VERSION, ZSTD_LEVEL,
};
use crate::cmdline::RefCompressArgs;
use crate::compress::{
    self, read_seq_sketch_compressed, read_string, read_uvarint, write_hashes, write_string,
    write_uvarint,
};
use crate::constants::*;
use crate::seeding::rev_mm_hash64;
use crate::types::*;
use fxhash::{FxHashMap, FxHashSet};
use log::*;
use rayon::prelude::*;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};
use std::time::Instant;

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

// --- adaptive present/absent subset coding ----------------------------------

/// Encode the sorted `present` indices into `domain` using the smallest of a
/// bitmask, present-Rice, or absent(complement)-Rice. Self-delimiting given the
/// `domain`, which the decoder knows from the reference DB.
pub(crate) fn encode_subset(out: &mut Vec<u8>, present: &[u64], domain: u64) -> io::Result<()> {
    // present-Rice (cheap: O(present))
    let mut p_rice = Vec::new();
    write_hashes(&mut p_rice, present)?;
    let bm_len = domain.div_ceil(8) as usize;

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

// --- single-substitution error k-mer encoding --------------------------------

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

/// Novel k-mers are processed in chunks of this size in `find_error_kmers`.
/// Exposed at module level so the caller can estimate work before invoking.
const NOVEL_CHUNK: usize = 8_000_000;

/// Max chunk-genome scan pairs before skipping error-kmer classification.
/// Each pair costs ~0.7 CPU-s (3.5 Mbp × 2 probes × 100 ns); 40 000 pairs ≈
/// 30 min at 4 threads.  Samples with very large novel sets (e.g. 292 M novel
/// k-mers × 6642 genomes) gain negligible compression from error reclassification
/// and should just skip it.
const MAX_ERROR_SCAN_PAIRS: usize = 40_000;

/// Locate `key`'s block and the bit mask to set/test within the blocked bloom
/// that prefilters the per-chunk novel-k-mer index. Build and query MUST derive
/// bits identically, so both go through here.
///
/// Four bits/key (vs the classic two) cuts the false-positive rate several-fold
/// at no extra memory access: every bit lives in the single 64-bit block that is
/// already loaded, so it costs only a couple of ALU ops. Two multiplicative
/// hashes supply four well-separated 6-bit fields. The block index and bit count
/// together give ~16 bits/key (see `build` below), for a false-positive rate
/// under ~0.5% instead of the prior ~5%. A false positive only wastes one
/// `index` lookup, so the lower rate trims that cache-miss tail.
#[inline]
fn bloom_locate(key: u32, n_block_mask: usize) -> (usize, u64) {
    let h1 = (key as u64).wrapping_mul(0x9e37_79b9_7f4a_7c15);
    let h2 = (key as u64).wrapping_mul(0xc2b2_ae3d_27d4_eb4f);
    let block = (h1 >> 6) as usize & n_block_mask;
    let bmask = (1u64 << (h1 & 63))
        | (1u64 << ((h1 >> 32) & 63))
        | (1u64 << (h2 & 63))
        | (1u64 << ((h2 >> 32) & 63));
    (block, bmask)
}

/// Try to express each `novel` hash as a single-base substitution of a hit
/// **representative** genome's k-mer (sequencing errors produce exactly such
/// k-mers). Reconstruction is deterministic and verified by hash equality, so a
/// match is never a false positive. Returns the matched entries grouped by genome
/// (already sorted by genomic position) and the set of consumed novel hashes.
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
    let k2 = k / 2;
    let shift_hi = 2 * k2 as u32;
    let mask_full = low_bits_mask(2 * k as u32);
    let mask_low = low_bits_mask(shift_hi);

    let mut sorted_hits = hits.to_vec();
    sorted_hits.sort_unstable();
    let eligible: Vec<u32> = sorted_hits
        .into_iter()
        .filter(|&g| idx.genomes[g as usize].is_rep)
        .collect();

    // Phase 1 (I/O): load all eligible genome sequences once so they can be
    // reused across every novel-k-mer chunk below without re-seeking.
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

    // Phase 2 (CPU): scan in chunks of NOVEL_CHUNK novel k-mers.
    let mut per_genome_matches: FxHashMap<u32, Vec<(u64, ErrorEntry)>> = FxHashMap::default();

    for chunk in novel.chunks(NOVEL_CHUNK) {
        // Index each novel hash under both half-k-mer keys: a single-base
        // substitution preserves at least one half, so the true genome
        // neighbour is guaranteed to share a key.
        let mut index: FxHashMap<u32, Vec<(u64, u64)>> = FxHashMap::default();
        index.reserve(chunk.len() * 4);
        for &h in chunk {
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

        // Blocked bloom prefilter over the index keys: ~4 keys per 64-bit block
        // (16 bits/key) with 4 bits set per key (see `bloom_locate`).
        let n_blocks = (index.len() / 4 + 1)
            .next_power_of_two()
            .min(128 * 1024 * 1024 / 8);
        let n_block_mask = n_blocks - 1;
        let mut bloom = vec![0u64; n_blocks];
        for &key in index.keys() {
            let (block, bmask) = bloom_locate(key, n_block_mask);
            bloom[block] |= bmask;
        }

        let chunk_results: Vec<(u32, Vec<(u64, ErrorEntry)>)> = genome_seqs
            .par_iter()
            .map(|(g, seq)| {
                let mut matches: Vec<(u64, ErrorEntry)> = Vec::new();
                let mut local_consumed: FxHashSet<u64> = FxHashSet::default();
                let mut base_global = 0u64;
                for (clen, contig) in &seq.contigs {
                    let clen = *clen;
                    if clen >= k {
                        let (mut f, _) = window_fr(contig, 0, k);
                        for start in 0..=(clen - k) {
                            if start > 0 {
                                let pos = start + k - 1;
                                let base = (contig[pos / 4] >> (2 * (pos % 4))) & 3;
                                f = ((f << 2) | base as u64) & mask_full;
                            }
                            for key in [(f >> shift_hi) as u32, (f & mask_low) as u32] {
                                let (block, bmask) = bloom_locate(key, n_block_mask);
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

        for (g, matches) in chunk_results {
            if !matches.is_empty() {
                per_genome_matches.entry(g).or_default().extend(matches);
            }
        }
    }

    // Phase 3 (sequential): cross-genome dedup in sorted eligible order.
    for &g in &eligible {
        let Some(matches) = per_genome_matches.remove(&g) else {
            continue;
        };
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

// --- compress_seq* entry points ----------------------------------------------

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
    compress_seq_with_screen_ani_and_error_kmers(
        inner,
        sketch,
        idx,
        ref_db_name,
        meta,
        ref_screen_ani,
        min_dense_kmers_for_error,
        true,
    )
}

pub fn compress_seq_with_screen_ani_and_error_kmers<W: Write>(
    inner: W,
    sketch: &SequencesSketch,
    idx: &RefIndex,
    ref_db_name: &str,
    meta: ReadSketchMeta,
    ref_screen_ani: f64,
    min_dense_kmers_for_error: usize,
    enable_error_kmers: bool,
) -> io::Result<()> {
    compress_seq_with_screen_ani_and_telemetry(
        inner,
        sketch,
        idx,
        ref_db_name,
        meta,
        ref_screen_ani,
        min_dense_kmers_for_error,
        enable_error_kmers,
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
    enable_error_kmers: bool,
) -> io::Result<Vec<RefCompressTelemetry>> {
    use super::ref_build::sparse_naive_ani;

    let total_start = Instant::now();
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
    drop(map);
    novel.shrink_to_fit();

    let error_start = Instant::now();
    let error_by_genome = if enable_error_kmers && idx.has_genome_seqs() && !novel.is_empty() {
        let error_hits: Vec<u32> = hits
            .iter()
            .copied()
            .filter(|g| {
                per_genome.get(g).map(|v| v.len()).unwrap_or(0) >= min_dense_kmers_for_error
            })
            .collect();
        let n_novel_chunks = novel.len().div_ceil(NOVEL_CHUNK);
        let work_pairs = n_novel_chunks.saturating_mul(error_hits.len());
        if work_pairs > MAX_ERROR_SCAN_PAIRS {
            info!(
                "ref-compress stage4 skip: {} novel chunks × {} eligible genomes = {} pairs \
                 exceeds limit {}; error-kmer scan skipped",
                n_novel_chunks,
                error_hits.len(),
                work_pairs,
                MAX_ERROR_SCAN_PAIRS
            );
            Vec::new()
        } else {
            let (by_genome, consumed) =
                find_error_kmers(idx, &error_hits, &novel, sketch.c, sketch.k)?;
            if !consumed.is_empty() {
                novel.retain(|h| !consumed.contains(h));
            }
            by_genome
        }
    } else {
        if !enable_error_kmers {
            info!("ref-compress stage4 skip: error-kmer encoding disabled");
        }
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

    // hit genomes: sorted global ids, delta-coded
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

// --- CLI helpers: inspect, verify, run_ref_compress -------------------------

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

fn verify_ref_sketch(
    path: &std::path::Path,
    original: &SequencesSketch,
    idx: &RefIndex,
) -> io::Result<()> {
    use super::decompress_seq;
    let r = BufReader::with_capacity(
        10_000_000,
        File::open(path).unwrap_or_else(|_| panic!("Could not open {:?}", path)),
    );
    let decoded = decompress_seq(r, idx)?;
    compare_seq_sketches(original, &decoded)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn original_path_candidates(ref_path: &std::path::Path, sample_file: &str) -> Vec<PathBuf> {
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

fn find_original_sketch_path(ref_path: &std::path::Path, sample_file: &str) -> io::Result<PathBuf> {
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
    super::init_logger(args.trace);
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
            use super::decompress_seq;
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
                !args.no_error_kmer,
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
