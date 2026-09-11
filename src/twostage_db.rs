//! Two-stage seekable genome database (`.syl2db`) for `profile --two-stage`.
//!
//! A standard sylph database (`.syldb`) is a single bincoded `Vec<GenomeSketch>`
//! that must be loaded in full to profile a sample. For two-stage profiling we
//! only ever need the *dense* k-mers of the handful of genomes a sample actually
//! contains; loading every genome's dense k-mers up front is wasteful for a large
//! reference. This module re-packs a `.syldb` into a two-stage seekable layout,
//! mirroring (in spirit) the `.sylref` format of wwood/sylph#2 but deliberately
//! simpler:
//!
//!   * **No k-mer dereplication.** Each genome keeps its *own complete* k-mer set;
//!     a k-mer shared by several genomes is stored in each of them.
//!   * **No shared/pooled hashes.** There is no conserved-k-mer pool.
//!
//! ## Layout
//!
//! ```text
//! [4]  magic  "SY2D"
//! [1]  version
//! [8]  checksum      (XXH64 of everything after the header, u64 LE)
//! [8]  index_offset  (u64 LE)
//! [8]  footer_offset (u64 LE)
//! ---- body ----
//! dense block 0, dense block 1, ...      (each: Golomb-Rice compressed)
//! ---- index ----
//! ScreenIndex                            (pooled-MPHF stage-1 screen index)
//! ---- footer ----
//! bincode(Footer)                        (per-genome metadata)
//! ```
//!
//! **Stage 1 (pooled MPHF, loaded fully).** The `ScreenIndex` is one minimal
//! perfect hash over the *distinct* sparse (`screen_c`) k-mers of all genomes,
//! with a multi-owner CSR mapping each sparse k-mer to the list of genomes that
//! carry it. A sample is screened by a single pass over its k-mers (work ∝
//! sample, not reference): each sample k-mer is looked up once and its coverage
//! pushed to every owning genome. The contained genomes are exactly those a
//! per-genome `get_stats` screen would find -- this is the "Path B" organisation
//! (see `experiments/7_mphf_screen_again`), but multi-owner because `.syl2db`
//! keeps shared k-mers.
//!
//! **Stage 2 (dense, Golomb-Rice, loaded on demand).** Each genome's *full*
//! `genome_kmers` and `pseudotax_tracked_nonused_kmers` are an independently
//! Golomb-Rice-coded block at a known offset. Only the genomes that pass the
//! stage-1 screen are decoded (and cached across samples) to reconstruct their
//! exact `GenomeSketch` for the dense profiling pass.
//!
//! ## Adding genomes
//!
//! Because each genome's dense block is independent and self-delimiting, growing a
//! database ([`run_db_add`]) does not have to re-encode anything: the existing
//! blocks are copied to the new file byte for byte and the new genomes' blocks are
//! appended after them. Only the stage-1 index has to be rebuilt, since its MPHF is
//! built over the pooled key set and cannot absorb new keys — and rebuilding it
//! needs every genome's *sparse* k-mers, which is why the copied blocks are still
//! decoded on the way past. Both writers share [`DbBuilder`], which streams blocks
//! out as they are produced instead of assembling the body in RAM.
//!
//! ## Integrity
//!
//! Profiling only ever touches a few blocks of the file, so a corrupt database is
//! not detected as a side effect of reading it (unlike the single-frame `.sylspc` /
//! `.sylspr` formats, whose zstd checksum is validated on every decode). The header
//! instead carries an XXH64 of the rest of the file, checked on demand by
//! [`TwoStageDb::verify_checksum`] — which `weebill inspect` does.
//!
//! ## Upstream sylph compatibility
//!
//! Upstream sylph writes the same format without that checksum (version 2), so its
//! header is 8 bytes shorter and its offsets sit 8 bytes earlier; everything after
//! the header is byte-identical in meaning. Both versions are read here, so a
//! sylph-built `.syl2db` — including a hosted database — works directly with
//! `profile --two-stage`, `db-add` and `inspect`, with only the integrity check
//! unavailable for a v2 file. Databases written here are version 3 (with the
//! checksum), which sylph 1.0 does not yet read.

use crate::cmdline::{DbAddArgs, DbConvertArgs};
use crate::constants::*;
use crate::types::*;
use boomphf::Mphf;
use fxhash::{FxHashMap, FxHashSet};
use log::*;
use rayon::prelude::*;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::os::unix::fs::FileExt;
use std::path::Path;
use std::sync::{Arc, Mutex};

const MAGIC: &[u8; 4] = b"SY2D";
/// Format version written here (dense Golomb-Rice blocks + pooled-MPHF stage-1 screen
/// index + whole-file checksum).
const VERSION: u8 = 3;
/// Upstream sylph's version of the same format, without the whole-file checksum.
/// Readable here (see [`parse_header`]); anything older must be rebuilt with
/// `db-convert`.
const VERSION_NO_CHECKSUM: u8 = 2;
/// magic (4) + version (1) + XXH64 of the rest of the file (8) + index offset (8)
/// + footer offset (8)
const HEADER_LEN: u64 = 29;
/// [`VERSION_NO_CHECKSUM`] header: as above without the checksum.
const HEADER_LEN_V2: u64 = 21;
/// boomphf construction gamma (space/speed trade-off), matching the ref-delta
/// sparse index.
const MPHF_GAMMA: f64 = 2.0;

// --- primitive integer / bit coding -----------------------------------------

fn write_uvarint(w: &mut Vec<u8>, mut x: u64) {
    loop {
        let mut byte = (x & 0x7f) as u8;
        x >>= 7;
        if x != 0 {
            byte |= 0x80;
        }
        w.push(byte);
        if x == 0 {
            break;
        }
    }
}

fn read_uvarint<R: Read>(r: &mut R) -> io::Result<u64> {
    let mut result: u64 = 0;
    let mut shift: u32 = 0;
    loop {
        let mut b = [0u8; 1];
        r.read_exact(&mut b)?;
        result |= ((b[0] & 0x7f) as u64) << shift;
        if b[0] & 0x80 == 0 {
            break;
        }
        shift += 7;
        if shift >= 64 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "uvarint overflow",
            ));
        }
    }
    Ok(result)
}

/// LSB-first bit writer accumulating into a byte buffer.
struct BitWriter {
    buf: Vec<u8>,
    cur: u8,
    nbits: u8,
}

impl BitWriter {
    fn new() -> Self {
        BitWriter {
            buf: Vec::new(),
            cur: 0,
            nbits: 0,
        }
    }
    #[inline]
    fn write_bit(&mut self, b: u32) {
        if b != 0 {
            self.cur |= 1 << self.nbits;
        }
        self.nbits += 1;
        if self.nbits == 8 {
            self.buf.push(self.cur);
            self.cur = 0;
            self.nbits = 0;
        }
    }
    #[inline]
    fn write_bits(&mut self, val: u64, n: u32) {
        for i in 0..n {
            self.write_bit(((val >> i) & 1) as u32);
        }
    }
    #[inline]
    fn write_unary(&mut self, q: u64) {
        for _ in 0..q {
            self.write_bit(1);
        }
        self.write_bit(0);
    }
    fn finish(mut self) -> Vec<u8> {
        if self.nbits > 0 {
            self.buf.push(self.cur);
        }
        self.buf
    }
}

/// LSB-first bit reader that decodes a word at a time: bytes are buffered into a
/// 64-bit accumulator so `read_bits`/`read_unary` extract many bits per
/// instruction (shift/mask, trailing_ones) instead of one bit per call.
struct BitReader<'a> {
    buf: &'a [u8],
    pos: usize,
    acc: u64,
    nbits: u32,
}

impl<'a> BitReader<'a> {
    fn new(buf: &'a [u8]) -> Self {
        BitReader {
            buf,
            pos: 0,
            acc: 0,
            nbits: 0,
        }
    }
    /// Pull bytes into the accumulator until it holds >= 56 bits (so it never
    /// exceeds 63, keeping shifts in range) or the input is exhausted.
    #[inline]
    fn refill(&mut self) {
        while self.nbits < 56 && self.pos < self.buf.len() {
            self.acc |= (self.buf[self.pos] as u64) << self.nbits;
            self.pos += 1;
            self.nbits += 8;
        }
    }
    #[inline]
    fn read_bits(&mut self, n: u32) -> io::Result<u64> {
        if n == 0 {
            return Ok(0);
        }
        if n > 32 {
            // Split so each half fits the (>=56 bit) accumulator comfortably.
            let lo = self.read_bits(32)?;
            let hi = self.read_bits(n - 32)?;
            return Ok(lo | (hi << 32));
        }
        if self.nbits < n {
            self.refill();
            if self.nbits < n {
                return Err(io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    "bitstream truncated",
                ));
            }
        }
        let v = self.acc & ((1u64 << n) - 1);
        self.acc >>= n;
        self.nbits -= n;
        Ok(v)
    }
    /// Unary = run of 1s terminated by a 0 (matches `BitWriter::write_unary`).
    #[inline]
    fn read_unary(&mut self) -> io::Result<u64> {
        let mut q = 0u64;
        loop {
            if self.nbits == 0 {
                self.refill();
                if self.nbits == 0 {
                    return Err(io::Error::new(
                        io::ErrorKind::UnexpectedEof,
                        "bitstream truncated",
                    ));
                }
            }
            let ones = (self.acc | (1u64 << self.nbits)).trailing_ones(); // <= nbits
            if ones >= self.nbits {
                // all buffered bits are 1s; consume them and continue
                q += self.nbits as u64;
                self.acc = 0;
                self.nbits = 0;
            } else {
                // `ones` 1-bits then the terminating 0
                q += ones as u64;
                self.acc >>= ones + 1;
                self.nbits -= ones + 1;
                return Ok(q);
            }
        }
    }
}

/// Sort + delta + Golomb-Rice encode a set of hashes onto `out`. Order is not
/// preserved (hash sets are order-independent); duplicates become zero gaps and
/// are preserved. The Rice parameter is chosen from the mean gap and written
/// inline, so the block is self-delimiting given the leading count.
fn write_hashes(out: &mut Vec<u8>, hashes: &[u64]) {
    let mut sorted = hashes.to_vec();
    sorted.sort_unstable();
    write_uvarint(out, sorted.len() as u64);
    if sorted.is_empty() {
        return;
    }
    let mut deltas = Vec::with_capacity(sorted.len());
    let mut prev = 0u64;
    for &h in &sorted {
        deltas.push(h - prev);
        prev = h;
    }
    // Rice parameter k ~ log2(mean gap): near-optimal for the geometric gap
    // distribution of uniformly random hashes.
    let sum: u128 = deltas.iter().map(|&d| d as u128).sum();
    let mean = (sum / deltas.len() as u128).max(1);
    let mut k = 0u32;
    while k < 63 && (1u128 << (k + 1)) <= mean {
        k += 1;
    }
    out.push(k as u8);
    let mut bw = BitWriter::new();
    for &d in &deltas {
        bw.write_unary(d >> k);
        if k > 0 {
            bw.write_bits(d & ((1u64 << k) - 1), k);
        }
    }
    let bits = bw.finish();
    write_uvarint(out, bits.len() as u64);
    out.extend_from_slice(&bits);
}

fn read_hashes<R: Read>(r: &mut R) -> io::Result<Vec<u64>> {
    let n = read_uvarint(r)? as usize;
    if n == 0 {
        return Ok(Vec::new());
    }
    let mut kb = [0u8; 1];
    r.read_exact(&mut kb)?;
    let k = kb[0] as u32;
    let blen = read_uvarint(r)? as usize;
    let mut bits = vec![0u8; blen];
    r.read_exact(&mut bits)?;
    let mut br = BitReader::new(&bits);
    let mut out = Vec::with_capacity(n);
    let mut prev = 0u64;
    for _ in 0..n {
        let q = br.read_unary()?;
        let low = if k > 0 { br.read_bits(k)? } else { 0 };
        let d = (q << k) | low;
        prev = prev.wrapping_add(d);
        out.push(prev);
    }
    Ok(out)
}

// --- footer (stage-1 sparse index + metadata) -------------------------------

/// Per-genome metadata. Everything here is loaded into memory when the database
/// is opened; the dense block at `dense_offset` is decoded lazily. The stage-1
/// sparse k-mers live in the pooled `ScreenIndex`, not here.
#[derive(serde::Serialize, serde::Deserialize, Clone, Debug, PartialEq)]
pub struct GenomeMeta {
    pub file_name: String,
    pub first_contig_name: String,
    pub gn_size: usize,
    pub min_spacing: usize,
    pub has_pseudotax: bool,
    pub dense_offset: u64,
}

#[derive(serde::Serialize, serde::Deserialize, Clone, Debug, PartialEq)]
pub struct Footer {
    /// Dense rate: every k-mer is kept in the dense blocks at this `-c`.
    pub c: usize,
    pub k: usize,
    /// Sparse stage-1 screen rate (`screen_c >= c`).
    pub screen_c: usize,
    pub genomes: Vec<GenomeMeta>,
}

// --- stage-1 pooled-MPHF screen index ---------------------------------------

#[inline]
fn sparse_fingerprint(h: u64) -> u32 {
    (h ^ (h >> 32)) as u32
}

/// FracMinHash threshold for the stage-1 screen: a k-mer is "sparse" iff its
/// hash is `< u64::MAX / screen_c` (the same rule as `subsample_view` and the
/// `.syl2db` build).
#[inline]
pub(crate) fn screen_threshold(screen_c: usize) -> u64 {
    u64::MAX / screen_c.max(1) as u64
}

/// Pooled stage-1 screen index: one MPHF over the distinct sparse k-mers of all
/// genomes, plus a multi-owner CSR (a k-mer may belong to several genomes).
/// Owners are a *multiset* -- a genome appears once per occurrence of the k-mer
/// in its sparse set -- so duplicate k-mers count exactly as the per-genome
/// `get_stats` loop would.
pub struct ScreenIndex {
    pub screen_c: usize,
    pub k: usize,
    mphf: Mphf<u64>,
    /// Per slot: fingerprint of the k-mer, to reject foreign (non-indexed) hashes
    /// that the MPHF would otherwise map to an arbitrary slot.
    fingerprints: Vec<u32>,
    /// CSR row offsets, length `n_slots + 1`.
    owner_offsets: Vec<u32>,
    /// CSR owner genome ids (flat); slot `s` owns `owners[off[s]..off[s+1]]`.
    owners: Vec<u32>,
    /// Per genome: number of sparse k-mers (the `n_kmers` ANI denominator).
    pub sparse_count: Vec<u32>,
}

impl ScreenIndex {
    /// Build from each genome's sparse (`screen_c`) k-mers. Owners are kept as a
    /// multiset so the screen reproduces the per-genome `get_stats` counts
    /// exactly (including any duplicate k-mers within a genome).
    pub fn build(sparse_per_genome: &[Vec<u64>], screen_c: usize, k: usize) -> ScreenIndex {
        let sparse_count: Vec<u32> = sparse_per_genome.iter().map(|v| v.len() as u32).collect();
        let total: usize = sparse_per_genome.iter().map(|v| v.len()).sum();

        // (k-mer, genome) pairs, sorted so equal k-mers form contiguous runs.
        let mut pairs: Vec<(u64, u32)> = Vec::with_capacity(total);
        for (g, v) in sparse_per_genome.iter().enumerate() {
            for &h in v {
                pairs.push((h, g as u32));
            }
        }
        pairs.sort_unstable();

        // Distinct keys for the MPHF.
        let mut keys: Vec<u64> = Vec::new();
        for &(h, _) in &pairs {
            if keys.last() != Some(&h) {
                keys.push(h);
            }
        }
        let mphf = Mphf::new_parallel(MPHF_GAMMA, &keys, Some(0));
        let n_slots = keys.len();

        // Per-slot owner counts -> CSR offsets.
        let mut fingerprints = vec![0u32; n_slots];
        let mut owner_offsets = vec![0u32; n_slots + 1];
        let mut i = 0;
        while i < pairs.len() {
            let h = pairs[i].0;
            let mut j = i;
            while j < pairs.len() && pairs[j].0 == h {
                j += 1;
            }
            let slot = mphf.hash(&h) as usize;
            fingerprints[slot] = sparse_fingerprint(h);
            owner_offsets[slot + 1] = (j - i) as u32; // count, prefix-summed below
            i = j;
        }
        for s in 0..n_slots {
            owner_offsets[s + 1] += owner_offsets[s];
        }

        // Fill owners using a per-slot write cursor.
        let mut owners = vec![0u32; total];
        let mut cursor: Vec<u32> = owner_offsets[..n_slots].to_vec();
        let mut i = 0;
        while i < pairs.len() {
            let h = pairs[i].0;
            let slot = mphf.hash(&h) as usize;
            let mut j = i;
            while j < pairs.len() && pairs[j].0 == h {
                owners[cursor[slot] as usize] = pairs[j].1;
                cursor[slot] += 1;
                j += 1;
            }
            i = j;
        }

        ScreenIndex {
            screen_c,
            k,
            mphf,
            fingerprints,
            owner_offsets,
            owners,
            sparse_count,
        }
    }

    pub fn num_genomes(&self) -> usize {
        self.sparse_count.len()
    }

    /// Single inverted pass over the sample: for each sample k-mer below the
    /// screen threshold with non-zero count, look it up and push its coverage to
    /// every owning genome. Returns `genome -> matched coverage counts`, exactly
    /// the per-genome `covs` a `get_stats(.., None, ..)` screen would collect.
    pub fn gather_hits(&self, sample: &SequencesSketch) -> FxHashMap<u32, Vec<u32>> {
        let thresh = screen_threshold(self.screen_c);
        let mut hits: FxHashMap<u32, Vec<u32>> = FxHashMap::default();
        for (&h, &cnt) in sample.kmer_counts.iter() {
            if h >= thresh || cnt == 0 {
                continue;
            }
            if let Some(slot) = self.mphf.try_hash(&h) {
                let slot = slot as usize;
                if slot < self.fingerprints.len()
                    && self.fingerprints[slot] == sparse_fingerprint(h)
                {
                    let lo = self.owner_offsets[slot] as usize;
                    let hi = self.owner_offsets[slot + 1] as usize;
                    for &g in &self.owners[lo..hi] {
                        hits.entry(g).or_default().push(cnt);
                    }
                }
            }
        }
        hits
    }

    /// Serialize the index into `out` (raw little-endian blocks).
    fn write_to_vec(&self, out: &mut Vec<u8>) -> io::Result<()> {
        let mphf_bytes = bincode::serialize(&self.mphf).map_err(io::Error::other)?;
        write_uvarint(out, mphf_bytes.len() as u64);
        out.extend_from_slice(&mphf_bytes);
        write_uvarint(out, self.fingerprints.len() as u64); // n_slots
        write_uvarint(out, self.owners.len() as u64);
        write_uvarint(out, self.sparse_count.len() as u64); // n_genomes
        for &fp in &self.fingerprints {
            out.extend_from_slice(&fp.to_le_bytes());
        }
        for &o in &self.owner_offsets {
            out.extend_from_slice(&o.to_le_bytes());
        }
        for &o in &self.owners {
            out.extend_from_slice(&o.to_le_bytes());
        }
        for &c in &self.sparse_count {
            out.extend_from_slice(&c.to_le_bytes());
        }
        Ok(())
    }

    /// Parse an index block produced by `write_to_vec`. `screen_c`/`k` come from
    /// the footer (not duplicated in the block).
    fn read(mut r: &[u8], screen_c: usize, k: usize) -> io::Result<ScreenIndex> {
        let mphf_len = read_uvarint(&mut r)? as usize;
        let mut mphf_bytes = vec![0u8; mphf_len];
        r.read_exact(&mut mphf_bytes)?;
        let mphf: Mphf<u64> = bincode::deserialize(&mphf_bytes)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let n_slots = read_uvarint(&mut r)? as usize;
        let n_owners = read_uvarint(&mut r)? as usize;
        let n_genomes = read_uvarint(&mut r)? as usize;

        let read_u32_vec = |r: &mut &[u8], n: usize| -> io::Result<Vec<u32>> {
            let mut v = Vec::with_capacity(n);
            for _ in 0..n {
                let mut buf = [0u8; 4];
                r.read_exact(&mut buf)?;
                v.push(u32::from_le_bytes(buf));
            }
            Ok(v)
        };
        let fingerprints = read_u32_vec(&mut r, n_slots)?;
        let owner_offsets = read_u32_vec(&mut r, n_slots + 1)?;
        let owners = read_u32_vec(&mut r, n_owners)?;
        let sparse_count = read_u32_vec(&mut r, n_genomes)?;

        Ok(ScreenIndex {
            screen_c,
            k,
            mphf,
            fingerprints,
            owner_offsets,
            owners,
            sparse_count,
        })
    }
}

// --- writing ----------------------------------------------------------------

/// Golomb-Rice encode one genome's dense region: its `genome_kmers` block, the
/// pseudotax flag, and (if present) its `pseudotax_tracked_nonused_kmers` block.
/// The inverse of the closure in [`TwoStageDb::decode_dense`].
fn encode_dense_block(gs: &GenomeSketch) -> Vec<u8> {
    let mut block = Vec::new();
    write_hashes(&mut block, &gs.genome_kmers);
    match &gs.pseudotax_tracked_nonused_kmers {
        Some(p) => {
            block.push(1);
            write_hashes(&mut block, p);
        }
        None => block.push(0),
    }
    block
}

/// The sparse (stage-1 screen) subset of a genome's k-mers: those below the
/// `screen_c` FracMinHash threshold. Duplicates are kept, so the pooled index
/// counts a repeated k-mer exactly as a per-genome `get_stats` screen would.
fn sparse_subset(genome_kmers: &[u64], screen_c: usize) -> Vec<u64> {
    let thresh = if screen_c == 0 {
        u64::MAX
    } else {
        u64::MAX / screen_c as u64
    };
    genome_kmers
        .iter()
        .copied()
        .filter(|&h| h < thresh)
        .collect()
}

/// Below this many sparse k-mers, a genome is warned about loudly: it has so few
/// k-mers in total that even the adaptive floor cannot be reached, so it both screens
/// unreliably and drags the whole database's effective screen rate down toward the
/// dense rate (see [`DbBuilder::finish`]).
const SPARSE_WARN_THRESHOLD: usize = 20;

/// The `n` smallest hashes of `hashes` (all of them if it has fewer), as a multiset
/// so repeated k-mers keep their multiplicity. O(n) partial selection.
fn smallest_n(hashes: &[u64], n: usize) -> Vec<u64> {
    if n == 0 {
        return Vec::new();
    }
    if n >= hashes.len() {
        return hashes.to_vec();
    }
    let mut v = hashes.to_vec();
    v.select_nth_unstable(n - 1);
    v.truncate(n);
    v
}

/// A genome's stage-1 sparse set: the FracMinHash subset at `screen_c`, unless that
/// falls short of `min_sparse_kmers`, in which case the `min_sparse_kmers` smallest
/// hashes are taken instead — a denser, genome-specific screen rate.
///
/// Without the floor a small genome is nearly invisible to the screen: at
/// `--screen-c 3000` a 50 kbp plasmid has ~17 sparse k-mers, so a couple of chance
/// matches dominate its screen ANI and it is dropped (or admitted) on noise. This is
/// the same rule, and the same default, as upstream sylph's `--min-sparse-kmers`.
///
/// In every case the result is a *prefix of the genome's dense hashes ordered by
/// value*: all of them, those below a threshold, or the smallest `min_sparse_kmers`.
/// That invariant is what lets [`TwoStageDb::raw_block_and_sparse`] rebuild the exact
/// same set from the stored per-genome count alone, so `db-add` reproduces a
/// from-scratch build without the selection rate having to be stored per genome.
fn select_sparse(genome_kmers: &[u64], screen_c: usize, min_sparse_kmers: usize) -> Vec<u64> {
    let nominal = sparse_subset(genome_kmers, screen_c);
    // `nominal.len() == genome_kmers.len()` means every dense k-mer is already in the
    // sparse set, so there is nothing denser to fall back to.
    if nominal.len() >= min_sparse_kmers || nominal.len() == genome_kmers.len() {
        return nominal;
    }
    let sparse = smallest_n(genome_kmers, min_sparse_kmers);
    debug!(
        "using a denser genome-specific stage-1 screen rate: {} sparse k-mers at -c {} is \
         below --min-sparse-kmers {}, taking the {} smallest of {} dense k-mers instead",
        nominal.len(),
        screen_c,
        min_sparse_kmers,
        sparse.len(),
        genome_kmers.len()
    );
    sparse
}

/// Streaming writer for the two-stage layout. Each genome's dense block is written
/// out as it is added, so the body is never held in RAM; only the per-genome
/// metadata and *sparse* k-mers are accumulated, because the pooled stage-1 index
/// is built over all of them at the end.
///
/// The header records the checksum and the index/footer offsets, none of which are
/// known until everything has been written, so it starts as placeholders that
/// [`DbBuilder::finish`] patches — hence the `Seek` bound.
pub struct DbBuilder<W: Write + Seek> {
    w: crate::checksum::HashingWriter<W>,
    /// Absolute file offset of the next byte to be written.
    pos: u64,
    c: usize,
    k: usize,
    /// Nominal stage-1 screen rate requested by the caller (`--screen-c`). The rate
    /// actually recorded may be denser; see [`DbBuilder::finish`].
    screen_c: usize,
    /// Floor on each genome's sparse k-mer count (`--min-sparse-kmers`).
    min_sparse_kmers: usize,
    genomes: Vec<GenomeMeta>,
    sparse_per_genome: Vec<Vec<u64>>,
}

impl<W: Write + Seek> DbBuilder<W> {
    pub fn new(
        mut w: W,
        c: usize,
        k: usize,
        screen_c: usize,
        min_sparse_kmers: usize,
    ) -> io::Result<DbBuilder<W>> {
        w.write_all(MAGIC)?;
        w.write_all(&[VERSION])?;
        w.write_all(&0u64.to_le_bytes())?; // checksum, patched by finish()
        w.write_all(&0u64.to_le_bytes())?; // index offset, patched by finish()
        w.write_all(&0u64.to_le_bytes())?; // footer offset, patched by finish()
        Ok(DbBuilder {
            // Only the body/index/footer are covered by the checksum: the header is
            // patched afterwards, and every field in it is validated on open anyway.
            w: crate::checksum::HashingWriter::new(w),
            pos: HEADER_LEN,
            c,
            k,
            screen_c,
            min_sparse_kmers: min_sparse_kmers.max(1),
            genomes: Vec::new(),
            sparse_per_genome: Vec::new(),
        })
    }

    /// Append a genome whose dense block is *already* encoded — copied verbatim out
    /// of an existing database, so it is never decoded and re-encoded. `sparse` is
    /// its stage-1 subset at this builder's `screen_c` and `meta` its metadata; the
    /// `dense_offset` in `meta` is ignored and replaced with the offset in this file.
    pub fn push_encoded(
        &mut self,
        block: &[u8],
        sparse: Vec<u64>,
        meta: &GenomeMeta,
    ) -> io::Result<()> {
        self.w.write_all(block)?;
        self.genomes.push(GenomeMeta {
            dense_offset: self.pos,
            ..meta.clone()
        });
        self.pos += block.len() as u64;
        self.sparse_per_genome.push(sparse);
        Ok(())
    }

    /// Append a genome sketch, encoding its dense block and deriving its stage-1
    /// sparse subset.
    pub fn push_sketch(&mut self, gs: &GenomeSketch) -> io::Result<()> {
        let block = encode_dense_block(gs);
        let sparse = select_sparse(&gs.genome_kmers, self.screen_c, self.min_sparse_kmers);
        if sparse.len() < SPARSE_WARN_THRESHOLD.min(self.min_sparse_kmers) {
            warn!(
                "genome '{}' (file {}) has only {} k-mers in total; its whole k-mer set is used \
                 as its stage-1 screen entry, and because it is (one of) the densest genome(s) \
                 in this database it drags the WHOLE database's stage-1 screen rate down toward \
                 the dense -c {}, making screening slower for every sample. Detection \
                 reliability at this size is inherently poor -- consider keeping tiny \
                 genomes/contigs/plasmids in a plain .syldb instead, or check that this \
                 genome/contig was sketched as intended.",
                gs.first_contig_name,
                gs.file_name,
                sparse.len(),
                self.c
            );
        }
        self.push_encoded(&block, sparse, &meta_of(gs, 0))
    }

    pub fn num_genomes(&self) -> usize {
        self.genomes.len()
    }

    /// Coarsest screen rate whose FracMinHash threshold still admits every sparse
    /// k-mer pushed so far, capped at the requested `screen_c`.
    ///
    /// Safe by the integer identity `floor(a / floor(a / b)) >= b`: with
    /// `c = u64::MAX / (max_hash + 1)`, `screen_threshold(c) = u64::MAX / c >=
    /// max_hash + 1 > max_hash`, so every stored k-mer passes `gather_hits`'
    /// early-exit filter.
    fn effective_screen_c(&self) -> usize {
        let max_sparse = self
            .sparse_per_genome
            .iter()
            .flat_map(|v| v.iter().copied())
            .max();
        match max_sparse {
            Some(h) => self
                .screen_c
                .min((u64::MAX / h.saturating_add(1)).max(1) as usize),
            None => self.screen_c,
        }
    }

    /// Build the pooled stage-1 index, write it plus the footer, and patch the
    /// header with the checksum and section offsets.
    ///
    /// The screen rate recorded in the file is the *effective* one: `gather_hits`
    /// skips sample k-mers at or above `screen_threshold(screen_c)`, so a genome that
    /// had to use a denser genome-specific rate (see [`select_sparse`]) would have its
    /// densest stored k-mers unmatchable if the nominal rate were recorded. The
    /// effective rate is therefore the coarsest rate whose threshold still admits
    /// every stored sparse k-mer, and never coarser than the requested one.
    pub fn finish(mut self) -> io::Result<()> {
        let index_offset = self.pos;
        let effective_screen_c = self.effective_screen_c();
        if effective_screen_c != self.screen_c {
            info!(
                "effective stage-1 screen -c adjusted from {} to {} to keep small genomes' \
                 screen k-mers matchable (see warnings above)",
                self.screen_c, effective_screen_c
            );
        }
        let screen_index = ScreenIndex::build(&self.sparse_per_genome, effective_screen_c, self.k);
        self.sparse_per_genome = Vec::new();
        let mut index_block: Vec<u8> = Vec::new();
        screen_index.write_to_vec(&mut index_block)?;
        self.w.write_all(&index_block)?;
        let footer_offset = index_offset + index_block.len() as u64;
        drop(index_block);

        let footer = Footer {
            c: self.c,
            k: self.k,
            screen_c: effective_screen_c,
            genomes: self.genomes,
        };
        let footer_bytes = bincode::serialize(&footer).map_err(io::Error::other)?;
        self.w.write_all(&footer_bytes)?;

        // Everything after the header is checksummed, so `inspect` can detect a
        // truncated or bit-rotted database that a seeking reader would otherwise
        // decode into a plausible-looking but wrong sketch.
        let (mut w, checksum) = self.w.finish();
        w.flush()?;
        w.seek(SeekFrom::Start(5))?;
        w.write_all(&checksum.to_le_bytes())?;
        w.write_all(&index_offset.to_le_bytes())?;
        w.write_all(&footer_offset.to_le_bytes())?;
        w.flush()?;
        Ok(())
    }
}

/// Per-genome footer metadata for a sketch, at a given dense-block offset.
fn meta_of(gs: &GenomeSketch, dense_offset: u64) -> GenomeMeta {
    GenomeMeta {
        file_name: gs.file_name.clone(),
        first_contig_name: gs.first_contig_name.clone(),
        gn_size: gs.gn_size,
        min_spacing: gs.min_spacing,
        has_pseudotax: gs.pseudotax_tracked_nonused_kmers.is_some(),
        dense_offset,
    }
}

/// Re-pack genome sketches into the two-stage seekable layout and write to `w`.
/// `screen_c` is the (coarser) stage-1 subsampling rate; it must be `>= c`.
/// `min_sparse_kmers` is the per-genome sparse floor (see [`select_sparse`]) and must
/// be `>= 1` -- 0 would let a genome's screen entry be empty, making it invisible to
/// the stage-1 screen forever. Dense blocks are Golomb-Rice coded.
pub fn write_two_stage_db<W: Write + Seek>(
    w: W,
    sketches: &[GenomeSketch],
    screen_c: usize,
    min_sparse_kmers: usize,
) -> io::Result<()> {
    let c = sketches.first().map(|s| s.c).unwrap_or(0);
    let k = sketches.first().map(|s| s.k).unwrap_or(0);
    let mut builder = DbBuilder::new(w, c, k, screen_c, min_sparse_kmers)?;
    for gs in sketches {
        builder.push_sketch(gs)?;
    }
    builder.finish()
}

// --- opened database --------------------------------------------------------

/// Backing store for the dense blocks.
///   * `File`  - positional `read_at` (pread) of just the requested block bytes.
///     No shared cursor, so concurrent reads from any number of threads need no
///     lock; only the touched block bytes (plus reclaimable OS page cache) cost
///     memory, so RSS stays low. This is the path used for `.syl2db` files.
///   * `Owned` - whole file in memory; for in-memory readers / tests.
enum DenseData {
    File(File),
    Owned(Vec<u8>),
}

impl DenseData {
    /// Whole-file bytes; only valid for the in-memory backing (used to parse the
    /// header/footer). The `File` backing is read positionally via `with_block`.
    #[inline]
    fn bytes(&self) -> &[u8] {
        match self {
            DenseData::Owned(v) => &v[..],
            DenseData::File(_) => unreachable!("File-backed db is read via with_block"),
        }
    }
}

/// An opened two-stage database. Construction loads only the stage-1 sparse
/// index (the bincoded footer); per-genome dense blocks are decoded on demand,
/// in parallel (each decode positionally reads its own block, no shared cursor
/// or lock).
pub struct TwoStageDb {
    pub c: usize,
    pub k: usize,
    pub screen_c: usize,
    /// XXH64 of the file after the header, as recorded when it was written. Checked
    /// by [`TwoStageDb::verify_checksum`], not on open: validating it costs a full
    /// read of a file that profiling otherwise only touches a few blocks of. `None`
    /// for an upstream sylph (v2) database, which carries no checksum.
    checksum: Option<u64>,
    /// Header length of the file as opened; the checksum covers everything after it.
    header_len: u64,
    /// File offset where the dense-block region ends (start of the screen index);
    /// used to bound the last genome's block for positional reads.
    index_offset: u64,
    genomes: Vec<GenomeMeta>,
    /// Pooled stage-1 screen index (Path B). Replaces the per-genome sparse
    /// sketches; querying a sample against it yields the contained genomes.
    pub screen_index: ScreenIndex,
    data: DenseData,
    cache: Mutex<FxHashMap<u32, Arc<GenomeSketch>>>,
}

/// A parsed `.syl2db` header. `checksum` is absent for upstream sylph's version 2,
/// whose header has no checksum field.
struct Header {
    len: u64,
    checksum: Option<u64>,
    index_offset: u64,
    footer_offset: u64,
}

/// Parse the magic + version header of either readable version.
///
/// Version 2 is upstream sylph's: the same layout minus the 8-byte whole-file
/// checksum, so its header is 8 bytes shorter and its offsets sit 8 bytes earlier.
/// Everything after the header -- dense blocks, screen index, footer -- is identical
/// between the two, so a sylph-built database is read here directly (and vice versa
/// once sylph learns to skip the extra 8 bytes); only the on-demand integrity check
/// is unavailable for a v2 file.
fn parse_header(hdr: &[u8]) -> io::Result<Header> {
    if hdr.len() < HEADER_LEN_V2 as usize || &hdr[0..4] != MAGIC {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "not a sylph two-stage database",
        ));
    }
    let version = hdr[4];
    let (len, checksum_bytes) = match version {
        VERSION => (HEADER_LEN, Some(5..13)),
        VERSION_NO_CHECKSUM => (HEADER_LEN_V2, None),
        other => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "two-stage database is version {}, but only versions {} and {} are readable; rebuild it with db-convert",
                    other, VERSION_NO_CHECKSUM, VERSION
                ),
            ))
        }
    };
    if (hdr.len() as u64) < len {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "two-stage database header is truncated",
        ));
    }
    let checksum = checksum_bytes.map(|r| u64::from_le_bytes(hdr[r].try_into().unwrap()));
    let offsets_at = len as usize - 16;
    let index_offset = u64::from_le_bytes(hdr[offsets_at..offsets_at + 8].try_into().unwrap());
    let footer_offset =
        u64::from_le_bytes(hdr[offsets_at + 8..offsets_at + 16].try_into().unwrap());
    Ok(Header {
        len,
        checksum,
        index_offset,
        footer_offset,
    })
}

/// Assemble a `TwoStageDb` from its parsed footer + screen index + backing store.
fn build_db(
    footer: Footer,
    header: &Header,
    screen_index: ScreenIndex,
    data: DenseData,
) -> TwoStageDb {
    TwoStageDb {
        c: footer.c,
        k: footer.k,
        screen_c: footer.screen_c,
        checksum: header.checksum,
        header_len: header.len,
        index_offset: header.index_offset,
        genomes: footer.genomes,
        screen_index,
        data,
        cache: Mutex::new(FxHashMap::default()),
    }
}

/// Parse the header + index + footer of a `.syl2db` already resident in `data`.
fn from_bytes(data: DenseData) -> io::Result<TwoStageDb> {
    let bytes = data.bytes();
    let header = parse_header(bytes)?;
    if header.index_offset > header.footer_offset || header.footer_offset as usize > bytes.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "two-stage database offsets out of range",
        ));
    }
    let footer: Footer = bincode::deserialize(&bytes[header.footer_offset as usize..])
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    let screen_index = ScreenIndex::read(
        &bytes[header.index_offset as usize..header.footer_offset as usize],
        footer.screen_c,
        footer.k,
    )?;
    Ok(build_db(footer, &header, screen_index, data))
}

/// Open a `.syl2db` from an in-memory reader (reads it all into memory).
pub fn open<R: Read>(mut r: R) -> io::Result<TwoStageDb> {
    let mut v = Vec::new();
    r.read_to_end(&mut v)?;
    from_bytes(DenseData::Owned(v))
}

impl TwoStageDb {
    pub fn len(&self) -> usize {
        self.genomes.len()
    }
    pub fn is_empty(&self) -> bool {
        self.genomes.is_empty()
    }

    /// Source file name of genome `g` (for `--screen-dump` and diagnostics).
    pub fn genome_file_name(&self, g: u32) -> &str {
        &self.genomes[g as usize].file_name
    }

    /// `(file_name, gn_size)` for every genome. Used by `profile --apply-unknown`
    /// to look up genome sizes without decoding any dense block.
    pub fn genome_sizes(&self) -> Vec<(&str, usize)> {
        self.genomes
            .iter()
            .map(|g| (g.file_name.as_str(), g.gn_size))
            .collect()
    }

    /// End offset of genome `g`'s region (start of the next genome's block, or
    /// the screen index for the last genome). Genomes are stored in ascending
    /// offset, and the index block immediately follows the last dense block.
    fn block_end(&self, g: u32) -> u64 {
        let gi = g as usize;
        if gi + 1 < self.genomes.len() {
            self.genomes[gi + 1].dense_offset
        } else {
            self.index_offset
        }
    }

    /// Run `f` on genome `g`'s whole on-disk region (`genome_kmers` block, the
    /// pseudotax flag, and the optional pseudotax block). For the file backing
    /// this is one positional read of just that region; for the owned backing it
    /// is a zero-copy slice. No shared cursor, so it is safe to call concurrently
    /// from many threads.
    fn with_block<T>(&self, g: u32, f: impl FnOnce(&[u8]) -> io::Result<T>) -> io::Result<T> {
        let start = self.genomes[g as usize].dense_offset as usize;
        match &self.data {
            DenseData::Owned(v) => f(&v[start..]),
            DenseData::File(file) => {
                let end = self.block_end(g) as usize;
                let mut buf = vec![0u8; end - start];
                file.read_exact_at(&mut buf, start as u64)?;
                f(&buf)
            }
        }
    }

    /// Decode genome `g`'s full dense `GenomeSketch` without touching the cache.
    /// Concurrent calls from different threads do not contend on any shared
    /// cursor (each does its own positional read). The two-stage pass-1 uses this
    /// to decode each survivor into a short-lived sketch, probe it, and drop it
    /// unless it passes -- so the discarded majority is never cached.
    pub fn decode_dense(&self, g: u32) -> io::Result<GenomeSketch> {
        let meta = &self.genomes[g as usize];
        let (genome_kmers, pseudotax) = self.with_block(g, |bytes| {
            let mut cur = bytes;
            let gk = read_hashes(&mut cur)?;
            let mut flag = [0u8; 1];
            cur.read_exact(&mut flag)?;
            let pt = if flag[0] != 0 {
                Some(read_hashes(&mut cur)?)
            } else {
                None
            };
            Ok((gk, pt))
        })?;
        Ok(GenomeSketch {
            genome_kmers,
            pseudotax_tracked_nonused_kmers: pseudotax,
            file_name: meta.file_name.clone(),
            first_contig_name: meta.first_contig_name.clone(),
            c: self.c,
            k: self.k,
            gn_size: meta.gn_size,
            min_spacing: meta.min_spacing,
        })
    }

    /// Metadata of genome `g`, as stored in the footer.
    pub fn genome_meta(&self, g: u32) -> &GenomeMeta {
        &self.genomes[g as usize]
    }

    /// Everything needed to copy genome `g` into another database: its dense region
    /// exactly as it sits on disk (so it is never re-encoded), plus its stage-1
    /// sparse k-mers.
    ///
    /// The sparse set has to be recovered by decoding the block, because the pooled
    /// `ScreenIndex` stores only the inverted k-mer → owners mapping and cannot hand
    /// the per-genome sparse sets back. One positional read serves both.
    ///
    /// It is recovered from the stored *count* rather than by re-subsampling at
    /// `screen_c`: a genome's sparse set is always the smallest `n` of its dense
    /// hashes (see [`select_sparse`]), and `n` is in the screen index, so this
    /// reproduces the original selection exactly — including for a small genome that
    /// was selected at a denser genome-specific rate, which re-subsampling at the
    /// database's screen rate would silently change.
    pub fn raw_block_and_sparse(&self, g: u32) -> io::Result<(Vec<u8>, Vec<u64>)> {
        let n_sparse = self.screen_index.sparse_count[g as usize] as usize;
        self.with_block(g, |bytes| {
            // The owned backing hands back everything from the block's start, so trim
            // to this genome's region; the file backing already reads exactly it.
            let len = match &self.data {
                DenseData::Owned(_) => {
                    (self.block_end(g) - self.genomes[g as usize].dense_offset) as usize
                }
                DenseData::File(_) => bytes.len(),
            };
            let block = &bytes[..len];
            let mut cur = block;
            let gk = read_hashes(&mut cur)?;
            Ok((block.to_vec(), smallest_n(&gk, n_sparse)))
        })
    }

    /// Re-hash the whole file and compare against the checksum in its header. This
    /// reads every byte, so it is on-demand (`weebill inspect`) rather than part of
    /// opening the database.
    ///
    /// `Ok(false)` means the file carries no checksum to check (an upstream sylph v2
    /// database), as opposed to `Ok(true)` for a verified one.
    pub fn verify_checksum(&self) -> io::Result<bool> {
        let Some(expected) = self.checksum else {
            return Ok(false);
        };
        let got = match &self.data {
            DenseData::Owned(v) => crate::checksum::hash_reader(&v[self.header_len as usize..])?,
            DenseData::File(file) => {
                let mut r = BufReader::with_capacity(1 << 20, file.try_clone()?);
                r.seek(SeekFrom::Start(self.header_len))?;
                crate::checksum::hash_reader(r)?
            }
        };
        if got != expected {
            return Err(crate::checksum::mismatch(
                "the two-stage database",
                expected,
                got,
            ));
        }
        Ok(true)
    }

    /// Decode genome `g`'s full dense `GenomeSketch`, caching it across calls.
    pub fn load_dense(&self, g: u32) -> io::Result<Arc<GenomeSketch>> {
        if let Some(a) = self.cache.lock().unwrap().get(&g) {
            return Ok(a.clone());
        }
        let sketch = Arc::new(self.decode_dense(g)?);
        self.cache.lock().unwrap().insert(g, sketch.clone());
        Ok(sketch)
    }
}

/// Open a `.syl2db` file from a path. Only the header + footer (the stage-1
/// sparse index) are read up front; dense blocks are read positionally on demand
/// during profiling, so opening is cheap and RSS stays low.
pub fn open_file(path: &str) -> io::Result<TwoStageDb> {
    let file = File::open(path)?;
    let mut hdr = [0u8; HEADER_LEN as usize];
    // A v2 (upstream sylph) database is 8 bytes shorter in the header, and could in
    // principle be a file of only HEADER_LEN_V2 bytes, so a short read is not fatal
    // here; `parse_header` rejects anything genuinely too short.
    let read = read_at_most(&file, &mut hdr, 0)?;
    let header = parse_header(&hdr[..read])?;
    let flen = file.metadata()?.len();
    if header.index_offset > header.footer_offset || header.footer_offset > flen {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "two-stage database offsets out of range",
        ));
    }
    let mut fbytes = vec![0u8; (flen - header.footer_offset) as usize];
    file.read_exact_at(&mut fbytes, header.footer_offset)?;
    let footer: Footer =
        bincode::deserialize(&fbytes).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    let mut ibytes = vec![0u8; (header.footer_offset - header.index_offset) as usize];
    file.read_exact_at(&mut ibytes, header.index_offset)?;
    let screen_index = ScreenIndex::read(&ibytes, footer.screen_c, footer.k)?;
    Ok(build_db(
        footer,
        &header,
        screen_index,
        DenseData::File(file),
    ))
}

/// Positional read of up to `buf.len()` bytes, returning how many were read.
fn read_at_most(file: &File, buf: &mut [u8], offset: u64) -> io::Result<usize> {
    let mut done = 0;
    while done < buf.len() {
        match file.read_at(&mut buf[done..], offset + done as u64)? {
            0 => break,
            n => done += n,
        }
    }
    Ok(done)
}

// --- CLI handler ------------------------------------------------------------

/// Load every genome sketch from a database, in either the legacy bincode
/// (`.syldb`) or compressed (`.syldbc`) encoding — detected by content, as
/// everywhere else, not by extension.
fn load_genome_sketches(path: &str) -> Vec<GenomeSketch> {
    let file = File::open(path).unwrap_or_else(|_| panic!("Could not open {}", path));
    let mut reader = BufReader::with_capacity(10_000_000, file);
    if crate::compress::peek_is_compressed(&mut reader).unwrap_or(false) {
        crate::compress::read_genome_sketches_compressed(&mut reader)
            .unwrap_or_else(|e| panic!("{} is not a valid compressed database sketch: {}", path, e))
    } else {
        bincode::deserialize_from(&mut reader)
            .unwrap_or_else(|_| panic!("{} is not a valid database sketch (.syldb/.syldbc)", path))
    }
}

pub fn run_db_convert(args: DbConvertArgs) {
    let level = if args.trace {
        log::LevelFilter::Trace
    } else if args.debug {
        log::LevelFilter::Debug
    } else {
        log::LevelFilter::Info
    };
    simple_logger::SimpleLogger::new()
        .with_level(level)
        .init()
        .unwrap();
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .ok();

    if args.files.is_empty() {
        error!("No genome database sketches (*.syldb) supplied; exiting");
        std::process::exit(1);
    }

    let mut sketches: Vec<GenomeSketch> = Vec::new();
    for f in &args.files {
        info!("Loading genome sketches from {}", f);
        sketches.extend(load_genome_sketches(f));
    }
    if sketches.is_empty() {
        error!("No genome sketches found in input; exiting");
        std::process::exit(1);
    }

    let c = sketches[0].c;
    let k = sketches[0].k;
    for s in &sketches {
        if s.c != c || s.k != k {
            error!("Input sketches have inconsistent -c/-k; exiting");
            std::process::exit(1);
        }
    }
    if sketches
        .iter()
        .any(|s| s.pseudotax_tracked_nonused_kmers.is_none())
    {
        error!(
            "Some input genomes were sketched with --disable-profiling (no profiling k-mers). \
             A two-stage database is for `profile`; re-sketch without --disable-profiling. Exiting"
        );
        std::process::exit(1);
    }
    if args.screen_c < c {
        error!(
            "--screen-c ({}) must be >= the database -c ({}); the screen can only be made sparser, never denser. Exiting",
            args.screen_c, c
        );
        std::process::exit(1);
    }
    if args.min_sparse_kmers < 1 {
        error!(
            "--min-sparse-kmers must be >= 1: a genome with an empty stage-1 screen entry \
             could never pass the screen. Exiting"
        );
        std::process::exit(1);
    }

    let out = if args.output.ends_with(TWO_STAGE_DB_SUFFIX) {
        args.output.clone()
    } else {
        format!("{}{}", args.output, TWO_STAGE_DB_SUFFIX)
    };
    if let Some(parent) = Path::new(&out).parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent).ok();
        }
    }
    info!(
        "Converting {} genomes (dense -c {}, stage-1 screen -c {}) -> {}",
        sketches.len(),
        c,
        args.screen_c,
        out
    );
    let w =
        BufWriter::new(File::create(&out).unwrap_or_else(|_| panic!("Could not create {}", out)));
    write_two_stage_db(w, &sketches, args.screen_c, args.min_sparse_kmers)
        .unwrap_or_else(|e| panic!("Failed to write {}: {}", out, e));
    info!("Wrote two-stage database to {}", out);
}

/// Resolve the output path for `db-add`, appending the suffix if absent.
fn suffixed_output(output: &str) -> String {
    if output.ends_with(TWO_STAGE_DB_SUFFIX) {
        output.to_string()
    } else {
        format!("{}{}", output, TWO_STAGE_DB_SUFFIX)
    }
}

pub fn run_db_add(args: DbAddArgs) {
    let level = if args.trace {
        log::LevelFilter::Trace
    } else if args.debug {
        log::LevelFilter::Debug
    } else {
        log::LevelFilter::Info
    };
    simple_logger::SimpleLogger::new()
        .with_level(level)
        .init()
        .unwrap();
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .ok();

    if args.files.is_empty() && args.genomes.is_empty() {
        error!(
            "No genomes to add: give database sketches (*.syldb/*.syldbc) as arguments and/or FASTA files with -g. Exiting"
        );
        std::process::exit(1);
    }

    // Open the existing database first (loading only its stage-1 index + footer), so
    // its -c/-k/--screen-c dictate what the new genomes must match and a bad
    // database is reported before anything is sketched or created.
    info!("Opening two-stage database {}", args.database);
    let db = open_file(&args.database).unwrap_or_else(|e| {
        error!("{} is not a valid two-stage database: {}", args.database, e);
        std::process::exit(1);
    });
    if !args.no_verify {
        info!("Verifying {} before growing it...", args.database);
        let verified = db.verify_checksum().unwrap_or_else(|e| {
            error!("{}. Refusing to grow a corrupt database; exiting", e);
            std::process::exit(1);
        });
        if !verified {
            warn!(
                "{} was written by upstream sylph and carries no whole-file checksum, so it \
                 cannot be verified; growing it regardless. The database written here will \
                 carry one.",
                args.database
            );
        }
    }
    info!(
        "{} holds {} genomes (dense -c {}, -k {}, stage-1 screen -c {})",
        args.database,
        db.len(),
        db.c,
        db.k,
        db.screen_c
    );

    // Collect the new genomes: pre-sketched databases, plus any FASTAs sketched
    // here at the existing database's -c/-k (so they are directly comparable).
    let mut new_sketches: Vec<GenomeSketch> = Vec::new();
    for f in &args.files {
        info!("Loading genome sketches from {}", f);
        new_sketches.extend(load_genome_sketches(f));
    }
    if !args.genomes.is_empty() {
        // Default the spacing to whatever the existing genomes used: a mismatch
        // changes which k-mers a genome contributes, so the added genomes would be
        // sketched on subtly different terms from the ones they are profiled against.
        let min_spacing = args.min_spacing_kmer.unwrap_or_else(|| {
            if db.is_empty() {
                30
            } else {
                db.genome_meta(0).min_spacing
            }
        });
        info!(
            "Sketching {} genome fasta(s) at -c {} -k {} --min-spacing {} (matching {})",
            args.genomes.len(),
            db.c,
            db.k,
            min_spacing,
            args.database
        );
        let sketched: Vec<GenomeSketch> = args
            .genomes
            .par_iter()
            .filter_map(|g| crate::sketch::sketch_genome(db.c, db.k, g, min_spacing, true))
            .collect();
        if sketched.len() < args.genomes.len() {
            warn!(
                "Only {} of {} genome fasta(s) could be sketched; the rest were skipped",
                sketched.len(),
                args.genomes.len()
            );
        }
        new_sketches.extend(sketched);
    }
    if new_sketches.is_empty() {
        error!("No genome sketches found in the inputs to add; exiting");
        std::process::exit(1);
    }

    // The new genomes must match the database exactly: FracMinHash lets a sketch be
    // made sparser but never denser, and a -c/-k mismatch would silently make the
    // added genomes' k-mers incomparable to the existing ones.
    for s in &new_sketches {
        if s.c != db.c || s.k != db.k {
            error!(
                "Genome '{}' was sketched with -c {} -k {}, but {} is -c {} -k {}; re-sketch the added genomes to match. Exiting",
                s.file_name, s.c, s.k, args.database, db.c, db.k
            );
            std::process::exit(1);
        }
        if s.pseudotax_tracked_nonused_kmers.is_none() {
            error!(
                "Genome '{}' was sketched with --disable-profiling (no profiling k-mers). A two-stage database is for `profile`; re-sketch without --disable-profiling. Exiting",
                s.file_name
            );
            std::process::exit(1);
        }
    }

    // A .syl2db identifies genomes by `file_name` (that is all a profile TSV
    // reports, and all the dense-sketch cache and --apply-unknown can key on), so a
    // duplicate name would make two genomes indistinguishable downstream.
    let mut existing: FxHashSet<&str> = FxHashSet::default();
    for g in 0..db.len() as u32 {
        existing.insert(db.genome_file_name(g));
    }
    let mut kept: Vec<GenomeSketch> = Vec::with_capacity(new_sketches.len());
    let mut seen_new: FxHashSet<String> = FxHashSet::default();
    let mut skipped = 0usize;
    for s in new_sketches {
        if existing.contains(s.file_name.as_str()) {
            if args.skip_existing {
                skipped += 1;
                continue;
            }
            error!(
                "Genome '{}' is already in {}. Pass --skip-existing to ignore genomes already present. Exiting",
                s.file_name, args.database
            );
            std::process::exit(1);
        }
        if !seen_new.insert(s.file_name.clone()) {
            error!(
                "Genome '{}' appears more than once in the genomes to add; each genome must be distinct. Exiting",
                s.file_name
            );
            std::process::exit(1);
        }
        kept.push(s);
    }
    if skipped > 0 {
        info!(
            "Skipping {} genome(s) already present in {}",
            skipped, args.database
        );
    }
    // Only reachable under --skip-existing (a duplicate is otherwise fatal above),
    // where the point is idempotence: re-running the same add is a no-op, not a
    // failure. An in-place run has nothing to do; an `-o` run still owes the caller
    // the output file, so it falls through and writes the copy.
    let new_sketches = kept;
    if new_sketches.is_empty() {
        if args.output.is_none() {
            info!(
                "Every genome given is already in {}; leaving it unchanged.",
                args.database
            );
            return;
        }
        info!(
            "Every genome given is already in {}; the output will be an unchanged copy of it.",
            args.database
        );
    }

    let out = match &args.output {
        Some(o) => suffixed_output(o),
        // Default is in-place: write a sibling temp file and rename over the
        // original only once it is complete, so an interrupted run cannot leave a
        // half-written database in place of a good one.
        None => args.database.clone(),
    };
    if let Some(parent) = Path::new(&out).parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent).ok();
        }
    }
    // Never truncate an existing output directly. Besides making replacement
    // atomic, this protects the source when `-o` is a hard link or symlink to it:
    // path canonicalization alone cannot identify hard links.
    let replace_existing = Path::new(&out).exists();
    let tmp = format!("{}.tmp{}", out, std::process::id());
    let write_path = if replace_existing {
        tmp.clone()
    } else {
        out.clone()
    };

    info!(
        "Adding {} genome(s) to {} existing -> {} ({} genomes total)",
        new_sketches.len(),
        db.len(),
        out,
        db.len() + new_sketches.len()
    );

    // Copy every existing dense block through verbatim (Golomb-Rice re-encoding
    // would be pure waste), then append the new ones. The stage-1 index cannot be
    // extended in place -- its MPHF is built over the pooled key set -- so it is
    // rebuilt from every genome's sparse k-mers, which is why the copied blocks are
    // still decoded on the way past.
    let result = (|| -> io::Result<()> {
        let replacement_permissions = if replace_existing {
            Some(std::fs::metadata(&out)?.permissions())
        } else {
            None
        };
        let file = File::create(&write_path)?;
        if let Some(permissions) = replacement_permissions {
            file.set_permissions(permissions)?;
        }
        let w = BufWriter::with_capacity(1 << 20, file);
        // The existing genomes' sparse sets are reproduced from their stored counts, so
        // the floor only applies to the genomes being added; they are selected at the
        // database's recorded screen rate, which is what the existing genomes were
        // (nominally) selected at.
        let mut builder = DbBuilder::new(w, db.c, db.k, db.screen_c, args.min_sparse_kmers.max(1))?;
        // Chunked so the reads/decodes of a chunk run in parallel while the writer
        // stays sequential (block order defines the footer offsets). Keeping at most
        // one decoded block per Rayon worker bounds transient dense-block memory by
        // the requested concurrency rather than an arbitrary genome count.
        let copy_chunk = rayon::current_num_threads().max(1);
        for chunk_start in (0..db.len()).step_by(copy_chunk) {
            let chunk_end = (chunk_start + copy_chunk).min(db.len());
            let decoded: Vec<io::Result<(Vec<u8>, Vec<u64>)>> = (chunk_start..chunk_end)
                .into_par_iter()
                .map(|g| db.raw_block_and_sparse(g as u32))
                .collect();
            for (i, d) in decoded.into_iter().enumerate() {
                let (block, sparse) = d?;
                let g = (chunk_start + i) as u32;
                builder.push_encoded(&block, sparse, db.genome_meta(g))?;
            }
            info!("Copied {}/{} existing genomes", chunk_end, db.len());
        }
        for gs in &new_sketches {
            builder.push_sketch(gs)?;
        }
        builder.finish()
    })();
    if let Err(e) = result {
        let _ = std::fs::remove_file(&write_path);
        error!("Failed to write {}: {}. Exiting", write_path, e);
        std::process::exit(1);
    }
    if replace_existing {
        if let Err(e) = std::fs::rename(&tmp, &out) {
            let _ = std::fs::remove_file(&tmp);
            error!("Could not replace {} with the grown database: {}", out, e);
            std::process::exit(1);
        }
    }
    info!(
        "Wrote two-stage database with {} genomes to {}",
        db.len() + new_sketches.len(),
        out
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    fn roundtrip_hashes(input: &[u64]) {
        let mut expected = input.to_vec();
        expected.sort_unstable();
        let mut buf = Vec::new();
        write_hashes(&mut buf, input);
        let mut r = &buf[..];
        assert_eq!(read_hashes(&mut r).unwrap(), expected);
        assert!(r.is_empty(), "read_hashes left trailing bytes");
    }

    #[test]
    fn hashes_roundtrip_various() {
        roundtrip_hashes(&[]);
        roundtrip_hashes(&[0]);
        roundtrip_hashes(&[42]);
        roundtrip_hashes(&[5, 5, 5]); // duplicates -> zero gaps
        roundtrip_hashes(&[u64::MAX, 0, 1, u64::MAX / 2]);
        roundtrip_hashes(&[10, 9, 8, 7, 6, 5, 4, 3, 2, 1]);
        // many uniformly-spread hashes (the realistic FracMinHash case)
        let mut v = Vec::new();
        let mut x = 0xdead_beef_u64;
        for _ in 0..5000 {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            v.push(x >> 8); // keep them in a sub-range like a fracminhash threshold
        }
        roundtrip_hashes(&v);
    }

    /// `write_two_stage_db` patches its header at the end and so needs `Seek`; a
    /// bare `Vec<u8>` is not seekable, hence the cursor.
    fn write_db_to_vec(sketches: &[GenomeSketch], screen_c: usize) -> Vec<u8> {
        write_db_to_vec_min_sparse(sketches, screen_c, 1)
    }

    fn write_db_to_vec_min_sparse(
        sketches: &[GenomeSketch],
        screen_c: usize,
        min_sparse_kmers: usize,
    ) -> Vec<u8> {
        let mut cur = std::io::Cursor::new(Vec::new());
        write_two_stage_db(&mut cur, sketches, screen_c, min_sparse_kmers).unwrap();
        cur.into_inner()
    }

    fn gsketch(name: &str, kmers: Vec<u64>, pt: Option<Vec<u64>>) -> GenomeSketch {
        GenomeSketch {
            genome_kmers: kmers,
            pseudotax_tracked_nonused_kmers: pt,
            file_name: name.to_string(),
            first_contig_name: format!("{}_c1", name),
            c: 50,
            k: 31,
            gn_size: 12345,
            min_spacing: 30,
        }
    }

    #[test]
    fn db_write_open_load_roundtrip() {
        // screen_c = 200 (coarser than c = 50): the sparse subset keeps hashes
        // below u64::MAX/200.
        let thresh = u64::MAX / 200;
        let g0_kmers: Vec<u64> = vec![1, 2, 3, thresh - 1, thresh + 10, thresh * 3, 9_000_000_000];
        let g1_kmers: Vec<u64> = vec![7, thresh + 1, thresh * 2, 123_456_789_000];
        let sketches = vec![
            gsketch("g0.fa", g0_kmers.clone(), Some(vec![100, 200, 300])),
            gsketch("g1.fa", g1_kmers.clone(), Some(vec![])),
        ];

        let db = open(std::io::Cursor::new(write_db_to_vec(&sketches, 200))).unwrap();

        assert_eq!(db.c, 50);
        assert_eq!(db.k, 31);
        assert_eq!(db.screen_c, 200);
        assert_eq!(db.len(), 2);

        // stage-1 index: per-genome sparse count = fracminhash subset size at screen_c
        let expect_sparse = |ks: &[u64]| -> Vec<u64> {
            let mut v: Vec<u64> = ks.iter().copied().filter(|&h| h < thresh).collect();
            v.sort_unstable();
            v
        };
        assert_eq!(db.screen_c, 200);
        assert_eq!(
            db.screen_index.sparse_count[0] as usize,
            expect_sparse(&g0_kmers).len()
        );
        assert_eq!(
            db.screen_index.sparse_count[1] as usize,
            expect_sparse(&g1_kmers).len()
        );

        // stage-2 dense block reconstructs the exact genome k-mers + pseudotax
        let d0 = db.load_dense(0).unwrap();
        let mut got = d0.genome_kmers.clone();
        got.sort_unstable();
        let mut exp = g0_kmers.clone();
        exp.sort_unstable();
        assert_eq!(got, exp);
        assert_eq!(d0.c, 50);
        assert_eq!(d0.k, 31);
        assert_eq!(d0.gn_size, 12345);
        assert_eq!(d0.file_name, "g0.fa");
        assert_eq!(
            d0.pseudotax_tracked_nonused_kmers,
            Some(vec![100, 200, 300])
        );

        let d1 = db.load_dense(1).unwrap();
        let mut got1 = d1.genome_kmers.clone();
        got1.sort_unstable();
        let mut exp1 = g1_kmers.clone();
        exp1.sort_unstable();
        assert_eq!(got1, exp1);
        assert_eq!(d1.pseudotax_tracked_nonused_kmers, Some(vec![]));

        // second load hits the cache and returns the same data
        let d0b = db.load_dense(0).unwrap();
        assert_eq!(d0b.genome_kmers, d0.genome_kmers);
    }

    fn sample_from(counts: &[(u64, u32)]) -> SequencesSketch {
        let mut s = SequencesSketch::new(String::new(), 50, 31, false, None, 0.0);
        for &(h, c) in counts {
            s.kmer_counts.insert(h, c);
        }
        s
    }

    /// `gather_hits` must reproduce, per genome, the exact coverage multiset that
    /// the per-genome `get_stats` loop collects: intersection of the genome's
    /// sparse k-mers with the sample, with duplicate k-mers counted per
    /// occurrence and shared k-mers credited to every owner.
    #[test]
    fn screen_index_matches_per_genome_intersection() {
        let screen_c = 100usize;
        let thresh = screen_threshold(screen_c);
        // All below threshold so every k-mer is "sparse". g0 has a duplicate (5);
        // 5 and 9 are shared across genomes.
        let sparse = vec![
            vec![5u64, 5, 9, 11],  // g0: duplicate 5
            vec![9u64, 11, 20],    // g1
            vec![5u64, 30, 40, 9], // g2
        ];
        for v in &sparse {
            assert!(v.iter().all(|&h| h < thresh));
        }
        let idx = ScreenIndex::build(&sparse, screen_c, 31);
        assert_eq!(idx.sparse_count, vec![4u32, 3, 4]);

        // Sample: some matching k-mers (with counts), a zero-count k-mer (ignored),
        // a k-mer above threshold (ignored), and a foreign k-mer (no owner).
        let sample = sample_from(&[
            (5, 7),
            (9, 3),
            (11, 2),
            (40, 9),
            (99, 0),         // zero count -> ignored
            (thresh + 1, 5), // above screen threshold -> ignored
            (123456, 4),     // foreign -> no owner
        ]);

        let hits = idx.gather_hits(&sample);

        // Brute-force per-genome reference (mirrors get_stats winner_map=None).
        let mut expected: FxHashMap<u32, Vec<u32>> = FxHashMap::default();
        for (g, v) in sparse.iter().enumerate() {
            for &h in v {
                if h < thresh {
                    if let Some(&c) = sample.kmer_counts.get(&h) {
                        if c != 0 {
                            expected.entry(g as u32).or_default().push(c);
                        }
                    }
                }
            }
        }
        let norm = |m: &FxHashMap<u32, Vec<u32>>| -> Vec<(u32, Vec<u32>)> {
            let mut out: Vec<(u32, Vec<u32>)> = m
                .iter()
                .map(|(&g, v)| {
                    let mut v = v.clone();
                    v.sort_unstable();
                    (g, v)
                })
                .collect();
            out.sort();
            out
        };
        assert_eq!(norm(&hits), norm(&expected));
        // g0 sees 5 twice (duplicate) + 9 + 11 -> covs {7,7,3,2}
        let mut g0 = hits[&0].clone();
        g0.sort_unstable();
        assert_eq!(g0, vec![2, 3, 7, 7]);
    }

    /// Growing a database by copying existing dense blocks verbatim and appending a
    /// new one must produce exactly what writing all the sketches at once produces:
    /// same genome metadata, same decoded dense k-mers, and a stage-1 index that
    /// screens a sample identically (including for the appended genome).
    #[test]
    fn db_add_matches_building_all_at_once() {
        let screen_c = 200usize;
        let thresh = screen_threshold(screen_c);
        let g0 = vec![1u64, 2, thresh - 1, thresh + 10, 9_000_000_000];
        let g1 = vec![7u64, thresh + 1, 123_456_789_000];
        let g2 = vec![3u64, 5, thresh - 2, thresh * 4];
        let s0 = gsketch("g0.fa", g0.clone(), Some(vec![100, 200]));
        let s1 = gsketch("g1.fa", g1.clone(), Some(vec![]));
        let s2 = gsketch("g2.fa", g2.clone(), Some(vec![7, 8, 9]));

        // Two genomes, then grow by copying their blocks through and appending g2.
        let base = open(std::io::Cursor::new(write_db_to_vec(
            &[s0.clone(), s1.clone()],
            screen_c,
        )))
        .unwrap();
        let mut cur = std::io::Cursor::new(Vec::new());
        let mut builder = DbBuilder::new(&mut cur, base.c, base.k, base.screen_c, 1).unwrap();
        for g in 0..base.len() as u32 {
            let (block, sparse) = base.raw_block_and_sparse(g).unwrap();
            builder
                .push_encoded(&block, sparse, base.genome_meta(g))
                .unwrap();
        }
        builder.push_sketch(&s2).unwrap();
        builder.finish().unwrap();
        let grown = open(std::io::Cursor::new(cur.into_inner())).unwrap();

        let all = open(std::io::Cursor::new(write_db_to_vec(
            &[s0, s1, s2],
            screen_c,
        )))
        .unwrap();

        assert_eq!(grown.len(), 3);
        assert_eq!(grown.c, all.c);
        assert_eq!(grown.k, all.k);
        assert_eq!(grown.screen_c, all.screen_c);
        for g in 0..3u32 {
            assert_eq!(grown.genome_meta(g), all.genome_meta(g));
            let a = grown.decode_dense(g).unwrap();
            let b = all.decode_dense(g).unwrap();
            assert_eq!(a, b, "dense block {} differs after db-add", g);
        }
        assert_eq!(
            grown.screen_index.sparse_count,
            all.screen_index.sparse_count
        );

        // The rebuilt stage-1 index must screen identically, including for the
        // appended genome's k-mers.
        let sample = sample_from(&[
            (1, 4),
            (2, 2),
            (3, 6),
            (5, 1),
            (7, 3),
            (thresh - 1, 5),
            (thresh - 2, 8),
        ]);
        let norm = |m: FxHashMap<u32, Vec<u32>>| -> Vec<(u32, Vec<u32>)> {
            let mut out: Vec<(u32, Vec<u32>)> = m
                .into_iter()
                .map(|(g, mut v)| {
                    v.sort_unstable();
                    (g, v)
                })
                .collect();
            out.sort();
            out
        };
        let from_grown = norm(grown.screen_index.gather_hits(&sample));
        assert_eq!(from_grown, norm(all.screen_index.gather_hits(&sample)));
        // g2 was the appended genome; it must be screened, not silently absent.
        assert!(
            from_grown.iter().any(|(g, _)| *g == 2),
            "appended genome missing from the rebuilt screen index"
        );
    }

    /// The index survives a serialize/parse round-trip through the file format.
    #[test]
    fn screen_index_roundtrips_through_db() {
        let thresh = u64::MAX / 200;
        let g0 = vec![5u64, 5, thresh - 1, thresh + 9]; // last is above screen thresh
        let g1 = vec![5u64, 7, thresh - 2];
        let sketches = vec![
            gsketch("g0.fa", g0.clone(), Some(vec![1])),
            gsketch("g1.fa", g1.clone(), Some(vec![2])),
        ];
        let db = open(std::io::Cursor::new(write_db_to_vec(&sketches, 200))).unwrap();

        let sample = sample_from(&[(5, 4), (7, 6), (thresh - 1, 1), (thresh - 2, 9)]);
        let hits = db.screen_index.gather_hits(&sample);
        // g0: 5 twice + (thresh-1) once -> {4,4,1}; g1: 5 + 7 + (thresh-2) -> {4,6,9}
        let mut g0h = hits[&0].clone();
        g0h.sort_unstable();
        assert_eq!(g0h, vec![1, 4, 4]);
        let mut g1h = hits[&1].clone();
        g1h.sort_unstable();
        assert_eq!(g1h, vec![4, 6, 9]);
    }

    /// Both header layouts are parsed, with the offsets read from the right place in
    /// each: upstream sylph's version 2 has no checksum field, so its offsets sit 8
    /// bytes earlier than in the version 3 written here.
    #[test]
    fn parses_both_header_versions() {
        let mut v3 = Vec::new();
        v3.extend_from_slice(MAGIC);
        v3.push(VERSION);
        v3.extend_from_slice(&0xdead_beef_u64.to_le_bytes()); // checksum
        v3.extend_from_slice(&111u64.to_le_bytes()); // index offset
        v3.extend_from_slice(&222u64.to_le_bytes()); // footer offset
        assert_eq!(v3.len(), HEADER_LEN as usize);
        let h = parse_header(&v3).unwrap();
        assert_eq!(h.len, HEADER_LEN);
        assert_eq!(h.checksum, Some(0xdead_beef));
        assert_eq!((h.index_offset, h.footer_offset), (111, 222));

        let mut v2 = Vec::new();
        v2.extend_from_slice(MAGIC);
        v2.push(VERSION_NO_CHECKSUM);
        v2.extend_from_slice(&111u64.to_le_bytes());
        v2.extend_from_slice(&222u64.to_le_bytes());
        assert_eq!(v2.len(), HEADER_LEN_V2 as usize);
        let h = parse_header(&v2).unwrap();
        assert_eq!(h.len, HEADER_LEN_V2);
        assert_eq!(h.checksum, None);
        assert_eq!((h.index_offset, h.footer_offset), (111, 222));
        // A v2 file is exactly HEADER_LEN_V2 long here, i.e. shorter than the buffer a
        // v3 read would fill; that must not be mistaken for truncation.
        assert!(parse_header(&v2[..v2.len() - 1]).is_err());

        // Unknown version, and a v3 header cut short.
        let mut future = v3.clone();
        future[4] = VERSION + 1;
        assert!(parse_header(&future).is_err());
        assert!(parse_header(&v3[..HEADER_LEN as usize - 1]).is_err());
        assert!(parse_header(b"NOPE").is_err());
    }

    /// A genome too small to reach `--min-sparse-kmers` at the nominal `--screen-c`
    /// gets a denser, genome-specific screen rate, and the database-wide rate recorded
    /// in the file is loosened just enough that its k-mers are still matchable.
    #[test]
    fn small_genome_gets_denser_screen_rate() {
        let screen_c = 3000usize;
        let thresh = screen_threshold(screen_c);
        // A big genome with plenty of sparse k-mers, and a small one with none at all
        // at screen_c (every hash above the nominal threshold).
        let big: Vec<u64> = (0..60u64).map(|i| thresh / 2 + i).collect();
        let small: Vec<u64> = (0..10u64).map(|i| thresh * 7 + i * 13).collect();
        let sketches = vec![
            gsketch("big.fa", big.clone(), Some(vec![1])),
            gsketch("small.fa", small.clone(), Some(vec![2])),
        ];

        // Plain FracMinHash selection leaves the small genome with no screen k-mers at
        // all, i.e. invisible to the stage-1 screen however deeply it is covered.
        assert!(sparse_subset(&small, screen_c).is_empty());

        // With a floor of 5 it stores its 5 smallest k-mers, and the recorded screen
        // rate is denser so those k-mers pass `gather_hits`' threshold.
        let on = open(std::io::Cursor::new(write_db_to_vec_min_sparse(
            &sketches, screen_c, 5,
        )))
        .unwrap();
        assert_eq!(on.screen_index.sparse_count[1], 5);
        assert!(
            on.screen_c < screen_c,
            "effective screen -c should have been densified, got {}",
            on.screen_c
        );
        // The big genome is unaffected: it clears the floor at the nominal rate, so it
        // keeps exactly its FracMinHash subset.
        assert_eq!(
            on.screen_index.sparse_count[0] as usize,
            sparse_subset(&big, screen_c).len()
        );
        let smallest_five = smallest_n(&small, 5);
        let hits = on
            .screen_index
            .gather_hits(&sample_from(&[(smallest_five[0], 4)]));
        assert_eq!(hits.get(&1).map(|v| v.as_slice()), Some(&[4u32][..]));
    }

    /// `db-add` must reproduce a densified genome's screen entry exactly. Its sparse
    /// set is recovered from the stored count, not by re-subsampling at the database's
    /// screen rate -- which, for a genome selected at a denser rate, would hand back a
    /// different set and change every genome's screen ANI denominator.
    #[test]
    fn db_add_preserves_densified_sparse_sets() {
        let screen_c = 3000usize;
        let thresh = screen_threshold(screen_c);
        let big: Vec<u64> = (0..60u64).map(|i| thresh / 2 + i).collect();
        let small: Vec<u64> = (0..10u64).map(|i| thresh * 7 + i * 13).collect();
        let s0 = gsketch("big.fa", big, Some(vec![1]));
        let s1 = gsketch("small.fa", small, Some(vec![2]));
        let s2 = gsketch(
            "extra.fa",
            (0..80u64).map(|i| thresh / 3 + i).collect(),
            None,
        );

        let base = open(std::io::Cursor::new(write_db_to_vec_min_sparse(
            &[s0.clone(), s1.clone()],
            screen_c,
            5,
        )))
        .unwrap();
        let mut cur = std::io::Cursor::new(Vec::new());
        let mut builder = DbBuilder::new(&mut cur, base.c, base.k, base.screen_c, 5).unwrap();
        for g in 0..base.len() as u32 {
            let (block, sparse) = base.raw_block_and_sparse(g).unwrap();
            builder
                .push_encoded(&block, sparse, base.genome_meta(g))
                .unwrap();
        }
        builder.push_sketch(&s2).unwrap();
        builder.finish().unwrap();
        let grown = open(std::io::Cursor::new(cur.into_inner())).unwrap();

        // The copied genomes keep exactly the screen entries they were built with.
        assert_eq!(
            &grown.screen_index.sparse_count[..2],
            &base.screen_index.sparse_count[..]
        );
        assert_eq!(grown.screen_c, base.screen_c);
        let sample = sample_from(&[(smallest_n(&s1.genome_kmers, 5)[0], 3)]);
        assert_eq!(
            grown.screen_index.gather_hits(&sample).get(&1),
            base.screen_index.gather_hits(&sample).get(&1)
        );
    }
}
