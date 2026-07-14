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
//! ## Integrity
//!
//! Profiling only ever touches a few blocks of the file, so a corrupt database is
//! not detected as a side effect of reading it (unlike the single-frame `.sylspc` /
//! `.sylspr` formats, whose zstd checksum is validated on every decode). The header
//! instead carries an XXH64 of the rest of the file, checked on demand by
//! [`TwoStageDb::verify_checksum`] — which `weebill inspect` does.

use crate::cmdline::DbConvertArgs;
use crate::constants::*;
use crate::types::*;
use boomphf::Mphf;
use fxhash::FxHashMap;
use log::*;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::os::unix::fs::FileExt;
use std::path::Path;
use std::sync::{Arc, Mutex};

const MAGIC: &[u8; 4] = b"SY2D";
/// Format version (dense Golomb-Rice blocks + pooled-MPHF stage-1 screen index).
/// Only this exact version is readable — an older database must be rebuilt with
/// `db-convert`.
const VERSION: u8 = 3;
/// magic (4) + version (1) + XXH64 of the rest of the file (8) + index offset (8)
/// + footer offset (8)
const HEADER_LEN: u64 = 29;
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

/// Re-pack genome sketches into the two-stage seekable layout and write to `w`.
/// `screen_c` is the (coarser) stage-1 subsampling rate; it must be `>= c`.
/// Dense blocks are Golomb-Rice coded.
pub fn write_two_stage_db<W: Write>(
    mut w: W,
    sketches: &[GenomeSketch],
    screen_c: usize,
) -> io::Result<()> {
    let c = sketches.first().map(|s| s.c).unwrap_or(0);
    let k = sketches.first().map(|s| s.k).unwrap_or(0);
    // FracMinHash threshold for the sparse subsample (same rule as subsample_view).
    let thresh = if screen_c == 0 {
        u64::MAX
    } else {
        u64::MAX / screen_c as u64
    };

    let mut body: Vec<u8> = Vec::new();
    let mut genomes: Vec<GenomeMeta> = Vec::with_capacity(sketches.len());
    let mut sparse_per_genome: Vec<Vec<u64>> = Vec::with_capacity(sketches.len());

    for gs in sketches {
        let dense_offset = HEADER_LEN + body.len() as u64;
        write_hashes(&mut body, &gs.genome_kmers);
        match &gs.pseudotax_tracked_nonused_kmers {
            Some(p) => {
                body.push(1);
                write_hashes(&mut body, p);
            }
            None => body.push(0),
        }
        sparse_per_genome.push(
            gs.genome_kmers
                .iter()
                .copied()
                .filter(|&h| h < thresh)
                .collect(),
        );
        genomes.push(GenomeMeta {
            file_name: gs.file_name.clone(),
            first_contig_name: gs.first_contig_name.clone(),
            gn_size: gs.gn_size,
            min_spacing: gs.min_spacing,
            has_pseudotax: gs.pseudotax_tracked_nonused_kmers.is_some(),
            dense_offset,
        });
    }

    // Pooled stage-1 screen index, then footer metadata.
    let screen_index = ScreenIndex::build(&sparse_per_genome, screen_c, k);
    drop(sparse_per_genome);
    let mut index_block: Vec<u8> = Vec::new();
    screen_index.write_to_vec(&mut index_block)?;

    let footer = Footer {
        c,
        k,
        screen_c,
        genomes,
    };
    let footer_bytes = bincode::serialize(&footer).map_err(io::Error::other)?;
    let index_offset = HEADER_LEN + body.len() as u64;
    let footer_offset = index_offset + index_block.len() as u64;

    // Everything after the header is checksummed, so `inspect` can detect a
    // truncated or bit-rotted database that a seeking reader would otherwise decode
    // into a plausible-looking but wrong sketch.
    let checksum = {
        let mut h = crate::checksum::HashingWriter::new(io::sink());
        h.write_all(&body)?;
        h.write_all(&index_block)?;
        h.write_all(&footer_bytes)?;
        h.finish().1
    };

    w.write_all(MAGIC)?;
    w.write_all(&[VERSION])?;
    w.write_all(&checksum.to_le_bytes())?;
    w.write_all(&index_offset.to_le_bytes())?;
    w.write_all(&footer_offset.to_le_bytes())?;
    w.write_all(&body)?;
    w.write_all(&index_block)?;
    w.write_all(&footer_bytes)?;
    Ok(())
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
    /// read of a file that profiling otherwise only touches a few blocks of.
    checksum: u64,
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

/// Parse the magic + version header; return `(checksum, index_offset, footer_offset)`.
fn parse_header(hdr: &[u8]) -> io::Result<(u64, u64, u64)> {
    if hdr.len() < HEADER_LEN as usize || &hdr[0..4] != MAGIC {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "not a sylph two-stage database",
        ));
    }
    if hdr[4] != VERSION {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "two-stage database is version {}, but only version {} is readable; rebuild it with db-convert",
                hdr[4], VERSION
            ),
        ));
    }
    let checksum = u64::from_le_bytes(hdr[5..13].try_into().unwrap());
    let index_offset = u64::from_le_bytes(hdr[13..21].try_into().unwrap());
    let footer_offset = u64::from_le_bytes(hdr[21..29].try_into().unwrap());
    Ok((checksum, index_offset, footer_offset))
}

/// Assemble a `TwoStageDb` from its parsed footer + screen index + backing store.
fn build_db(
    footer: Footer,
    checksum: u64,
    index_offset: u64,
    screen_index: ScreenIndex,
    data: DenseData,
) -> TwoStageDb {
    TwoStageDb {
        c: footer.c,
        k: footer.k,
        screen_c: footer.screen_c,
        checksum,
        index_offset,
        genomes: footer.genomes,
        screen_index,
        data,
        cache: Mutex::new(FxHashMap::default()),
    }
}

/// Parse the header + index + footer of a `.syl2db` already resident in `data`.
fn from_bytes(data: DenseData) -> io::Result<TwoStageDb> {
    let bytes = data.bytes();
    let (checksum, index_offset, footer_offset) = parse_header(bytes)?;
    if index_offset > footer_offset || footer_offset as usize > bytes.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "two-stage database offsets out of range",
        ));
    }
    let footer: Footer = bincode::deserialize(&bytes[footer_offset as usize..])
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    let screen_index = ScreenIndex::read(
        &bytes[index_offset as usize..footer_offset as usize],
        footer.screen_c,
        footer.k,
    )?;
    Ok(build_db(footer, checksum, index_offset, screen_index, data))
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

    /// Re-hash the whole file and compare against the checksum in its header. This
    /// reads every byte, so it is on-demand (`weebill inspect`) rather than part of
    /// opening the database.
    pub fn verify_checksum(&self) -> io::Result<()> {
        let got = match &self.data {
            DenseData::Owned(v) => crate::checksum::hash_reader(&v[HEADER_LEN as usize..])?,
            DenseData::File(file) => {
                let mut r = BufReader::with_capacity(1 << 20, file.try_clone()?);
                r.seek(SeekFrom::Start(HEADER_LEN))?;
                crate::checksum::hash_reader(r)?
            }
        };
        if got != self.checksum {
            return Err(crate::checksum::mismatch(
                "the two-stage database",
                self.checksum,
                got,
            ));
        }
        Ok(())
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
    file.read_exact_at(&mut hdr, 0)?;
    let (checksum, index_offset, footer_offset) = parse_header(&hdr)?;
    let flen = file.metadata()?.len();
    if index_offset > footer_offset || footer_offset > flen {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "two-stage database offsets out of range",
        ));
    }
    let mut fbytes = vec![0u8; (flen - footer_offset) as usize];
    file.read_exact_at(&mut fbytes, footer_offset)?;
    let footer: Footer =
        bincode::deserialize(&fbytes).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    let mut ibytes = vec![0u8; (footer_offset - index_offset) as usize];
    file.read_exact_at(&mut ibytes, index_offset)?;
    let screen_index = ScreenIndex::read(&ibytes, footer.screen_c, footer.k)?;
    Ok(build_db(
        footer,
        checksum,
        index_offset,
        screen_index,
        DenseData::File(file),
    ))
}

// --- CLI handler ------------------------------------------------------------

fn load_genome_sketches(path: &str) -> Vec<GenomeSketch> {
    let file = File::open(path).unwrap_or_else(|_| panic!("Could not open {}", path));
    let reader = BufReader::with_capacity(10_000_000, file);
    bincode::deserialize_from(reader)
        .unwrap_or_else(|_| panic!("{} is not a valid database sketch (.syldb)", path))
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
    write_two_stage_db(w, &sketches, args.screen_c)
        .unwrap_or_else(|e| panic!("Failed to write {}: {}", out, e));
    info!("Wrote two-stage database to {}", out);
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

        let mut buf = Vec::new();
        write_two_stage_db(&mut buf, &sketches, 200).unwrap();
        let db = open(std::io::Cursor::new(buf)).unwrap();

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
        let mut buf = Vec::new();
        write_two_stage_db(&mut buf, &sketches, 200).unwrap();
        let db = open(std::io::Cursor::new(buf)).unwrap();

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
}
