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
//!     a k-mer shared by several genomes is stored in each of them. We do not
//!     assign each k-mer a single owner.
//!   * **No shared/pooled hashes.** There is no conserved-k-mer pool.
//!
//! ## Layout
//!
//! ```text
//! [4]  magic  "SY2D"
//! [1]  version
//! [8]  footer_offset (u64 LE)
//! ---- body ----
//! dense block 0, dense block 1, ...      (each: Golomb-Rice compressed)
//! ---- footer ----
//! bincode(Footer)                        (the sparse stage-1 index + metadata)
//! ```
//!
//! **Stage 1 (sparse, bincoded, loaded fully).** The footer is a bincoded
//! `Footer` holding, per genome, a FracMinHash subsample of its k-mers at a
//! coarser rate `screen_c` (the same thresholding `subsample_view` uses, so the
//! sparse set is a genuine sparser sketch). Querying a sample against these tiny
//! sketches yields the genomes it contains.
//!
//! **Stage 2 (dense, Golomb-Rice, loaded on demand).** Each genome's *full*
//! `genome_kmers` and `pseudotax_tracked_nonused_kmers` are an independently
//! Golomb-Rice-coded block at a known offset. Only the genomes that pass the
//! stage-1 screen are decoded (and cached across samples) to reconstruct their
//! exact `GenomeSketch` for the dense profiling pass.

use crate::cmdline::DbConvertArgs;
use crate::constants::*;
use crate::types::*;
use fxhash::FxHashMap;
use log::*;
use memmap2::Mmap;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::os::unix::fs::FileExt;
use std::path::Path;
use std::sync::{Arc, Mutex};

const MAGIC: &[u8; 4] = b"SY2D";
/// Dense blocks Golomb-Rice coded.
const VERSION_RICE: u8 = 1;
/// Dense blocks Elias-Fano coded.
const VERSION_EF: u8 = 2;
/// magic (4) + version (1) + footer offset (8)
const HEADER_LEN: u64 = 13;

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

struct BitReader<'a> {
    buf: &'a [u8],
    byte: usize,
    bit: u8,
}

impl<'a> BitReader<'a> {
    fn new(buf: &'a [u8]) -> Self {
        BitReader {
            buf,
            byte: 0,
            bit: 0,
        }
    }
    #[inline]
    fn read_bit(&mut self) -> io::Result<u32> {
        if self.byte >= self.buf.len() {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "bitstream truncated",
            ));
        }
        let b = (self.buf[self.byte] >> self.bit) & 1;
        self.bit += 1;
        if self.bit == 8 {
            self.bit = 0;
            self.byte += 1;
        }
        Ok(b as u32)
    }
    #[inline]
    fn read_bits(&mut self, n: u32) -> io::Result<u64> {
        let mut v = 0u64;
        for i in 0..n {
            v |= (self.read_bit()? as u64) << i;
        }
        Ok(v)
    }
    #[inline]
    fn read_unary(&mut self) -> io::Result<u64> {
        let mut q = 0u64;
        while self.read_bit()? == 1 {
            q += 1;
        }
        Ok(q)
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

/// Sort + Elias-Fano encode a set of hashes onto `out`. The sorted sequence is
/// split into high bits (a unary bucket bitvector, ~2 bits/element) and `l` low
/// bits stored verbatim, where `l = floor(log2(universe / n))`. Self-delimiting
/// given the leading count. Duplicates are preserved (zero high-gap, repeated 1).
fn write_hashes_ef(out: &mut Vec<u8>, hashes: &[u64]) {
    let mut sorted = hashes.to_vec();
    sorted.sort_unstable();
    let n = sorted.len();
    write_uvarint(out, n as u64);
    if n == 0 {
        return;
    }
    let universe = sorted[n - 1].saturating_add(1);
    // l = floor(log2(universe / n)), clamped to [0, 63].
    let l: u32 = {
        let q = universe / n as u64;
        if q < 2 {
            0
        } else {
            63 - q.leading_zeros()
        }
    };
    out.push(l as u8);
    let mut bw = BitWriter::new();
    // High bitvector: for element i emit (high_i - high_{i-1}) zeros then a 1.
    let mut prev_high = 0u64;
    for &v in &sorted {
        let high = v >> l;
        for _ in 0..(high - prev_high) {
            bw.write_bit(0);
        }
        bw.write_bit(1);
        prev_high = high;
    }
    // Low bits: l bits per element, in order.
    if l > 0 {
        let mask = (1u64 << l) - 1;
        for &v in &sorted {
            bw.write_bits(v & mask, l);
        }
    }
    let bits = bw.finish();
    write_uvarint(out, bits.len() as u64);
    out.extend_from_slice(&bits);
}

fn read_hashes_ef<R: Read>(r: &mut R) -> io::Result<Vec<u64>> {
    let n = read_uvarint(r)? as usize;
    if n == 0 {
        return Ok(Vec::new());
    }
    let mut lb = [0u8; 1];
    r.read_exact(&mut lb)?;
    let l = lb[0] as u32;
    let blen = read_uvarint(r)? as usize;
    let mut bits = vec![0u8; blen];
    r.read_exact(&mut bits)?;
    let mut br = BitReader::new(&bits);
    // Recover high parts by scanning the bucket bitvector: 0 advances the bucket,
    // 1 emits an element at the current bucket value.
    let mut highs = Vec::with_capacity(n);
    let mut cur = 0u64;
    while highs.len() < n {
        if br.read_bit()? == 1 {
            highs.push(cur);
        } else {
            cur += 1;
        }
    }
    let mut out = Vec::with_capacity(n);
    if l > 0 {
        for &h in &highs {
            let low = br.read_bits(l)?;
            out.push((h << l) | low);
        }
    } else {
        out.extend_from_slice(&highs);
    }
    Ok(out)
}

/// Dispatch dense-block encode/decode on the database codec.
fn write_dense_block(out: &mut Vec<u8>, hashes: &[u64], ef: bool) {
    if ef {
        write_hashes_ef(out, hashes)
    } else {
        write_hashes(out, hashes)
    }
}

fn read_dense_block<R: Read>(r: &mut R, ef: bool) -> io::Result<Vec<u64>> {
    if ef {
        read_hashes_ef(r)
    } else {
        read_hashes(r)
    }
}

// --- footer (stage-1 sparse index + metadata) -------------------------------

/// Per-genome metadata plus the stage-1 sparse k-mer subsample. Everything here
/// is loaded into memory when the database is opened; the dense block at
/// `dense_offset` is decoded lazily.
#[derive(serde::Serialize, serde::Deserialize, Clone, Debug, PartialEq)]
pub struct GenomeMeta {
    pub file_name: String,
    pub first_contig_name: String,
    pub gn_size: usize,
    pub min_spacing: usize,
    pub has_pseudotax: bool,
    pub dense_offset: u64,
    /// FracMinHash subsample of `genome_kmers` at `screen_c` (stage-1 screen).
    pub sparse_kmers: Vec<u64>,
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

// --- writing ----------------------------------------------------------------

/// Re-pack genome sketches into the two-stage seekable layout (Golomb-Rice dense
/// blocks). See [`write_two_stage_db_codec`] to choose the dense codec.
pub fn write_two_stage_db<W: Write>(
    w: W,
    sketches: &[GenomeSketch],
    screen_c: usize,
) -> io::Result<()> {
    write_two_stage_db_codec(w, sketches, screen_c, false)
}

/// Re-pack genome sketches into the two-stage seekable layout and write to `w`.
/// `screen_c` is the (coarser) stage-1 subsampling rate; it must be `>= c`.
/// When `use_ef`, dense blocks are Elias-Fano coded (version 2) instead of
/// Golomb-Rice (version 1).
pub fn write_two_stage_db_codec<W: Write>(
    mut w: W,
    sketches: &[GenomeSketch],
    screen_c: usize,
    use_ef: bool,
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

    for gs in sketches {
        let dense_offset = HEADER_LEN + body.len() as u64;
        write_dense_block(&mut body, &gs.genome_kmers, use_ef);
        match &gs.pseudotax_tracked_nonused_kmers {
            Some(p) => {
                body.push(1);
                write_dense_block(&mut body, p, use_ef);
            }
            None => body.push(0),
        }
        let sparse_kmers: Vec<u64> = gs
            .genome_kmers
            .iter()
            .copied()
            .filter(|&h| h < thresh)
            .collect();
        genomes.push(GenomeMeta {
            file_name: gs.file_name.clone(),
            first_contig_name: gs.first_contig_name.clone(),
            gn_size: gs.gn_size,
            min_spacing: gs.min_spacing,
            has_pseudotax: gs.pseudotax_tracked_nonused_kmers.is_some(),
            dense_offset,
            sparse_kmers,
        });
    }

    let footer = Footer {
        c,
        k,
        screen_c,
        genomes,
    };
    let footer_bytes = bincode::serialize(&footer).map_err(io::Error::other)?;
    let footer_offset = HEADER_LEN + body.len() as u64;

    let version = if use_ef { VERSION_EF } else { VERSION_RICE };
    w.write_all(MAGIC)?;
    w.write_all(&[version])?;
    w.write_all(&footer_offset.to_le_bytes())?;
    w.write_all(&body)?;
    w.write_all(&footer_bytes)?;
    Ok(())
}

// --- opened database --------------------------------------------------------

/// Backing store for the dense blocks.
///   * `Mmap`  - whole file memory-mapped; a block is a zero-copy slice. Used for
///     multi-threaded runs (parallel page-fault decode), at the cost of resident
///     file pages counting toward RSS.
///   * `File`  - positional `read_at` of just the requested block bytes. Used for
///     single-threaded runs: no parallelism to gain from mmap, and it keeps RSS
///     low (only the block, plus reclaimable OS page cache).
///   * `Owned` - whole file in memory; for in-memory readers / tests.
enum DenseData {
    Mmap(Mmap),
    File(File),
    Owned(Vec<u8>),
}

impl DenseData {
    /// Whole-file bytes; only valid for the in-memory backings (used to parse the
    /// header/footer). The `File` backing is read positionally via `with_block`.
    #[inline]
    fn bytes(&self) -> &[u8] {
        match self {
            DenseData::Mmap(m) => &m[..],
            DenseData::Owned(v) => &v[..],
            DenseData::File(_) => unreachable!("File-backed db is read via with_block"),
        }
    }
}

/// An opened two-stage database. Construction loads only the stage-1 sparse
/// index (the bincoded footer); per-genome dense blocks are decoded on demand,
/// in parallel (each decode reads its own slice of the memory map, no shared
/// cursor or lock), and cached.
pub struct TwoStageDb {
    pub c: usize,
    pub k: usize,
    pub screen_c: usize,
    /// Dense blocks are Elias-Fano (version 2) rather than Golomb-Rice (v1).
    ef: bool,
    /// File offset where the dense-block region ends (start of the footer); used
    /// to bound the last genome's block for positional reads.
    footer_offset: u64,
    genomes: Vec<GenomeMeta>,
    /// Stage-1 sparse sketches, one per genome, in database order. These stand in
    /// for the full database during the cheap screen.
    pub screen_sketches: Vec<GenomeSketch>,
    data: DenseData,
    cache: Mutex<FxHashMap<u32, Arc<GenomeSketch>>>,
}

/// Parse the magic + version header; return the dense codec flag and the footer
/// offset.
fn parse_header(hdr: &[u8]) -> io::Result<(bool, u64)> {
    if hdr.len() < HEADER_LEN as usize || &hdr[0..4] != MAGIC {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "not a sylph two-stage database",
        ));
    }
    let ef = match hdr[4] {
        VERSION_RICE => false,
        VERSION_EF => true,
        _ => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "unsupported two-stage database version",
            ))
        }
    };
    let footer_offset = u64::from_le_bytes(hdr[5..13].try_into().unwrap());
    Ok((ef, footer_offset))
}

/// Assemble a `TwoStageDb` from its parsed footer + backing store.
fn build_db(footer: Footer, ef: bool, footer_offset: u64, data: DenseData) -> TwoStageDb {
    // Build the stage-1 screen sketches from the sparse subsamples. The screen
    // pass runs with pseudotax disabled, so an empty tracked-kmer set is a safe
    // placeholder that keeps the `--disable-profiling` guard in `contain` happy;
    // the true pseudotax k-mers come from the dense blocks.
    let screen_sketches: Vec<GenomeSketch> = footer
        .genomes
        .iter()
        .map(|m| GenomeSketch {
            genome_kmers: m.sparse_kmers.clone(),
            pseudotax_tracked_nonused_kmers: Some(Vec::new()),
            file_name: m.file_name.clone(),
            first_contig_name: m.first_contig_name.clone(),
            c: footer.screen_c,
            k: footer.k,
            gn_size: m.gn_size,
            min_spacing: m.min_spacing,
        })
        .collect();
    TwoStageDb {
        c: footer.c,
        k: footer.k,
        screen_c: footer.screen_c,
        ef,
        footer_offset,
        genomes: footer.genomes,
        screen_sketches,
        data,
        cache: Mutex::new(FxHashMap::default()),
    }
}

/// Parse the header + footer of a `.syl2db` already resident in `data`.
fn from_bytes(data: DenseData) -> io::Result<TwoStageDb> {
    let bytes = data.bytes();
    let (ef, footer_offset) = parse_header(bytes)?;
    if footer_offset as usize > bytes.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "two-stage database footer offset out of range",
        ));
    }
    let footer: Footer = bincode::deserialize(&bytes[footer_offset as usize..])
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    Ok(build_db(footer, ef, footer_offset, data))
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

    /// End offset of genome `g`'s region (start of the next genome's block, or
    /// the footer for the last genome). Genomes are stored in ascending offset.
    fn block_end(&self, g: u32) -> u64 {
        let gi = g as usize;
        if gi + 1 < self.genomes.len() {
            self.genomes[gi + 1].dense_offset
        } else {
            self.footer_offset
        }
    }

    /// Run `f` on genome `g`'s whole on-disk region (`genome_kmers` block, the
    /// pseudotax flag, and the optional pseudotax block). For the mmap/owned
    /// backings this is a zero-copy slice; for the file backing it is one
    /// positional read of just that region. No shared cursor, so it is safe to
    /// call concurrently from many threads.
    fn with_block<T>(&self, g: u32, f: impl FnOnce(&[u8]) -> io::Result<T>) -> io::Result<T> {
        let start = self.genomes[g as usize].dense_offset as usize;
        match &self.data {
            DenseData::Mmap(m) => f(&m[start..]),
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
    /// cursor (mmap slice or positional read). The two-stage pass-1 uses this to
    /// decode each survivor into a short-lived sketch, probe it, and drop it
    /// unless it passes -- so the discarded majority is never cached.
    pub fn decode_dense(&self, g: u32) -> io::Result<GenomeSketch> {
        let meta = &self.genomes[g as usize];
        let (genome_kmers, pseudotax) = self.with_block(g, |bytes| {
            let mut cur = bytes;
            let gk = read_dense_block(&mut cur, self.ef)?;
            let mut flag = [0u8; 1];
            cur.read_exact(&mut flag)?;
            let pt = if flag[0] != 0 {
                Some(read_dense_block(&mut cur, self.ef)?)
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

/// Open a `.syl2db` file from a path. With `mmap`, the whole file is memory
/// mapped for lock-free parallel decode (higher RSS); otherwise blocks are read
/// positionally on demand (lower RSS, preferred when single-threaded).
pub fn open_file(path: &str, mmap: bool) -> io::Result<TwoStageDb> {
    let file = File::open(path)?;
    if mmap {
        // Safety: the database is treated as read-only; we never mutate the map
        // and the file is not expected to change underneath us during a run.
        let m = unsafe { Mmap::map(&file)? };
        return from_bytes(DenseData::Mmap(m));
    }
    let mut hdr = [0u8; HEADER_LEN as usize];
    file.read_exact_at(&mut hdr, 0)?;
    let (ef, footer_offset) = parse_header(&hdr)?;
    let flen = file.metadata()?.len();
    if footer_offset > flen {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "two-stage database footer offset out of range",
        ));
    }
    let mut fbytes = vec![0u8; (flen - footer_offset) as usize];
    file.read_exact_at(&mut fbytes, footer_offset)?;
    let footer: Footer = bincode::deserialize(&fbytes)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    Ok(build_db(footer, ef, footer_offset, DenseData::File(file)))
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
        "Converting {} genomes (dense -c {}, stage-1 screen -c {}, dense codec {}) -> {}",
        sketches.len(),
        c,
        args.screen_c,
        if args.ef { "elias-fano" } else { "golomb-rice" },
        out
    );
    let w =
        BufWriter::new(File::create(&out).unwrap_or_else(|_| panic!("Could not create {}", out)));
    write_two_stage_db_codec(w, &sketches, args.screen_c, args.ef)
        .unwrap_or_else(|e| panic!("Failed to write {}: {}", out, e));
    info!("Wrote two-stage database to {}", out);
}

#[cfg(test)]
mod tests {
    use super::*;

    fn roundtrip_hashes(input: &[u64]) {
        let mut expected = input.to_vec();
        expected.sort_unstable();
        // Golomb-Rice
        let mut buf = Vec::new();
        write_hashes(&mut buf, input);
        let mut r = &buf[..];
        assert_eq!(read_hashes(&mut r).unwrap(), expected);
        assert!(r.is_empty(), "read_hashes left trailing bytes");
        // Elias-Fano
        let mut buf = Vec::new();
        write_hashes_ef(&mut buf, input);
        let mut r = &buf[..];
        assert_eq!(read_hashes_ef(&mut r).unwrap(), expected);
        assert!(r.is_empty(), "read_hashes_ef left trailing bytes");
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

        // stage-1 sparse sketch = fracminhash subset at screen_c
        let expect_sparse = |ks: &[u64]| -> Vec<u64> {
            let mut v: Vec<u64> = ks.iter().copied().filter(|&h| h < thresh).collect();
            v.sort_unstable();
            v
        };
        let mut s0 = db.screen_sketches[0].genome_kmers.clone();
        s0.sort_unstable();
        assert_eq!(s0, expect_sparse(&g0_kmers));
        assert_eq!(db.screen_sketches[0].c, 200);

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
}
