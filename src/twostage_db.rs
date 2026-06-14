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
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::Path;
use std::sync::{Arc, Mutex};

const MAGIC: &[u8; 4] = b"SY2D";
const VERSION: u8 = 1;
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

/// Re-pack genome sketches into the two-stage seekable layout and write to `w`.
/// `screen_c` is the (coarser) stage-1 subsampling rate; it must be `>= c`.
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

    w.write_all(MAGIC)?;
    w.write_all(&[VERSION])?;
    w.write_all(&footer_offset.to_le_bytes())?;
    w.write_all(&body)?;
    w.write_all(&footer_bytes)?;
    Ok(())
}

// --- opened, seekable database ----------------------------------------------

pub trait ReadSeek: Read + Seek + Send {}
impl<T: Read + Seek + Send> ReadSeek for T {}

/// An opened two-stage database. Construction loads only the stage-1 sparse
/// index (the bincoded footer); per-genome dense blocks are decoded on demand
/// and cached.
pub struct TwoStageDb {
    pub c: usize,
    pub k: usize,
    pub screen_c: usize,
    genomes: Vec<GenomeMeta>,
    /// Stage-1 sparse sketches, one per genome, in database order. These stand in
    /// for the full database during the cheap screen.
    pub screen_sketches: Vec<GenomeSketch>,
    reader: Mutex<Box<dyn ReadSeek>>,
    cache: Mutex<FxHashMap<u32, Arc<GenomeSketch>>>,
}

/// Open a `.syl2db`, loading the stage-1 sparse index.
pub fn open<R: Read + Seek + Send + 'static>(mut r: R) -> io::Result<TwoStageDb> {
    let mut hdr = [0u8; HEADER_LEN as usize];
    r.read_exact(&mut hdr)?;
    if &hdr[0..4] != MAGIC {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "not a sylph two-stage database",
        ));
    }
    if hdr[4] != VERSION {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "unsupported two-stage database version",
        ));
    }
    let footer_offset = u64::from_le_bytes(hdr[5..13].try_into().unwrap());

    r.seek(SeekFrom::Start(footer_offset))?;
    let mut fbytes = Vec::new();
    r.read_to_end(&mut fbytes)?;
    let footer: Footer =
        bincode::deserialize(&fbytes).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

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

    Ok(TwoStageDb {
        c: footer.c,
        k: footer.k,
        screen_c: footer.screen_c,
        genomes: footer.genomes,
        screen_sketches,
        reader: Mutex::new(Box::new(r)),
        cache: Mutex::new(FxHashMap::default()),
    })
}

impl TwoStageDb {
    pub fn len(&self) -> usize {
        self.genomes.len()
    }
    pub fn is_empty(&self) -> bool {
        self.genomes.is_empty()
    }

    /// Decode genome `g`'s full dense `GenomeSketch` (cached across calls).
    pub fn load_dense(&self, g: u32) -> io::Result<Arc<GenomeSketch>> {
        if let Some(a) = self.cache.lock().unwrap().get(&g) {
            return Ok(a.clone());
        }
        let meta = &self.genomes[g as usize];
        let (genome_kmers, pseudotax) = {
            let mut rd = self.reader.lock().unwrap();
            rd.seek(SeekFrom::Start(meta.dense_offset))?;
            let mut r: &mut dyn ReadSeek = &mut **rd;
            let gk = read_hashes(&mut r)?;
            let mut flag = [0u8; 1];
            r.read_exact(&mut flag)?;
            let pt = if flag[0] != 0 {
                Some(read_hashes(&mut r)?)
            } else {
                None
            };
            (gk, pt)
        };
        let sketch = Arc::new(GenomeSketch {
            genome_kmers,
            pseudotax_tracked_nonused_kmers: pseudotax,
            file_name: meta.file_name.clone(),
            first_contig_name: meta.first_contig_name.clone(),
            c: self.c,
            k: self.k,
            gn_size: meta.gn_size,
            min_spacing: meta.min_spacing,
        });
        self.cache.lock().unwrap().insert(g, sketch.clone());
        Ok(sketch)
    }
}

/// Open a `.syl2db` file from a path.
pub fn open_file(path: &str) -> io::Result<TwoStageDb> {
    let r = BufReader::with_capacity(10_000_000, File::open(path)?);
    open(r)
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
        let mut buf = Vec::new();
        write_hashes(&mut buf, input);
        let mut r = &buf[..];
        let decoded = read_hashes(&mut r).unwrap();
        let mut expected = input.to_vec();
        expected.sort_unstable();
        assert_eq!(decoded, expected);
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
