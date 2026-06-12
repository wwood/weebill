//! Compressed sketch (de)serialization.
//!
//! Sketches are sets of FracMinHash hash values: `u64`s drawn (approximately)
//! uniformly from `[0, u64::MAX / c)`. Stored naively (raw little-endian `u64`
//! per hash, as `bincode` does) they barely compress, because a set of uniform
//! integers is high entropy.
//!
//! This module stores the same information far more compactly by exploiting the
//! structure of the data:
//!   1. *Struct-of-arrays*: hashes and their multiplicities are written as two
//!      separate streams instead of interleaved pairs.
//!   2. *Sort + delta + Golomb–Rice*: each hash set is sorted and stored as the
//!      gaps between successive values. Because the hashes are uniform, the gaps
//!      are geometrically distributed, for which Golomb–Rice coding is
//!      near-optimal: each gap is split as `gap = (q << k) | r`, the quotient
//!      `q` is written in unary and the `k`-bit remainder `r` verbatim. The
//!      Rice parameter `k` is chosen per block to minimise the encoded size and
//!      stored in one byte. This reaches the information-theoretic entropy of
//!      the set without a separate generic-compression pass — which matters
//!      because that pass (gzip/DEFLATE) is both the slowest part of writing and
//!      buys nothing on an already near-random bitstream.
//!   3. Small per-hash multiplicities (read sketches only) are stored as
//!      varints in their own stream.
//!
//! The on-disk container is a four-byte `SYLZ` magic, a version byte and a
//! sketch-type byte, followed by the encoded payload. The magic makes the
//! format self-describing and detectable at read time independently of the file
//! extension: a legacy `bincode` file begins with a little-endian length, which
//! could only equal the `SYLZ` magic at ~1.5 billion entries — impossible for a
//! FracMinHash sketch — so the two formats never collide.

use crate::types::*;
use fxhash::FxHashMap;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};

const MAGIC: &[u8; 4] = b"SYLZ";
const VERSION: u8 = 2;
const TYPE_GENOME_DB: u8 = 1;
const TYPE_SEQ_SAMPLE: u8 = 2;

/// Returns `true` if the buffered reader is positioned at the start of a
/// compressed sylph sketch, without consuming any bytes.
///
/// Detection is by the four-byte `SYLZ` magic that prefixes the file. The
/// legacy `bincode` format begins with a little-endian length, which could only
/// match this signature at ~1.5 billion entries — far beyond any possible
/// FracMinHash sketch — so legacy sketches are never misclassified and always
/// fall back to the `bincode` reader.
pub fn peek_is_compressed<R: BufRead>(reader: &mut R) -> io::Result<bool> {
    let buf = reader.fill_buf()?;
    Ok(buf.len() >= MAGIC.len() && &buf[..MAGIC.len()] == MAGIC)
}

// --- bit I/O ----------------------------------------------------------------

/// MSB-first bit writer over a byte sink. Bits accumulate until whole bytes can
/// be flushed; `finish` pads the final partial byte with zero bits so that the
/// stream is byte-aligned and byte-oriented fields can follow.
struct BitWriter<'a, W: Write> {
    w: &'a mut W,
    acc: u128,
    nbits: u32,
}

impl<'a, W: Write> BitWriter<'a, W> {
    fn new(w: &'a mut W) -> Self {
        BitWriter { w, acc: 0, nbits: 0 }
    }

    #[inline]
    fn flush_bytes(&mut self) -> io::Result<()> {
        while self.nbits >= 8 {
            let byte = (self.acc >> (self.nbits - 8)) as u8;
            self.w.write_all(&[byte])?;
            self.nbits -= 8;
            self.acc &= (1u128 << self.nbits) - 1;
        }
        Ok(())
    }

    /// Write the low `n` bits of `val` (n <= 64).
    #[inline]
    fn write_bits(&mut self, val: u64, n: u32) -> io::Result<()> {
        if n == 0 {
            return Ok(());
        }
        let masked = if n >= 64 { val } else { val & ((1u64 << n) - 1) };
        self.acc = (self.acc << n) | masked as u128;
        self.nbits += n;
        self.flush_bytes()
    }

    #[inline]
    fn write_bit(&mut self, bit: u64) -> io::Result<()> {
        self.write_bits(bit & 1, 1)
    }

    /// Write `q` one-bits followed by a terminating zero-bit (unary code).
    #[inline]
    fn write_unary(&mut self, mut q: u64) -> io::Result<()> {
        while q >= 32 {
            self.write_bits(0xFFFF_FFFF, 32)?;
            q -= 32;
        }
        if q > 0 {
            self.write_bits((1u64 << q) - 1, q as u32)?;
        }
        self.write_bit(0)
    }

    /// Pad to the next byte boundary and flush.
    fn finish(mut self) -> io::Result<()> {
        if self.nbits > 0 {
            let pad = 8 - self.nbits;
            self.acc <<= pad;
            self.nbits += pad;
            self.flush_bytes()?;
        }
        Ok(())
    }
}

/// MSB-first bit reader over a byte source. Reads at most one byte ahead of the
/// bits consumed, so after `align` the underlying reader is positioned exactly
/// at the next byte boundary (matching `BitWriter::finish`'s padding).
struct BitReader<'a, R: Read> {
    r: &'a mut R,
    cur: u8,
    nleft: u32,
}

impl<'a, R: Read> BitReader<'a, R> {
    fn new(r: &'a mut R) -> Self {
        BitReader { r, cur: 0, nleft: 0 }
    }

    #[inline]
    fn read_bit(&mut self) -> io::Result<u64> {
        if self.nleft == 0 {
            let mut byte = [0u8; 1];
            self.r.read_exact(&mut byte)?;
            self.cur = byte[0];
            self.nleft = 8;
        }
        self.nleft -= 1;
        Ok(((self.cur >> self.nleft) & 1) as u64)
    }

    #[inline]
    fn read_bits(&mut self, n: u32) -> io::Result<u64> {
        let mut v = 0u64;
        for _ in 0..n {
            v = (v << 1) | self.read_bit()?;
        }
        Ok(v)
    }

    #[inline]
    fn read_unary(&mut self) -> io::Result<u64> {
        let mut q = 0u64;
        while self.read_bit()? == 1 {
            q += 1;
            if q > 4096 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "unary code too long; corrupt sketch",
                ));
            }
        }
        Ok(q)
    }

    /// Discard any buffered bits, returning to a byte boundary.
    fn align(&mut self) {
        self.nleft = 0;
    }
}

// --- primitive encoders -----------------------------------------------------

fn write_uvarint<W: Write>(w: &mut W, mut x: u64) -> io::Result<()> {
    loop {
        let mut byte = (x & 0x7f) as u8;
        x >>= 7;
        if x != 0 {
            byte |= 0x80;
        }
        w.write_all(&[byte])?;
        if x == 0 {
            return Ok(());
        }
    }
}

fn read_uvarint<R: Read>(r: &mut R) -> io::Result<u64> {
    let mut result: u64 = 0;
    let mut shift: u32 = 0;
    loop {
        let mut byte = [0u8; 1];
        r.read_exact(&mut byte)?;
        result |= ((byte[0] & 0x7f) as u64) << shift;
        if byte[0] & 0x80 == 0 {
            return Ok(result);
        }
        shift += 7;
        if shift >= 64 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "varint overflows u64",
            ));
        }
    }
}

fn write_string<W: Write>(w: &mut W, s: &str) -> io::Result<()> {
    write_uvarint(w, s.len() as u64)?;
    w.write_all(s.as_bytes())
}

fn read_string<R: Read>(r: &mut R) -> io::Result<String> {
    let len = read_uvarint(r)? as usize;
    let mut buf = vec![0u8; len];
    r.read_exact(&mut buf)?;
    String::from_utf8(buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn write_bool<W: Write>(w: &mut W, b: bool) -> io::Result<()> {
    w.write_all(&[b as u8])
}

fn read_bool<R: Read>(r: &mut R) -> io::Result<bool> {
    let mut byte = [0u8; 1];
    r.read_exact(&mut byte)?;
    Ok(byte[0] != 0)
}

fn write_f64<W: Write>(w: &mut W, x: f64) -> io::Result<()> {
    w.write_all(&x.to_le_bytes())
}

fn read_f64<R: Read>(r: &mut R) -> io::Result<f64> {
    let mut buf = [0u8; 8];
    r.read_exact(&mut buf)?;
    Ok(f64::from_le_bytes(buf))
}

#[inline]
fn floor_log2(x: u64) -> u32 {
    63 - x.leading_zeros()
}

/// Total number of bits the sorted block's gaps occupy for Rice parameter `k`.
fn rice_cost(sorted: &[u64], k: u32) -> u128 {
    let mut prev = 0u64;
    let mut bits: u128 = 0;
    for &h in sorted {
        let gap = h - prev;
        prev = h;
        bits += (gap >> k) as u128 + 1 + k as u128;
    }
    bits
}

/// Sort + delta + Golomb–Rice encode a list of hash values. Order is not
/// preserved (hash sets are order-independent for sylph's purposes); duplicates,
/// if any, are preserved as zero gaps. The block is byte-aligned on completion.
fn write_hashes<W: Write>(w: &mut W, hashes: &[u64]) -> io::Result<()> {
    write_uvarint(w, hashes.len() as u64)?;
    if hashes.is_empty() {
        return Ok(());
    }
    let mut sorted = hashes.to_vec();
    sorted.sort_unstable();

    // Choose the Rice parameter k. The mean gap ~ max/n; the optimum is near
    // floor(log2(mean)), so evaluate that and its neighbours and keep the best.
    let n = sorted.len() as u64;
    let max = *sorted.last().unwrap();
    let mean = (max / n).max(1);
    let k0 = floor_log2(mean);
    let mut best_k = k0;
    let mut best_cost = u128::MAX;
    for k in k0.saturating_sub(1)..=(k0 + 1).min(63) {
        let cost = rice_cost(&sorted, k);
        if cost < best_cost {
            best_cost = cost;
            best_k = k;
        }
    }
    w.write_all(&[best_k as u8])?;

    let mut bw = BitWriter::new(w);
    let mut prev = 0u64;
    for &h in &sorted {
        let gap = h - prev;
        prev = h;
        bw.write_unary(gap >> best_k)?;
        bw.write_bits(gap, best_k)?;
    }
    bw.finish()
}

fn read_hashes<R: Read>(r: &mut R) -> io::Result<Vec<u64>> {
    let n = read_uvarint(r)? as usize;
    if n == 0 {
        return Ok(Vec::new());
    }
    let mut kbuf = [0u8; 1];
    r.read_exact(&mut kbuf)?;
    let k = kbuf[0] as u32;

    let mut br = BitReader::new(r);
    let mut out = Vec::with_capacity(n);
    let mut prev = 0u64;
    for _ in 0..n {
        let q = br.read_unary()?;
        let rem = br.read_bits(k)?;
        let gap = (q << k) | rem;
        prev = prev.wrapping_add(gap);
        out.push(prev);
    }
    br.align();
    Ok(out)
}

// --- sketch encoders --------------------------------------------------------

fn write_header<W: Write>(w: &mut W, sketch_type: u8) -> io::Result<()> {
    w.write_all(MAGIC)?;
    w.write_all(&[VERSION, sketch_type])
}

fn read_and_check_header<R: Read>(r: &mut R, expected_type: u8) -> io::Result<()> {
    let mut magic = [0u8; 4];
    r.read_exact(&mut magic)?;
    if &magic != MAGIC {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "not a compressed sylph sketch (bad magic)",
        ));
    }
    let mut meta = [0u8; 2];
    r.read_exact(&mut meta)?;
    if meta[0] != VERSION {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("unsupported compressed sketch version {}", meta[0]),
        ));
    }
    if meta[1] != expected_type {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "compressed sketch is of an unexpected type",
        ));
    }
    Ok(())
}

fn write_genome_sketch<W: Write>(w: &mut W, s: &GenomeSketch) -> io::Result<()> {
    write_string(w, &s.file_name)?;
    write_string(w, &s.first_contig_name)?;
    write_uvarint(w, s.c as u64)?;
    write_uvarint(w, s.k as u64)?;
    write_uvarint(w, s.gn_size as u64)?;
    write_uvarint(w, s.min_spacing as u64)?;
    write_hashes(w, &s.genome_kmers)?;
    match &s.pseudotax_tracked_nonused_kmers {
        Some(kmers) => {
            write_bool(w, true)?;
            write_hashes(w, kmers)?;
        }
        None => write_bool(w, false)?,
    }
    Ok(())
}

fn read_genome_sketch<R: Read>(r: &mut R) -> io::Result<GenomeSketch> {
    let file_name = read_string(r)?;
    let first_contig_name = read_string(r)?;
    let c = read_uvarint(r)? as usize;
    let k = read_uvarint(r)? as usize;
    let gn_size = read_uvarint(r)? as usize;
    let min_spacing = read_uvarint(r)? as usize;
    let genome_kmers = read_hashes(r)?;
    let pseudotax_tracked_nonused_kmers = if read_bool(r)? {
        Some(read_hashes(r)?)
    } else {
        None
    };
    Ok(GenomeSketch {
        genome_kmers,
        pseudotax_tracked_nonused_kmers,
        file_name,
        first_contig_name,
        c,
        k,
        gn_size,
        min_spacing,
    })
}

/// Write a database (a list of genome sketches) in the compressed format.
pub fn write_genome_sketches_compressed<W: Write>(
    inner: W,
    sketches: &[GenomeSketch],
) -> io::Result<()> {
    let mut w = BufWriter::new(inner);
    write_header(&mut w, TYPE_GENOME_DB)?;
    write_uvarint(&mut w, sketches.len() as u64)?;
    for s in sketches {
        write_genome_sketch(&mut w, s)?;
    }
    w.flush()
}

/// Read a database (a list of genome sketches) from the compressed format.
pub fn read_genome_sketches_compressed<R: Read>(inner: R) -> io::Result<Vec<GenomeSketch>> {
    let mut r = BufReader::with_capacity(1 << 22, inner);
    read_and_check_header(&mut r, TYPE_GENOME_DB)?;
    let n = read_uvarint(&mut r)? as usize;
    let mut out = Vec::with_capacity(n);
    for _ in 0..n {
        out.push(read_genome_sketch(&mut r)?);
    }
    Ok(out)
}

/// Write a sample (read) sketch in the compressed format.
pub fn write_seq_sketch_compressed<W: Write>(inner: W, s: &SequencesSketch) -> io::Result<()> {
    let mut w = BufWriter::new(inner);
    write_header(&mut w, TYPE_SEQ_SAMPLE)?;
    write_uvarint(&mut w, s.c as u64)?;
    write_uvarint(&mut w, s.k as u64)?;
    write_string(&mut w, &s.file_name)?;
    match &s.sample_name {
        Some(name) => {
            write_bool(&mut w, true)?;
            write_string(&mut w, name)?;
        }
        None => write_bool(&mut w, false)?,
    }
    write_bool(&mut w, s.paired)?;
    write_f64(&mut w, s.mean_read_length)?;

    // Struct-of-arrays: Rice-coded sorted keys, then their counts as varints.
    let mut pairs: Vec<(u64, u32)> = s.kmer_counts.iter().map(|(&k, &v)| (k, v)).collect();
    pairs.sort_unstable_by_key(|p| p.0);
    let keys: Vec<u64> = pairs.iter().map(|p| p.0).collect();
    write_hashes(&mut w, &keys)?;
    for (_, count) in &pairs {
        write_uvarint(&mut w, *count as u64)?;
    }

    w.flush()
}

/// Read a sample (read) sketch from the compressed format.
pub fn read_seq_sketch_compressed<R: Read>(inner: R) -> io::Result<SequencesSketch> {
    let mut r = BufReader::with_capacity(1 << 22, inner);
    read_and_check_header(&mut r, TYPE_SEQ_SAMPLE)?;
    let c = read_uvarint(&mut r)? as usize;
    let k = read_uvarint(&mut r)? as usize;
    let file_name = read_string(&mut r)?;
    let sample_name = if read_bool(&mut r)? {
        Some(read_string(&mut r)?)
    } else {
        None
    };
    let paired = read_bool(&mut r)?;
    let mean_read_length = read_f64(&mut r)?;

    let keys = read_hashes(&mut r)?;
    let mut kmer_counts = FxHashMap::with_capacity_and_hasher(keys.len(), Default::default());
    for &key in &keys {
        let count = read_uvarint(&mut r)? as u32;
        kmer_counts.insert(key, count);
    }

    Ok(SequencesSketch {
        kmer_counts,
        c,
        k,
        file_name,
        sample_name,
        paired,
        mean_read_length,
    })
}
