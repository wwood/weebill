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
//!      separate streams instead of interleaved pairs, so the small, repetitive
//!      counts group together and compress well.
//!   2. *Sort + delta + varint*: each hash set is sorted and stored as the
//!      gaps between successive values, encoded as LEB128 varints. The gaps are
//!      ~`u64::MAX / (c * n)` on average, requiring far fewer than 64 bits each.
//!   3. A generic `gzip`/DEFLATE pass on top of the varint stream.
//!
//! The on-disk container is a four-byte `SYLZ` magic, a version byte and a
//! sketch-type byte, followed by a gzip stream holding the payload. The magic
//! makes the format self-describing and detectable at read time independently
//! of the file extension: a legacy `bincode` file begins with a little-endian
//! length, which could only equal the `SYLZ` magic at ~1.5 billion entries —
//! impossible for a FracMinHash sketch — so the two formats never collide.

use crate::types::*;
use fxhash::FxHashMap;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};

const MAGIC: &[u8; 4] = b"SYLZ";
const VERSION: u8 = 1;
const TYPE_GENOME_DB: u8 = 1;
const TYPE_SEQ_SAMPLE: u8 = 2;

/// Returns `true` if the buffered reader is positioned at the start of a gzip
/// compressed sylph sketch, without consuming any bytes.
///
/// Detection is by the four-byte `SYLZ` magic that prefixes the file (ahead of
/// the gzip stream). The legacy `bincode` format begins with a little-endian
/// length, which could only match this signature at ~1.5 billion entries — far
/// beyond any possible FracMinHash sketch — so legacy sketches are never
/// misclassified and always fall back to the `bincode` reader.
pub fn peek_is_compressed<R: BufRead>(reader: &mut R) -> io::Result<bool> {
    let buf = reader.fill_buf()?;
    Ok(buf.len() >= MAGIC.len() && &buf[..MAGIC.len()] == MAGIC)
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

/// Sort + delta + varint encode a list of hash values. Order is not preserved
/// (hash sets are order-independent for sylph's purposes); duplicates, if any,
/// are preserved as zero gaps.
fn write_hashes<W: Write>(w: &mut W, hashes: &[u64]) -> io::Result<()> {
    write_uvarint(w, hashes.len() as u64)?;
    let mut sorted = hashes.to_vec();
    sorted.sort_unstable();
    let mut prev = 0u64;
    for &h in &sorted {
        write_uvarint(w, h.wrapping_sub(prev))?;
        prev = h;
    }
    Ok(())
}

fn read_hashes<R: Read>(r: &mut R) -> io::Result<Vec<u64>> {
    let n = read_uvarint(r)? as usize;
    let mut out = Vec::with_capacity(n);
    let mut prev = 0u64;
    for _ in 0..n {
        prev = prev.wrapping_add(read_uvarint(r)?);
        out.push(prev);
    }
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
    mut inner: W,
    sketches: &[GenomeSketch],
) -> io::Result<()> {
    write_header(&mut inner, TYPE_GENOME_DB)?;
    let enc = GzEncoder::new(inner, Compression::default());
    let mut w = BufWriter::new(enc);
    write_uvarint(&mut w, sketches.len() as u64)?;
    for s in sketches {
        write_genome_sketch(&mut w, s)?;
    }
    w.into_inner()
        .map_err(|e| io::Error::new(io::ErrorKind::Other, e.to_string()))?
        .finish()?;
    Ok(())
}

/// Read a database (a list of genome sketches) from the compressed format.
pub fn read_genome_sketches_compressed<R: Read>(mut inner: R) -> io::Result<Vec<GenomeSketch>> {
    read_and_check_header(&mut inner, TYPE_GENOME_DB)?;
    let mut r = BufReader::new(GzDecoder::new(inner));
    let n = read_uvarint(&mut r)? as usize;
    let mut out = Vec::with_capacity(n);
    for _ in 0..n {
        out.push(read_genome_sketch(&mut r)?);
    }
    Ok(out)
}

/// Write a sample (read) sketch in the compressed format.
pub fn write_seq_sketch_compressed<W: Write>(mut inner: W, s: &SequencesSketch) -> io::Result<()> {
    write_header(&mut inner, TYPE_SEQ_SAMPLE)?;
    let enc = GzEncoder::new(inner, Compression::default());
    let mut w = BufWriter::new(enc);
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

    // Struct-of-arrays: sorted, delta-coded keys followed by their counts.
    let mut pairs: Vec<(u64, u32)> = s.kmer_counts.iter().map(|(&k, &v)| (k, v)).collect();
    pairs.sort_unstable_by_key(|p| p.0);
    write_uvarint(&mut w, pairs.len() as u64)?;
    let mut prev = 0u64;
    for (k, _) in &pairs {
        write_uvarint(&mut w, k.wrapping_sub(prev))?;
        prev = *k;
    }
    for (_, count) in &pairs {
        write_uvarint(&mut w, *count as u64)?;
    }

    w.into_inner()
        .map_err(|e| io::Error::new(io::ErrorKind::Other, e.to_string()))?
        .finish()?;
    Ok(())
}

/// Read a sample (read) sketch from the compressed format.
pub fn read_seq_sketch_compressed<R: Read>(mut inner: R) -> io::Result<SequencesSketch> {
    read_and_check_header(&mut inner, TYPE_SEQ_SAMPLE)?;
    let mut r = BufReader::new(GzDecoder::new(inner));
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

    let n = read_uvarint(&mut r)? as usize;
    let mut keys = Vec::with_capacity(n);
    let mut prev = 0u64;
    for _ in 0..n {
        prev = prev.wrapping_add(read_uvarint(&mut r)?);
        keys.push(prev);
    }
    let mut kmer_counts =
        FxHashMap::with_capacity_and_hasher(n, Default::default());
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
