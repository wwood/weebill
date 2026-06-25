//! Decompression of reference-delta `.sylspr` sample sketches.

use super::ref_build::RefIndex;
use super::{
    substituted_hash, window_fr, SCHEME_ABSENT_RICE, SCHEME_BITMASK, SCHEME_PRESENT_RICE,
    SKETCH_MAGIC, SKETCH_VERSION,
};
use crate::compress::{read_hashes, read_string, read_uvarint};
use crate::types::*;
use fxhash::FxHashMap;
use std::io::{self, Read};

// --- adaptive present/absent subset decoding ---------------------------------

pub(crate) fn decode_subset<R: Read>(r: &mut R, domain: u64) -> io::Result<Vec<u64>> {
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

// --- decompress_seq / decompress_seq_with_meta -------------------------------

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
                let (f, rr) = window_fr(&seq.contigs[ci].1, local, k);
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
