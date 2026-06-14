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
//!   * A k-mer in two or more reps is ambiguous → the shared **pool**.
//!   * A k-mer in no rep but exactly one strain is owned by that strain.
//!   * A k-mer in no rep and several strains → the shared pool.
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
//!   * **Stage 2 (dense, loaded on demand):** each genome's full distinctive set
//!     is an independently Golomb–Rice-coded block at a known offset. Only the hit
//!     genomes' blocks are decoded (and cached across samples). The shared **pool**
//!     is the one exception — it is conserved across samples and is loaded once
//!     when the index is opened.
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
    self, read_genome_sketches_compressed, read_hashes, read_seq_sketch_compressed, read_string,
    read_uvarint, write_hashes, write_string, write_uvarint,
};
use crate::constants::*;
use crate::types::*;
use fxhash::FxHashMap;
use log::*;
use rayon::prelude::*;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::Path;
use std::sync::{Arc, Mutex};

const REFDB_MAGIC: &[u8; 4] = b"SYLR";
const REFDB_VERSION: u8 = 2;
const SKETCH_MAGIC: &[u8; 4] = b"SYLD"; // reference-Delta sample
const SKETCH_VERSION: u8 = 1;

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
    /// Hashes shared by ≥2 genomes (after rep preference), sorted and deduplicated.
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
        let owner = if a.rep_count == 1 {
            a.rep_id
        } else if a.rep_count >= 2 {
            POOL
        } else if a.strain_count == 1 {
            a.strain_id
        } else {
            POOL
        };
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
/// k-mers); `sparse_div <= 1` keeps all of them.
pub fn write_refdb<W: Write>(mut w: W, db: &RefDb, sparse_div: u64) -> io::Result<()> {
    let ng = db.genomes.len();
    let sp_ids = species_ids(&db.genomes);

    // Body: [sparse section][pool block][dense block 0][dense block 1]...
    let mut body: Vec<u8> = Vec::new();

    let sparse_off = body.len() as u64; // 0
    let mut sparse_count = vec![0usize; ng];
    for g in 0..ng {
        let mut cnt = 0usize;
        for &h in &db.distinctive[g] {
            if sparse_div <= 1 || h % sparse_div == 0 {
                body.extend_from_slice(&h.to_le_bytes());
                cnt += 1;
            }
        }
        sparse_count[g] = cnt;
    }

    let pool_off = body.len() as u64;
    write_hashes(&mut body, &db.pool)?;

    let mut dense_off = vec![0u64; ng];
    for g in 0..ng {
        dense_off[g] = body.len() as u64;
        write_hashes(&mut body, &db.distinctive[g])?;
    }

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
    let footer_comp = zstd::stream::encode_all(&footer[..], 9)?;
    let footer_offset = HEADER_LEN + body.len() as u64;

    w.write_all(REFDB_MAGIC)?;
    w.write_all(&[REFDB_VERSION])?;
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
    /// stage-1: sparse distinctive hash -> owning genome id.
    sparse_map: FxHashMap<u64, u32>,
    /// shared pool, loaded once.
    pool: Vec<u64>,
    pool_map: FxHashMap<u64, u32>,
    reader: Mutex<Box<dyn ReadSeek>>,
    cache: Mutex<FxHashMap<u32, Arc<Vec<u64>>>>,
}

/// Open a `.sylref`, loading stage 1 (sparse index) and the shared pool.
pub fn open_ref_index<R: Read + Seek + Send + 'static>(mut r: R) -> io::Result<RefIndex> {
    let mut hdr = [0u8; HEADER_LEN as usize];
    r.read_exact(&mut hdr)?;
    if &hdr[0..4] != REFDB_MAGIC {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "not a sylph reference DB"));
    }
    if hdr[4] != REFDB_VERSION {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "unsupported reference DB version"));
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
    let pool_offset = read_uvarint(&mut f)?;
    let _pool_domain = read_uvarint(&mut f)? as usize;
    let mut genomes = Vec::with_capacity(ng);
    for _ in 0..ng {
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
            dense_offset,
            dense_domain,
        });
    }

    // stage 1: sparse section (raw little-endian u64, grouped by genome)
    let total_sparse: usize = sparse_count.iter().sum();
    r.seek(SeekFrom::Start(sparse_offset))?;
    let mut sbuf = vec![0u8; total_sparse * 8];
    r.read_exact(&mut sbuf)?;
    let mut sparse_map: FxHashMap<u64, u32> =
        FxHashMap::with_capacity_and_hasher(total_sparse, Default::default());
    let mut pos = 0usize;
    for (g, &cnt) in sparse_count.iter().enumerate() {
        for _ in 0..cnt {
            let h = u64::from_le_bytes(sbuf[pos..pos + 8].try_into().unwrap());
            pos += 8;
            sparse_map.insert(h, g as u32);
        }
    }

    // shared pool (loaded once)
    r.seek(SeekFrom::Start(pool_offset))?;
    let pool = read_hashes(&mut r)?;
    let mut pool_map: FxHashMap<u64, u32> =
        FxHashMap::with_capacity_and_hasher(pool.len(), Default::default());
    for (i, &h) in pool.iter().enumerate() {
        pool_map.insert(h, i as u32);
    }

    Ok(RefIndex {
        c,
        k,
        sparse_div,
        fingerprint,
        genomes,
        sparse_map,
        pool,
        pool_map,
        reader: Mutex::new(Box::new(r)),
        cache: Mutex::new(FxHashMap::default()),
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
            if let Some(&g) = self.sparse_map.get(&h) {
                *counts.entry(g).or_insert(0) += 1;
            }
        }
        counts
            .into_iter()
            .filter(|&(_, c)| c >= SPARSE_MIN_HITS)
            .map(|(g, _)| g)
            .collect()
    }

    /// Stage-2: the genome's dense distinctive block, decoded on demand and cached.
    fn load_genome(&self, g: u32) -> io::Result<Arc<Vec<u64>>> {
        if let Some(a) = self.cache.lock().unwrap().get(&g) {
            return Ok(a.clone());
        }
        let off = self.genomes[g as usize].dense_offset;
        let arr = {
            let mut rd = self.reader.lock().unwrap();
            rd.seek(SeekFrom::Start(off))?;
            let mut r: &mut dyn ReadSeek = &mut **rd;
            read_hashes(&mut r)?
        };
        let arc = Arc::new(arr);
        self.cache.lock().unwrap().insert(g, arc.clone());
        Ok(arc)
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
        _ => Err(io::Error::new(io::ErrorKind::InvalidData, "bad subset scheme")),
    }
}

// --- compressing / decompressing a read sketch ------------------------------

/// Compress a read sketch against the reference index. Only the sample's hit
/// genomes' dense blocks are loaded; the pool is already resident in `idx`.
pub fn compress_seq<W: Write>(
    inner: W,
    sketch: &SequencesSketch,
    idx: &RefIndex,
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
        if let Some(&pidx) = idx.pool_map.get(&h) {
            pool_hits.push(pidx as u64);
        } else if let Some(&(g, i)) = map.get(&h) {
            per_genome.entry(g).or_default().push(i as u64);
        } else {
            novel.push(h);
        }
    }

    let mut payload = Vec::new();
    payload.extend_from_slice(SKETCH_MAGIC);
    payload.push(SKETCH_VERSION);
    payload.extend_from_slice(&idx.fingerprint.to_le_bytes());
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

    // hit genomes: sorted global ids, delta-coded (strains of a species are
    // contiguous, so a sample's hits cluster into small gaps)
    let mut hit_ids: Vec<u32> = per_genome.keys().copied().collect();
    hit_ids.sort_unstable();
    write_uvarint(&mut payload, hit_ids.len() as u64)?;
    let mut prev = 0u64;
    for &g in &hit_ids {
        write_uvarint(&mut payload, g as u64 - prev)?;
        prev = g as u64;
        let v = per_genome.get_mut(&g).unwrap();
        v.sort_unstable();
        encode_subset(&mut payload, v, idx.genomes[g as usize].dense_domain as u64)?;
    }

    // pool
    pool_hits.sort_unstable();
    write_uvarint(&mut payload, pool_hits.len() as u64)?;
    if !pool_hits.is_empty() {
        encode_subset(&mut payload, &pool_hits, idx.pool.len() as u64)?;
    }

    // novel hashes (Rice)
    write_hashes(&mut payload, &novel)?;

    // counts, in ascending-hash order (reproducible on decode)
    let mut keys: Vec<u64> = sketch.kmer_counts.keys().copied().collect();
    keys.sort_unstable();
    for h in &keys {
        write_uvarint(&mut payload, sketch.kmer_counts[h] as u64)?;
    }

    let mut enc = zstd::stream::write::Encoder::new(inner, ZSTD_LEVEL)?;
    enc.write_all(&payload)?;
    enc.finish()?;
    Ok(())
}

/// Decompress a reference-delta read sketch. Only the genomes referenced by the
/// sample (plus the resident pool) are loaded from the index.
pub fn decompress_seq<R: Read>(inner: R, idx: &RefIndex) -> io::Result<SequencesSketch> {
    let raw = zstd::stream::decode_all(inner)?;
    let mut r = &raw[..];
    let mut magic = [0u8; 4];
    r.read_exact(&mut magic)?;
    if &magic != SKETCH_MAGIC {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "not a reference-delta sketch"));
    }
    let mut ver = [0u8; 1];
    r.read_exact(&mut ver)?;
    if ver[0] != SKETCH_VERSION {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "unsupported reference-delta version"));
    }
    let mut fp = [0u8; 8];
    r.read_exact(&mut fp)?;
    if u64::from_le_bytes(fp) != idx.fingerprint {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "reference DB does not match the one used to compress this sample",
        ));
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

    hashes.sort_unstable();
    let mut kmer_counts = FxHashMap::with_capacity_and_hasher(hashes.len(), Default::default());
    for &h in &hashes {
        let count = read_uvarint(&mut r)? as u32;
        kmer_counts.insert(h, count);
    }

    Ok(SequencesSketch {
        kmer_counts,
        c,
        k,
        file_name,
        sample_name,
        paired: paired[0] != 0,
        mean_read_length,
    })
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

/// Load genome sketches from a `.syldb`, transparently handling the legacy
/// bincode and the compressed formats.
fn load_genome_sketches(path: &str) -> Vec<GenomeSketch> {
    let file = File::open(path).unwrap_or_else(|_| panic!("Could not open {}", path));
    let mut reader = BufReader::with_capacity(10_000_000, file);
    if compress::peek_is_compressed(&mut reader).unwrap_or(false) {
        read_genome_sketches_compressed(&mut reader)
            .unwrap_or_else(|_| panic!("{} is not a valid database sketch", path))
    } else {
        bincode::deserialize_from(&mut reader)
            .unwrap_or_else(|_| panic!("{} is not a valid database sketch", path))
    }
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
        map.insert(cols[0].trim().to_string(), (cols[1].trim().to_string(), is_rep));
    }
    map
}

pub fn run_ref_build(args: RefBuildArgs) {
    init_logger(args.trace);
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
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

    let mut sketches: Vec<GenomeSketch> = Vec::new();
    for f in &args.files {
        info!("Loading genome sketches from {}", f);
        sketches.extend(load_genome_sketches(f));
    }
    if sketches.is_empty() {
        error!("No genome sketches found; exiting");
        std::process::exit(1);
    }
    info!("Building k-mer dereplicated reference from {} genomes...", sketches.len());
    let db = build_refdb(&sketches, &taxonomy);
    let n_dist: usize = db.distinctive.iter().map(|d| d.len()).sum();
    let sparse_div = args.sparse_div.max(1);
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

    let out = format!("{}{}", args.output, REF_DB_SUFFIX);
    if let Some(parent) = Path::new(&out).parent() {
        std::fs::create_dir_all(parent).ok();
    }
    let w = BufWriter::new(File::create(&out).unwrap_or_else(|_| panic!("Could not create {}", out)));
    write_refdb(w, &db, sparse_div).unwrap_or_else(|e| panic!("Failed to write {}: {}", out, e));
    info!("Wrote reference database to {}", out);
}

/// Open a `.sylref` for querying/compression.
fn open_refdb_file(path: &str) -> RefIndex {
    let r = BufReader::with_capacity(
        10_000_000,
        File::open(path).unwrap_or_else(|_| panic!("Could not open {}", path)),
    );
    open_ref_index(r).unwrap_or_else(|e| panic!("{} is not a valid reference DB: {}", path, e))
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
    std::fs::create_dir_all(&args.output_dir)
        .expect("Could not create output directory; exiting");

    info!("Loading reference database {} (stage-1 sparse index + pool)", args.ref_db);
    let idx = open_refdb_file(&args.ref_db);

    let outdir = Path::new(&args.output_dir);
    let counter = Mutex::new(0usize);
    if args.decompress {
        args.files.par_iter().for_each(|f| {
            let r = BufReader::with_capacity(10_000_000, File::open(f).unwrap_or_else(|_| panic!("Could not open {}", f)));
            let sketch = decompress_seq(r, &idx).unwrap_or_else(|e| {
                error!("Failed to decompress {}: {}", f, e);
                std::process::exit(1);
            });
            let base = Path::new(f).file_name().unwrap().to_str().unwrap();
            let stem = base.strip_suffix(REF_SAMPLE_SUFFIX).unwrap_or(base);
            let out = outdir.join(format!("{}{}", stem, SAMPLE_FILE_SUFFIX));
            let mut w = BufWriter::new(File::create(&out).unwrap_or_else(|_| panic!("Could not create {:?}", out)));
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
            let w = BufWriter::new(File::create(&out).unwrap_or_else(|_| panic!("Could not create {:?}", out)));
            compress_seq(w, &sketch, &idx)
                .unwrap_or_else(|e| panic!("Failed to compress {}: {}", f, e));
            let mut c = counter.lock().unwrap();
            *c += 1;
            info!("Compressed {} -> {:?}", f, out);
        });
    }
    info!("Done ({} sample(s)).", *counter.lock().unwrap());
}
