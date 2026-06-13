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
//! This is the same idea as sylph's profiling k-mer reassignment, but materialised
//! as a compression reference: each genome keeps only its *distinctive* hashes,
//! and conserved/shared hashes live in one pool so they remain referenceable
//! (dropping them would push a large fraction of reads back to full-price coding).
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
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::Path;
use std::sync::Mutex;

const REFDB_MAGIC: &[u8; 4] = b"SYLR";
const REFDB_VERSION: u8 = 1;
const SKETCH_MAGIC: &[u8; 4] = b"SYLD"; // reference-Delta sample
const SKETCH_VERSION: u8 = 1;

const SCHEME_BITMASK: u8 = 0;
const SCHEME_PRESENT_RICE: u8 = 1;
const SCHEME_ABSENT_RICE: u8 = 2;

/// zstd level for the read-sketch payload (matches the normal sketch format).
const ZSTD_LEVEL: i32 = 3;

/// Metadata for one reference genome, in the dereplicated build order.
#[derive(Clone, Debug, PartialEq)]
pub struct RefGenome {
    pub file_name: String,
    pub species: String,
    pub is_rep: bool,
}

/// A k-mer dereplicated reference database: per-genome distinctive hash sets plus
/// one shared pool, with genomes ordered so a species' strains are contiguous.
#[derive(Clone, Debug, PartialEq)]
pub struct RefDb {
    pub c: usize,
    pub k: usize,
    pub genomes: Vec<RefGenome>,
    /// `distinctive[g]` = sorted, deduplicated hashes owned uniquely by genome `g`.
    pub distinctive: Vec<Vec<u64>>,
    /// Hashes shared by ≥2 genomes (after rep preference), sorted and deduplicated.
    pub pool: Vec<u64>,
    /// Digest over the full reference contents, used to reject decoding a sample
    /// against the wrong DB. Computed once at build/load time (see `compute_fingerprint`).
    pub fingerprint: u64,
}

/// A digest over the *entire* reference contents (every distinctive and pool hash,
/// in order, with length separators). A sample is encoded as indices into these
/// arrays, so any difference that could change a decoded hash must change the
/// digest — mixing only lengths and boundary hashes would let a different DB with
/// matching boundaries pass and silently produce a wrong `SequencesSketch`.
fn compute_fingerprint(c: usize, k: usize, distinctive: &[Vec<u64>], pool: &[u64]) -> u64 {
    let mut h: u64 = 1469598103934665603; // FNV-1a offset
    let mut mix = |x: u64| {
        h ^= x;
        h = h.wrapping_mul(1099511628211);
    };
    mix(c as u64);
    mix(k as u64);
    mix(distinctive.len() as u64);
    for d in distinctive.iter() {
        mix(0x5359_4c44_4953_5400); // "SYLDIST\0" domain separator
        mix(d.len() as u64);
        for &x in d {
            mix(x);
        }
    }
    mix(0x5359_4c50_4f4f_4c00); // "SYLPOOL\0" domain separator
    mix(pool.len() as u64);
    for &x in pool {
        mix(x);
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

    let fingerprint = compute_fingerprint(c, k, &distinctive, &pool);
    RefDb {
        c,
        k,
        genomes,
        distinctive,
        pool,
        fingerprint,
    }
}

/// Lookup table: hash -> (owner genome id or `POOL`, index within that array).
pub struct RefLookup {
    map: FxHashMap<u64, (u32, u32)>,
}

impl RefDb {
    pub fn build_lookup(&self) -> RefLookup {
        let mut map: FxHashMap<u64, (u32, u32)> =
            FxHashMap::with_capacity_and_hasher(self.pool.len(), Default::default());
        for (g, d) in self.distinctive.iter().enumerate() {
            for (i, &h) in d.iter().enumerate() {
                map.insert(h, (g as u32, i as u32));
            }
        }
        for (i, &h) in self.pool.iter().enumerate() {
            map.insert(h, (POOL, i as u32));
        }
        RefLookup { map }
    }
}

// --- reference DB serialization (zstd-framed) -------------------------------

pub fn write_refdb<W: Write>(inner: W, db: &RefDb) -> io::Result<()> {
    let mut w = zstd::stream::write::Encoder::new(inner, 9)?;
    w.write_all(REFDB_MAGIC)?;
    w.write_all(&[REFDB_VERSION])?;
    write_uvarint(&mut w, db.c as u64)?;
    write_uvarint(&mut w, db.k as u64)?;
    write_uvarint(&mut w, db.genomes.len() as u64)?;
    for (g, genome) in db.genomes.iter().enumerate() {
        write_string(&mut w, &genome.file_name)?;
        write_string(&mut w, &genome.species)?;
        w.write_all(&[genome.is_rep as u8])?;
        write_hashes(&mut w, &db.distinctive[g])?;
    }
    write_hashes(&mut w, &db.pool)?;
    w.finish()?;
    Ok(())
}

pub fn read_refdb<R: Read>(inner: R) -> io::Result<RefDb> {
    let mut r = zstd::stream::read::Decoder::new(inner)?;
    let mut magic = [0u8; 4];
    r.read_exact(&mut magic)?;
    if &magic != REFDB_MAGIC {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "not a sylph reference DB"));
    }
    let mut ver = [0u8; 1];
    r.read_exact(&mut ver)?;
    if ver[0] != REFDB_VERSION {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "unsupported reference DB version"));
    }
    let c = read_uvarint(&mut r)? as usize;
    let k = read_uvarint(&mut r)? as usize;
    let ng = read_uvarint(&mut r)? as usize;
    let mut genomes = Vec::with_capacity(ng);
    let mut distinctive = Vec::with_capacity(ng);
    for _ in 0..ng {
        let file_name = read_string(&mut r)?;
        let species = read_string(&mut r)?;
        let mut rep = [0u8; 1];
        r.read_exact(&mut rep)?;
        genomes.push(RefGenome {
            file_name,
            species,
            is_rep: rep[0] != 0,
        });
        distinctive.push(read_hashes(&mut r)?);
    }
    let pool = read_hashes(&mut r)?;
    let fingerprint = compute_fingerprint(c, k, &distinctive, &pool);
    Ok(RefDb {
        c,
        k,
        genomes,
        distinctive,
        pool,
        fingerprint,
    })
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

/// Compress a read sketch against the reference DB. The DB (or one with the same
/// fingerprint) is required to decompress.
pub fn compress_seq<W: Write>(
    inner: W,
    sketch: &SequencesSketch,
    db: &RefDb,
    lookup: &RefLookup,
) -> io::Result<()> {
    let ng = db.genomes.len();
    let mut per_genome: Vec<Vec<u64>> = vec![Vec::new(); ng];
    let mut pool_hits: Vec<u64> = Vec::new();
    let mut novel: Vec<u64> = Vec::new();
    for &h in sketch.kmer_counts.keys() {
        match lookup.map.get(&h) {
            Some(&(POOL, idx)) => pool_hits.push(idx as u64),
            Some(&(g, idx)) => per_genome[g as usize].push(idx as u64),
            None => novel.push(h),
        }
    }

    let mut payload = Vec::new();
    payload.extend_from_slice(SKETCH_MAGIC);
    payload.push(SKETCH_VERSION);
    payload.extend_from_slice(&db.fingerprint.to_le_bytes());
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

    // hit genomes: sorted ids, delta-coded (strains of a species are contiguous,
    // so a sample's hits cluster into small gaps)
    let mut hit_ids: Vec<u32> = (0..ng as u32).filter(|&g| !per_genome[g as usize].is_empty()).collect();
    hit_ids.sort_unstable();
    write_uvarint(&mut payload, hit_ids.len() as u64)?;
    let mut prev = 0u64;
    for &g in &hit_ids {
        write_uvarint(&mut payload, g as u64 - prev)?;
        prev = g as u64;
        let idx = &mut per_genome[g as usize];
        idx.sort_unstable();
        encode_subset(&mut payload, idx, db.distinctive[g as usize].len() as u64)?;
    }

    // pool
    pool_hits.sort_unstable();
    write_uvarint(&mut payload, pool_hits.len() as u64)?;
    if !pool_hits.is_empty() {
        encode_subset(&mut payload, &pool_hits, db.pool.len() as u64)?;
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

/// Decompress a reference-delta read sketch using its reference DB.
pub fn decompress_seq<R: Read>(inner: R, db: &RefDb) -> io::Result<SequencesSketch> {
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
    if u64::from_le_bytes(fp) != db.fingerprint {
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
        let indices = decode_subset(&mut r, db.distinctive[g as usize].len() as u64)?;
        for i in indices {
            hashes.push(db.distinctive[g as usize][i as usize]);
        }
    }
    let npool = read_uvarint(&mut r)? as usize;
    if npool > 0 {
        let indices = decode_subset(&mut r, db.pool.len() as u64)?;
        for i in indices {
            hashes.push(db.pool[i as usize]);
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
    info!(
        "Reference: {} genomes, {} distinctive k-mers, {} shared-pool k-mers",
        db.genomes.len(),
        n_dist,
        db.pool.len()
    );

    let out = format!("{}{}", args.output, REF_DB_SUFFIX);
    if let Some(parent) = Path::new(&out).parent() {
        std::fs::create_dir_all(parent).ok();
    }
    let w = BufWriter::new(File::create(&out).unwrap_or_else(|_| panic!("Could not create {}", out)));
    write_refdb(w, &db).unwrap_or_else(|e| panic!("Failed to write {}: {}", out, e));
    info!("Wrote reference database to {}", out);
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

    info!("Loading reference database {}", args.ref_db);
    let db = {
        let r = BufReader::with_capacity(10_000_000, File::open(&args.ref_db).unwrap_or_else(|_| panic!("Could not open {}", args.ref_db)));
        read_refdb(r).unwrap_or_else(|e| panic!("{} is not a valid reference DB: {}", args.ref_db, e))
    };

    let outdir = Path::new(&args.output_dir);
    let counter = Mutex::new(0usize);
    if args.decompress {
        args.files.par_iter().for_each(|f| {
            let r = BufReader::with_capacity(10_000_000, File::open(f).unwrap_or_else(|_| panic!("Could not open {}", f)));
            let sketch = decompress_seq(r, &db).unwrap_or_else(|e| {
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
        info!("Building reference lookup...");
        let lookup = db.build_lookup();
        args.files.par_iter().for_each(|f| {
            let sketch = load_seq_sketch(f);
            let base = Path::new(f).file_name().unwrap().to_str().unwrap();
            let stem = base
                .strip_suffix(SAMPLE_FILE_SUFFIX)
                .or_else(|| base.strip_suffix(SAMPLE_COMP_FILE_SUFFIX))
                .unwrap_or(base);
            let out = outdir.join(format!("{}{}", stem, REF_SAMPLE_SUFFIX));
            let w = BufWriter::new(File::create(&out).unwrap_or_else(|_| panic!("Could not create {:?}", out)));
            compress_seq(w, &sketch, &db, &lookup)
                .unwrap_or_else(|e| panic!("Failed to compress {}: {}", f, e));
            let mut c = counter.lock().unwrap();
            *c += 1;
            info!("Compressed {} -> {:?}", f, out);
        });
    }
    info!("Done ({} sample(s)).", *counter.lock().unwrap());
}
