//! Recovering the genome order a lost `.sylref` was built with.
//!
//! A reference's genome ids are reproducible — they come from the build order
//! (species, representatives first, then file name) — but the *ownership* of a
//! k-mer contested by fewer than `--pool-min-genomes` same-tier genomes is
//! decided by the order those genomes were fed to the build, which is whatever
//! order the sketches happened to sit in their `.syldb`. Rebuilding from
//! re-sketched genomes therefore reproduces the genome ids but not necessarily
//! the ownership, and a sample compressed against the original (`.sylspr`) then
//! refuses to decode against the rebuild.
//!
//! The order is recoverable because the evidence is still in the samples:
//!
//!   * **`pairs`** routes every genome's k-mers exactly as a build does and
//!     records, instead of a reference, the *contested groups*: for each k-mer
//!     class owned by 2 or more same-tier genomes below the pool threshold, which
//!     genomes are in the running. The winner is the group member that came first
//!     in the build — one unknown bit (or, for larger groups, one choice) per
//!     group, and nothing else about the order matters.
//!   * **`probe`** reads each `.sylspr`'s novel-hash section — the hashes no
//!     reference genome explained — which needs no reference at all, because the
//!     header records the byte lengths of the sections in front of it. A
//!     contested group's k-mer turning up as *novel* in a sample proves the
//!     group's winner is a genome that sample did not hit; combined with the
//!     sample's genome content (from its profile) that identifies the winner.
//!   * **`fingerprint`** scores a candidate order without building a reference,
//!     so a residual ambiguity can be searched cheaply. `len`, `first` and `last`
//!     of each genome's distinctive array — all the fingerprint sees — are
//!     recoverable from the group data, since min and max over a union are the
//!     min and max of the parts.
//!
//! The answer is then checked the only way that counts: build with
//! `ref-build --genome-order` and require the fingerprint to match the one every
//! `.sylspr` records.

use crate::cmdline::{RefRecoverArgs, RefRecoverMode};
use crate::compress::{read_hashes, read_uvarint, write_uvarint};
use crate::constants::*;
use crate::refdelta::ref_build::{route_genomes, Routed};
use fxhash::FxHashMap;
use log::*;
use rayon::prelude::*;
use std::fs::File;
use std::io::{self, BufWriter, Read, Write};
use std::path::Path;

const GRAPH_MAGIC: &[u8; 4] = b"SYLG";
const GRAPH_VERSION: u8 = 1;

/// The largest contested group this tool represents. A group has at most
/// `--pool-min-genomes - 1` members (anything larger goes to the shared pool), so
/// this bounds the pool threshold rather than the reference size.
const MAX_GROUP: usize = 8;

// --- the contested-group graph ----------------------------------------------

/// What a genome owns no matter the order: the count, and the extremes of the
/// hash set, which is all the fingerprint looks at.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct GenomeFixed {
    pub file_name: String,
    pub species: String,
    pub is_rep: bool,
    pub count: u64,
    pub min: u64,
    pub max: u64,
}

/// One class of k-mers contested by `members` (genome ids, ascending). All of
/// them go to whichever member came first in the build order.
#[derive(Clone, Debug, PartialEq)]
pub struct Group {
    pub members: Vec<u32>,
    pub count: u64,
    pub min: u64,
    pub max: u64,
    /// A sample of the group's hashes, for `probe` to look for in novel sections.
    pub probe: Vec<u64>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Graph {
    pub c: usize,
    pub k: usize,
    pub pool_min_genomes: u32,
    pub genomes: Vec<GenomeFixed>,
    pub groups: Vec<Group>,
    pub pool_count: u64,
    pub pool_min: u64,
    pub pool_max: u64,
}

/// Running min/max/count over a hash set, so a genome's or group's extremes can
/// be accumulated without holding the hashes.
#[derive(Clone, Copy, Debug)]
struct Extent {
    count: u64,
    min: u64,
    max: u64,
}

impl Default for Extent {
    fn default() -> Self {
        Extent {
            count: 0,
            min: u64::MAX,
            max: 0,
        }
    }
}

impl Extent {
    fn add(&mut self, h: u64) {
        self.count += 1;
        self.min = self.min.min(h);
        self.max = self.max.max(h);
    }
    fn merge(&mut self, o: &Extent) {
        if o.count == 0 {
            return;
        }
        self.count += o.count;
        self.min = self.min.min(o.min);
        self.max = self.max.max(o.max);
    }
}

/// Per-k-mer tally, mirroring `ref_build`'s but keeping every contender rather
/// than only the winner: which genomes win is exactly what is unknown here.
struct GroupAccum {
    rep: heapless_set::Set,
    strain: heapless_set::Set,
}

/// A tiny ascending set of genome ids that stops growing once the group is bound
/// for the shared pool, so a k-mer in thousands of genomes costs nothing.
mod heapless_set {
    use super::MAX_GROUP;

    #[derive(Clone, Copy)]
    pub struct Set {
        ids: [u32; MAX_GROUP],
        pub len: u8,
        pub count: u32,
    }

    impl Default for Set {
        fn default() -> Self {
            Set {
                ids: [0; MAX_GROUP],
                len: 0,
                count: 0,
            }
        }
    }

    impl Set {
        pub fn add(&mut self, id: u32) {
            self.count += 1;
            if (self.len as usize) >= MAX_GROUP {
                return;
            }
            // insert ascending, ignoring duplicates (a genome is deduped upstream)
            let mut i = 0usize;
            while i < self.len as usize && self.ids[i] < id {
                i += 1;
            }
            if i < self.len as usize && self.ids[i] == id {
                return;
            }
            let n = self.len as usize;
            self.ids.copy_within(i..n, i + 1);
            self.ids[i] = id;
            self.len += 1;
        }
        pub fn slice(&self) -> &[u32] {
            &self.ids[..self.len as usize]
        }
    }
}

/// Build the contested-group graph. Routing is a build's own pass 1, so the
/// genome ids here are exactly a reference's, whatever order the sketches are in.
pub fn build_graph(routed: &Routed, pool_min_genomes: u32) -> Graph {
    let ng = routed.genomes.len();
    let per_shard: Vec<ShardGroups> = (0..routed.partitions)
        .into_par_iter()
        .map(|pi| {
            shard_groups(routed, pi, pool_min_genomes, ng)
                .unwrap_or_else(|e| panic!("Failed to read scratch shard {}: {}", pi, e))
        })
        .collect();

    let mut fixed = vec![Extent::default(); ng];
    let mut groups: FxHashMap<Vec<u32>, (Extent, Vec<u64>)> = FxHashMap::default();
    let mut pool = Extent::default();
    for (f, g, p) in per_shard {
        for (i, e) in f.iter().enumerate() {
            fixed[i].merge(e);
        }
        for (members, (e, probe)) in g {
            let slot = groups
                .entry(members)
                .or_insert_with(|| (Extent::default(), Vec::new()));
            slot.0.merge(&e);
            slot.1.extend_from_slice(&probe);
        }
        pool.merge(&p);
    }

    let mut groups: Vec<Group> = groups
        .into_iter()
        .map(|(members, (e, mut probe))| {
            probe.sort_unstable();
            probe.dedup();
            probe.truncate(PROBE_PER_GROUP);
            Group {
                members,
                count: e.count,
                min: e.min,
                max: e.max,
                probe,
            }
        })
        .collect();
    // Stable order so the graph file is reproducible run to run.
    groups.sort_by(|a, b| a.members.cmp(&b.members));

    Graph {
        c: routed.c,
        k: routed.k,
        pool_min_genomes,
        genomes: routed
            .genomes
            .iter()
            .enumerate()
            .map(|(i, g)| GenomeFixed {
                file_name: g.file_name.clone(),
                species: g.species.clone(),
                is_rep: g.is_rep,
                count: fixed[i].count,
                min: fixed[i].min,
                max: fixed[i].max,
            })
            .collect(),
        groups,
        pool_count: pool.count,
        pool_min: pool.min,
        pool_max: pool.max,
    }
}

/// How many of a group's hashes to keep for probing. The probe map is the one
/// structure that scales with the reference rather than with the genome count,
/// so it is capped per group rather than kept whole.
const PROBE_PER_GROUP: usize = 16;

/// One shard's contribution: what each genome owns outright, the contested
/// groups it saw (with a few of their hashes for probing), and the pool.
type ShardGroups = (Vec<Extent>, FxHashMap<Vec<u32>, (Extent, Vec<u64>)>, Extent);

fn shard_groups(
    routed: &Routed,
    pi: usize,
    pool_min_genomes: u32,
    ng: usize,
) -> io::Result<ShardGroups> {
    use std::io::BufReader;
    let mut r = BufReader::with_capacity(1 << 20, File::open(routed.shard_path(pi))?);
    let mut acc: FxHashMap<u64, GroupAccum> = FxHashMap::default();
    // Genomes arrive by file id; the graph is expressed in genome ids so that it
    // does not depend on the order this particular run happened to see.
    while let Some(fid) = crate::refdelta::ref_build::read_shard_fid(&mut r)? {
        let mut rep = [0u8; 1];
        r.read_exact(&mut rep)?;
        let is_rep = rep[0] != 0;
        let gid = routed.remap[fid as usize];
        for h in read_hashes(&mut r)? {
            let e = acc.entry(h).or_insert_with(|| GroupAccum {
                rep: Default::default(),
                strain: Default::default(),
            });
            if is_rep {
                e.rep.add(gid);
            } else {
                e.strain.add(gid);
            }
        }
    }

    let mut fixed = vec![Extent::default(); ng];
    let mut groups: FxHashMap<Vec<u32>, (Extent, Vec<u64>)> = FxHashMap::default();
    let mut pool = Extent::default();
    for (h, a) in acc {
        // Representatives outrank strains: a k-mer in any representative is
        // settled among the representatives alone.
        let set = if a.rep.count > 0 { &a.rep } else { &a.strain };
        if set.count >= pool_min_genomes {
            pool.add(h);
            continue;
        }
        let members = set.slice();
        if members.len() == 1 {
            fixed[members[0] as usize].add(h);
            continue;
        }
        let slot = groups
            .entry(members.to_vec())
            .or_insert_with(|| (Extent::default(), Vec::new()));
        slot.0.add(h);
        if slot.1.len() < PROBE_PER_GROUP {
            slot.1.push(h);
        }
    }
    Ok((fixed, groups, pool))
}

// --- graph file I/O ----------------------------------------------------------

fn write_string<W: Write>(w: &mut W, s: &str) -> io::Result<()> {
    write_uvarint(w, s.len() as u64)?;
    w.write_all(s.as_bytes())
}

fn read_string<R: Read>(r: &mut R) -> io::Result<String> {
    let n = read_uvarint(r)? as usize;
    let mut v = vec![0u8; n];
    r.read_exact(&mut v)?;
    String::from_utf8(v).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e.to_string()))
}

pub fn write_graph<W: Write>(w: W, g: &Graph) -> io::Result<()> {
    let mut buf: Vec<u8> = Vec::new();
    buf.extend_from_slice(GRAPH_MAGIC);
    buf.push(GRAPH_VERSION);
    write_uvarint(&mut buf, g.c as u64)?;
    write_uvarint(&mut buf, g.k as u64)?;
    write_uvarint(&mut buf, g.pool_min_genomes as u64)?;
    write_uvarint(&mut buf, g.pool_count)?;
    buf.extend_from_slice(&g.pool_min.to_le_bytes());
    buf.extend_from_slice(&g.pool_max.to_le_bytes());
    write_uvarint(&mut buf, g.genomes.len() as u64)?;
    for gen in &g.genomes {
        write_string(&mut buf, &gen.file_name)?;
        write_string(&mut buf, &gen.species)?;
        buf.push(gen.is_rep as u8);
        write_uvarint(&mut buf, gen.count)?;
        buf.extend_from_slice(&gen.min.to_le_bytes());
        buf.extend_from_slice(&gen.max.to_le_bytes());
    }
    write_uvarint(&mut buf, g.groups.len() as u64)?;
    for grp in &g.groups {
        write_uvarint(&mut buf, grp.members.len() as u64)?;
        let mut prev = 0u64;
        for &m in &grp.members {
            write_uvarint(&mut buf, m as u64 - prev)?;
            prev = m as u64;
        }
        write_uvarint(&mut buf, grp.count)?;
        buf.extend_from_slice(&grp.min.to_le_bytes());
        buf.extend_from_slice(&grp.max.to_le_bytes());
        write_uvarint(&mut buf, grp.probe.len() as u64)?;
        for &h in &grp.probe {
            buf.extend_from_slice(&h.to_le_bytes());
        }
    }
    let mut enc = zstd::stream::write::Encoder::new(w, 3)?;
    enc.include_checksum(true)?;
    enc.write_all(&buf)?;
    enc.finish()?;
    Ok(())
}

pub fn read_graph<R: Read>(r: R) -> io::Result<Graph> {
    let buf = zstd::stream::decode_all(r)?;
    let mut r = &buf[..];
    let mut magic = [0u8; 4];
    r.read_exact(&mut magic)?;
    if &magic != GRAPH_MAGIC {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "not a weebill recovery graph",
        ));
    }
    let mut ver = [0u8; 1];
    r.read_exact(&mut ver)?;
    if ver[0] != GRAPH_VERSION {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "recovery graph is version {}, expected {}; rebuild it with `ref-recover pairs`",
                ver[0], GRAPH_VERSION
            ),
        ));
    }
    let c = read_uvarint(&mut r)? as usize;
    let k = read_uvarint(&mut r)? as usize;
    let pool_min_genomes = read_uvarint(&mut r)? as u32;
    let pool_count = read_uvarint(&mut r)?;
    let mut b = [0u8; 8];
    r.read_exact(&mut b)?;
    let pool_min = u64::from_le_bytes(b);
    r.read_exact(&mut b)?;
    let pool_max = u64::from_le_bytes(b);
    let ng = read_uvarint(&mut r)? as usize;
    let mut genomes = Vec::with_capacity(ng);
    for _ in 0..ng {
        let file_name = read_string(&mut r)?;
        let species = read_string(&mut r)?;
        let mut rep = [0u8; 1];
        r.read_exact(&mut rep)?;
        let count = read_uvarint(&mut r)?;
        r.read_exact(&mut b)?;
        let min = u64::from_le_bytes(b);
        r.read_exact(&mut b)?;
        let max = u64::from_le_bytes(b);
        genomes.push(GenomeFixed {
            file_name,
            species,
            is_rep: rep[0] != 0,
            count,
            min,
            max,
        });
    }
    let n_groups = read_uvarint(&mut r)? as usize;
    let mut groups = Vec::with_capacity(n_groups);
    for _ in 0..n_groups {
        let nm = read_uvarint(&mut r)? as usize;
        let mut members = Vec::with_capacity(nm);
        let mut prev = 0u64;
        for _ in 0..nm {
            prev += read_uvarint(&mut r)?;
            members.push(prev as u32);
        }
        let count = read_uvarint(&mut r)?;
        r.read_exact(&mut b)?;
        let min = u64::from_le_bytes(b);
        r.read_exact(&mut b)?;
        let max = u64::from_le_bytes(b);
        let np = read_uvarint(&mut r)? as usize;
        let mut probe = Vec::with_capacity(np);
        for _ in 0..np {
            r.read_exact(&mut b)?;
            probe.push(u64::from_le_bytes(b));
        }
        groups.push(Group {
            members,
            count,
            min,
            max,
            probe,
        });
    }
    Ok(Graph {
        c,
        k,
        pool_min_genomes,
        genomes,
        groups,
        pool_count,
        pool_min,
        pool_max,
    })
}

// --- analytic fingerprint ----------------------------------------------------

/// The fingerprint the reference would carry if built with `order` (genome id ->
/// its position in the build), without building it. Mirrors `ref_build`'s
/// `fingerprint` exactly: only each genome's k-mer count and hash extremes, plus
/// the (order-independent) pool, feed it.
pub fn fingerprint_for_order(g: &Graph, rank: &[u32]) -> u64 {
    let ng = g.genomes.len();
    let mut count: Vec<u64> = g.genomes.iter().map(|x| x.count).collect();
    let mut lo: Vec<u64> = g.genomes.iter().map(|x| x.min).collect();
    let mut hi: Vec<u64> = g.genomes.iter().map(|x| x.max).collect();
    for grp in &g.groups {
        // The winner is the member that comes first in this candidate order.
        let winner = *grp
            .members
            .iter()
            .min_by_key(|&&m| rank[m as usize])
            .expect("a contested group always has members") as usize;
        count[winner] += grp.count;
        lo[winner] = lo[winner].min(grp.min);
        hi[winner] = hi[winner].max(grp.max);
    }

    let mut h: u64 = 1469598103934665603;
    let mut mix = |x: u64| {
        h ^= x;
        h = h.wrapping_mul(1099511628211);
    };
    mix(g.c as u64);
    mix(g.k as u64);
    mix(ng as u64);
    mix(g.pool_count);
    for i in 0..ng {
        mix(i as u64);
        mix(count[i]);
        if count[i] > 0 {
            mix(lo[i]);
            mix(hi[i]);
        }
    }
    if g.pool_count > 0 {
        mix(g.pool_min);
        mix(g.pool_max);
    }
    h
}

// --- CLI ---------------------------------------------------------------------

pub fn run_ref_recover(args: RefRecoverArgs) {
    super::init_logger(args.trace);
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads.max(1))
        .build_global()
        .ok();
    match args.mode {
        RefRecoverMode::Pairs => run_pairs(args),
        RefRecoverMode::Probe => run_probe(args),
        RefRecoverMode::Fingerprint => run_fingerprint(args),
        RefRecoverMode::Search => run_search(args),
    }
}

fn run_pairs(args: RefRecoverArgs) {
    if args.files.is_empty() {
        error!("ref-recover pairs needs the genome database sketches (*.syldb) the reference was built from; exiting");
        std::process::exit(1);
    }
    let pool_min_genomes = args.pool_min_genomes.max(2);
    if pool_min_genomes as usize > MAX_GROUP + 1 {
        error!(
            "--pool-min-genomes {} would allow contested groups of {} genomes; this tool represents at most {}",
            pool_min_genomes,
            pool_min_genomes - 1,
            MAX_GROUP
        );
        std::process::exit(1);
    }
    let taxonomy = match &args.taxonomy {
        Some(p) => crate::refdelta::ref_build::parse_taxonomy_file(p),
        None => FxHashMap::default(),
    };
    let out = args.output.clone().unwrap_or_else(|| {
        error!("ref-recover pairs needs -o/--output; exiting");
        std::process::exit(1);
    });
    let routed = route_genomes(
        &args.files,
        &taxonomy,
        None, // the graph is order-independent by construction
        args.partitions.max(1),
        args.tmp_dir.as_deref(),
        &out,
        "sylref_recover",
    );
    let graph = build_graph(&routed, pool_min_genomes);
    routed.cleanup();

    let contested: u64 = graph.groups.iter().map(|g| g.count).sum();
    let pairs = graph.groups.iter().filter(|g| g.members.len() == 2).count();
    info!(
        "{} genomes, {} contested groups ({} of them pairs) over {} k-mers, {} pool k-mers",
        graph.genomes.len(),
        graph.groups.len(),
        pairs,
        contested,
        graph.pool_count
    );
    info!(
        "Ambiguity: {} unknown group(s); a genome touches at most {} of them",
        graph.groups.len(),
        max_group_degree(&graph)
    );

    if let Some(parent) = Path::new(&out).parent() {
        std::fs::create_dir_all(parent).ok();
    }
    let f = File::create(&out).unwrap_or_else(|e| panic!("Could not create {}: {}", out, e));
    write_graph(BufWriter::new(f), &graph)
        .unwrap_or_else(|e| panic!("Failed to write {}: {}", out, e));
    info!("Wrote recovery graph to {}", out);

    if let Some(tsv) = &args.tsv {
        write_graph_tsv(tsv, &graph).unwrap_or_else(|e| panic!("Failed to write {}: {}", tsv, e));
        info!("Wrote graph TSV to {}", tsv);
    }
}

fn max_group_degree(g: &Graph) -> usize {
    let mut deg = vec![0usize; g.genomes.len()];
    for grp in &g.groups {
        for &m in &grp.members {
            deg[m as usize] += 1;
        }
    }
    deg.into_iter().max().unwrap_or(0)
}

/// Two record types, keyed by the first column, so one file carries both the
/// genome table (needed to write a `--genome-order`) and the unknowns.
fn write_graph_tsv(path: &str, g: &Graph) -> io::Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    writeln!(
        w,
        "#record\tgenome_id\tfile_name\tspecies\tis_rep\tfixed_kmers"
    )?;
    writeln!(
        w,
        "#record\tgroup_id\tmembers\tcontested_kmers\tprobe_hashes"
    )?;
    for (i, gen) in g.genomes.iter().enumerate() {
        writeln!(
            w,
            "genome\t{}\t{}\t{}\t{}\t{}",
            i, gen.file_name, gen.species, gen.is_rep, gen.count
        )?;
    }
    for (i, grp) in g.groups.iter().enumerate() {
        let members: Vec<String> = grp.members.iter().map(|m| m.to_string()).collect();
        writeln!(
            w,
            "group\t{}\t{}\t{}\t{}",
            i,
            members.join(","),
            grp.count,
            grp.probe.len()
        )?;
    }
    w.flush()
}

fn load_graph(path: &str) -> Graph {
    let f = File::open(path).unwrap_or_else(|e| panic!("Could not open {}: {}", path, e));
    read_graph(std::io::BufReader::new(f))
        .unwrap_or_else(|e| panic!("{} is not a valid recovery graph: {}", path, e))
}

fn run_probe(args: RefRecoverArgs) {
    let graph_path = args.graph.clone().unwrap_or_else(|| {
        error!("ref-recover probe needs --graph from `ref-recover pairs`; exiting");
        std::process::exit(1);
    });
    let graph = load_graph(&graph_path);

    // hash -> group id, as one sorted array: the probe set is bounded per group,
    // so this stays proportional to the number of unknowns, not to the reference.
    let mut probe: Vec<(u64, u32)> = Vec::new();
    for (gi, grp) in graph.groups.iter().enumerate() {
        for &h in &grp.probe {
            probe.push((h, gi as u32));
        }
    }
    probe.sort_unstable();
    info!(
        "Probing {} sample(s) for {} hashes across {} contested groups",
        args.files.len(),
        probe.len(),
        graph.groups.len()
    );

    let rows: Vec<String> = args
        .files
        .par_iter()
        .filter_map(|path| match probe_one(path, &probe) {
            Ok(Some(rows)) => Some(rows),
            Ok(None) => None,
            Err(e) => {
                warn!("{}: {}", path, e);
                None
            }
        })
        .collect();

    let mut w: Box<dyn Write> = match &args.output {
        Some(p) => Box::new(BufWriter::new(
            File::create(p).unwrap_or_else(|e| panic!("Could not create {}: {}", p, e)),
        )),
        None => Box::new(BufWriter::new(std::io::stdout())),
    };
    writeln!(
        w,
        "sample_file\treference_fingerprint\thit_genomes\tnovel_hashes\tgroup_id\tnovel_group_hashes"
    )
    .ok();
    for r in &rows {
        write!(w, "{}", r).ok();
    }
    w.flush().ok();
}

/// Read one `.sylspr`'s novel hashes and report which contested groups they hit.
/// Nothing here needs the reference: the header records the byte lengths of the
/// hit and pool sections, so the novel section can be reached by skipping them.
fn probe_one(path: &str, probe: &[(u64, u32)]) -> io::Result<Option<String>> {
    if !path.ends_with(REF_SAMPLE_SUFFIX) {
        // A collection usually holds plain sketches too — samples too diverse to
        // have been reference-compressed. They constrain nothing; skip quietly.
        return Ok(None);
    }
    let f = File::open(path)?;
    let payload = zstd::stream::decode_all(std::io::BufReader::with_capacity(1 << 16, f))?;
    let mut r = &payload[..];
    let mut head = [0u8; 5];
    r.read_exact(&mut head)?;
    if &head[..4] != crate::refdelta::SKETCH_MAGIC {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "not a reference-delta sketch",
        ));
    }
    let mut fp = [0u8; 8];
    r.read_exact(&mut fp)?;
    let fingerprint = u64::from_le_bytes(fp);
    let _ref_db = read_string(&mut r)?;
    let _c = read_uvarint(&mut r)?;
    let _k = read_uvarint(&mut r)?;
    let _sample_file = read_string(&mut r)?;
    let mut tag = [0u8; 1];
    r.read_exact(&mut tag)?;
    if tag[0] != 0 {
        let _ = read_string(&mut r)?;
    }
    let mut skip1 = [0u8; 9]; // paired flag + mean read length
    r.read_exact(&mut skip1)?;
    let _num_reads = read_uvarint(&mut r)?;
    let hit_genomes = read_uvarint(&mut r)?;
    let _assigned = read_uvarint(&mut r)?;
    let _pool = read_uvarint(&mut r)?;
    let n_novel = read_uvarint(&mut r)?;
    let hit_len = read_uvarint(&mut r)? as usize;
    let pool_len = read_uvarint(&mut r)? as usize;
    let novel_len = read_uvarint(&mut r)? as usize;
    let _count_len = read_uvarint(&mut r)?;
    let _err_count = read_uvarint(&mut r)?;
    let _err_len = read_uvarint(&mut r)?;

    let mut skip = vec![0u8; hit_len + pool_len];
    r.read_exact(&mut skip)?;
    let mut novel_bytes = vec![0u8; novel_len];
    r.read_exact(&mut novel_bytes)?;
    let novel = read_hashes(&mut &novel_bytes[..])?;
    if novel.len() as u64 != n_novel {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "novel section decoded {} hashes but the header says {}",
                novel.len(),
                n_novel
            ),
        ));
    }

    let mut per_group: FxHashMap<u32, u32> = FxHashMap::default();
    for h in novel {
        if let Ok(i) = probe.binary_search_by_key(&h, |&(k, _)| k) {
            *per_group.entry(probe[i].1).or_insert(0) += 1;
        }
    }
    let mut out = String::new();
    if per_group.is_empty() {
        // Still emit the sample, so a run's coverage is visible.
        out.push_str(&format!(
            "{}\t{:016x}\t{}\t{}\t\t\n",
            path, fingerprint, hit_genomes, n_novel
        ));
        return Ok(Some(out));
    }
    let mut ids: Vec<u32> = per_group.keys().copied().collect();
    ids.sort_unstable();
    for gi in ids {
        out.push_str(&format!(
            "{}\t{:016x}\t{}\t{}\t{}\t{}\n",
            path, fingerprint, hit_genomes, n_novel, gi, per_group[&gi]
        ));
    }
    Ok(Some(out))
}

fn run_fingerprint(args: RefRecoverArgs) {
    let graph_path = args.graph.clone().unwrap_or_else(|| {
        error!("ref-recover fingerprint needs --graph from `ref-recover pairs`; exiting");
        std::process::exit(1);
    });
    let graph = load_graph(&graph_path);
    let order_path = args.genome_order.clone().unwrap_or_else(|| {
        error!("ref-recover fingerprint needs --genome-order; exiting");
        std::process::exit(1);
    });
    let rank = rank_from_order_file(&graph, &order_path);
    println!("{:016x}", fingerprint_for_order(&graph, &rank));
}

/// `rank[genome_id]` = that genome's position in the candidate build order.
fn rank_from_order_file(g: &Graph, path: &str) -> Vec<u32> {
    let text = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("Could not read --genome-order file {}: {}", path, e));
    let mut by_name: FxHashMap<&str, u32> = FxHashMap::default();
    for (i, gen) in g.genomes.iter().enumerate() {
        by_name.insert(gen.file_name.as_str(), i as u32);
    }
    // Tolerate a different path spelling, as `ref-build` does.
    let mut by_base: FxHashMap<&str, u32> = FxHashMap::default();
    for (i, gen) in g.genomes.iter().enumerate() {
        if let Some(b) = Path::new(gen.file_name.as_str())
            .file_name()
            .and_then(|s| s.to_str())
        {
            by_base.entry(b).or_insert(i as u32);
        }
    }
    let mut rank = vec![u32::MAX; g.genomes.len()];
    let mut pos = 0u32;
    for line in text.lines() {
        let name = line.trim();
        if name.is_empty() || name.starts_with('#') {
            continue;
        }
        let gid = by_name
            .get(name)
            .or_else(|| {
                Path::new(name)
                    .file_name()
                    .and_then(|b| b.to_str())
                    .and_then(|b| by_base.get(b))
            })
            .copied()
            .unwrap_or_else(|| panic!("--genome-order lists {}, which is not in the graph", name));
        rank[gid as usize] = pos;
        pos += 1;
    }
    if let Some(i) = rank.iter().position(|&r| r == u32::MAX) {
        panic!(
            "--genome-order does not list {} (genome id {})",
            g.genomes[i].file_name, i
        );
    }
    rank
}

// --- searching the residual --------------------------------------------------

/// A total build order consistent with "winner precedes the rest" for every
/// group, or `None` if those constraints contradict each other. Ties break by
/// genome id, so the same constraints always give the same order.
fn order_from_winners(ng: usize, groups: &[Group], winners: &[u32]) -> Option<Vec<u32>> {
    let mut succ: Vec<Vec<u32>> = vec![Vec::new(); ng];
    let mut indeg = vec![0u32; ng];
    for (gi, grp) in groups.iter().enumerate() {
        let w = winners[gi] as usize;
        for &m in &grp.members {
            if m != winners[gi] {
                succ[w].push(m);
                indeg[m as usize] += 1;
            }
        }
    }
    let mut heap: std::collections::BinaryHeap<std::cmp::Reverse<u32>> = (0..ng as u32)
        .filter(|&g| indeg[g as usize] == 0)
        .map(std::cmp::Reverse)
        .collect();
    let mut rank = vec![u32::MAX; ng];
    let mut pos = 0u32;
    while let Some(std::cmp::Reverse(g)) = heap.pop() {
        rank[g as usize] = pos;
        pos += 1;
        // A node is final once popped, so its successor list can be taken.
        for m in std::mem::take(&mut succ[g as usize]) {
            indeg[m as usize] -= 1;
            if indeg[m as usize] == 0 {
                heap.push(std::cmp::Reverse(m));
            }
        }
    }
    if (pos as usize) != ng {
        return None; // the evidence is cyclic under this assignment
    }
    Some(rank)
}

fn parse_solved(path: &str, n_groups: usize) -> Vec<Option<u32>> {
    let text = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("Could not read --solved file {}: {}", path, e));
    let mut out = vec![None; n_groups];
    for (lineno, line) in text.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut it = line.split('\t');
        let gi: usize = match it.next().and_then(|x| x.parse().ok()) {
            Some(v) => v,
            None => {
                error!(
                    "--solved line {}: expected <group_id><TAB><genome_id>",
                    lineno + 1
                );
                std::process::exit(1);
            }
        };
        let w: u32 = match it.next().and_then(|x| x.parse().ok()) {
            Some(v) => v,
            None => {
                error!(
                    "--solved line {}: expected <group_id><TAB><genome_id>",
                    lineno + 1
                );
                std::process::exit(1);
            }
        };
        if gi >= n_groups {
            error!(
                "--solved line {}: group {} is beyond the {} groups in the graph",
                lineno + 1,
                gi,
                n_groups
            );
            std::process::exit(1);
        }
        out[gi] = Some(w);
    }
    out
}

fn run_search(args: RefRecoverArgs) {
    let graph_path = args.graph.clone().unwrap_or_else(|| {
        error!("ref-recover search needs --graph from `ref-recover pairs`; exiting");
        std::process::exit(1);
    });
    let graph = load_graph(&graph_path);
    let target_hex = args.target.clone().unwrap_or_else(|| {
        error!("ref-recover search needs --target, the fingerprint the samples record; exiting");
        std::process::exit(1);
    });
    let target = u64::from_str_radix(target_hex.trim().trim_start_matches("0x"), 16)
        .unwrap_or_else(|_| {
            error!("--target {} is not 16 hex digits; exiting", target_hex);
            std::process::exit(1);
        });
    let out = args.output.clone().unwrap_or_else(|| {
        error!("ref-recover search needs -o/--output for the recovered order; exiting");
        std::process::exit(1);
    });

    let n_groups = graph.groups.len();
    let solved = match &args.solved {
        Some(p) => parse_solved(p, n_groups),
        None => vec![None; n_groups],
    };
    // Validate the supplied winners really are members, so a stale --solved file
    // fails loudly instead of quietly searching the wrong space.
    for (gi, w) in solved.iter().enumerate() {
        if let Some(w) = w {
            if !graph.groups[gi].members.contains(w) {
                error!(
                    "--solved says group {} is won by genome {}, which is not one of its members ({:?}); exiting",
                    gi, w, graph.groups[gi].members
                );
                std::process::exit(1);
            }
        }
    }

    let open: Vec<usize> = (0..n_groups).filter(|&i| solved[i].is_none()).collect();
    let combinations: u128 = open
        .iter()
        .map(|&i| graph.groups[i].members.len() as u128)
        .product::<u128>()
        .max(1);
    info!(
        "{} contested group(s): {} settled by evidence, {} to search ({} combination(s))",
        n_groups,
        n_groups - open.len(),
        open.len(),
        combinations
    );
    if combinations > args.max_candidates as u128 {
        error!(
            "{} combinations exceeds --max-candidates {}. Probe more samples to settle groups first: each solved group divides this number.",
            combinations, args.max_candidates
        );
        std::process::exit(1);
    }

    let ng = graph.genomes.len();
    let mut winners: Vec<u32> = (0..n_groups)
        .map(|i| solved[i].unwrap_or(graph.groups[i].members[0]))
        .collect();
    let mut digits = vec![0usize; open.len()];
    let mut tried: u64 = 0;
    let mut cyclic: u64 = 0;
    loop {
        for (d, &gi) in digits.iter().zip(open.iter()) {
            winners[gi] = graph.groups[gi].members[*d];
        }
        tried += 1;
        match order_from_winners(ng, &graph.groups, &winners) {
            Some(rank) => {
                if fingerprint_for_order(&graph, &rank) == target {
                    let mut by_rank: Vec<(u32, usize)> =
                        rank.iter().enumerate().map(|(g, &r)| (r, g)).collect();
                    by_rank.sort_unstable();
                    let f = File::create(&out)
                        .unwrap_or_else(|e| panic!("Could not create {}: {}", out, e));
                    let mut w = BufWriter::new(f);
                    for (_, g) in by_rank {
                        writeln!(w, "{}", graph.genomes[g].file_name).ok();
                    }
                    w.flush().ok();
                    info!(
                        "Matched {:016x} after {} candidate(s); wrote the order to {}",
                        target, tried, out
                    );
                    info!(
                        "Rebuild with: weebill ref-build <db> --genome-order {}",
                        out
                    );
                    return;
                }
            }
            None => cyclic += 1,
        }
        // odometer over the open groups
        let mut i = 0usize;
        loop {
            if i == digits.len() {
                error!(
                    "No assignment of the {} open group(s) reproduces {:016x} ({} tried, {} rejected as contradictory). The graph itself must differ from the lost reference: check --taxonomy, --pool-min-genomes, -c/-k and that the genome set is exactly the original one.",
                    open.len(), target, tried, cyclic
                );
                std::process::exit(1);
            }
            digits[i] += 1;
            if digits[i] < graph.groups[open[i]].members.len() {
                break;
            }
            digits[i] = 0;
            i += 1;
        }
    }
}
