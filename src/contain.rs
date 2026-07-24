use crate::cmdline::*;
use crate::constants::*;
use crate::inference::*;
use crate::merge::merge_sketches;
use crate::sketch::*;
use crate::twostage_db::{ScreenIndex, TwoStageDb};
use crate::types::*;
use fxhash::FxHashMap;
use log::*;
use rayon::prelude::*;
use statrs::distribution::{DiscreteCDF, Poisson};
use std::fs::File;
use std::io;
use std::io::prelude::*;
use std::io::BufReader;
use std::io::BufWriter;
use std::path::Path;
use std::sync::Mutex;

fn print_ani_result(ani_result: &AniResult, pseudotax: bool, writer: &mut Box<dyn Write + Send>) {
    let print_final_ani = format!("{:.2}", f64::min(ani_result.final_est_ani * 100., 100.));
    let lambda_print;
    if let AdjustStatus::Lambda(lambda) = ani_result.lambda {
        lambda_print = format!("{:.3}", lambda);
    } else if ani_result.lambda == AdjustStatus::High {
        lambda_print = "HIGH".to_string();
    } else {
        lambda_print = "LOW".to_string();
    }
    let low_ani = ani_result.ani_ci.0;
    let high_ani = ani_result.ani_ci.1;
    let low_lambda = ani_result.lambda_ci.0;
    let high_lambda = ani_result.lambda_ci.1;

    let ci_ani;
    if low_ani.is_none() || high_ani.is_none() {
        ci_ani = "NA-NA".to_string();
    } else {
        ci_ani = format!(
            "{:.2}-{:.2}",
            low_ani.unwrap() * 100.,
            high_ani.unwrap() * 100.
        );
    }

    let ci_lambda;
    if low_lambda.is_none() || high_lambda.is_none() {
        ci_lambda = "NA-NA".to_string();
    } else {
        ci_lambda = format!("{:.2}-{:.2}", low_lambda.unwrap(), high_lambda.unwrap());
    }

    //"Sample_file\tQuery_file\tTaxonomic_abundance\tSequence_abundance\tAdjusted_ANI\tEff_cov\tANI_5-95_percentile\tEff_lambda\tLambda_5-95_percentile\tMedian_cov\tMean_cov_geq1\tContainment_ind\tNaive_ANI\tContig_name",

    if !pseudotax {
        writeln!(
            writer,
            "{}\t{}\t{}\t{:.3}\t{}\t{}\t{}\t{:.0}\t{:.3}\t{}/{}\t{:.2}\t{}",
            ani_result.seq_name,
            ani_result.gn_name,
            print_final_ani,
            ani_result.final_est_cov,
            ci_ani,
            lambda_print,
            ci_lambda,
            ani_result.median_cov,
            ani_result.mean_cov,
            ani_result.containment_index.0,
            ani_result.containment_index.1,
            ani_result.naive_ani * 100.,
            ani_result.contig_name,
        )
        .expect("Error writing to file");
    } else {
        writeln!(
            writer,
            "{}\t{}\t{:.4}\t{:.4}\t{}\t{:.3}\t{}\t{}\t{}\t{:.0}\t{:.3}\t{}/{}\t{:.2}\t{}\t{}",
            ani_result.seq_name,
            ani_result.gn_name,
            ani_result.rel_abund.unwrap(),
            ani_result.seq_abund.unwrap(),
            print_final_ani,
            ani_result.final_est_cov,
            ci_ani,
            lambda_print,
            ci_lambda,
            ani_result.median_cov,
            ani_result.mean_cov,
            ani_result.containment_index.0,
            ani_result.containment_index.1,
            ani_result.naive_ani * 100.,
            ani_result.kmers_lost.unwrap(),
            ani_result.contig_name,
        )
        .expect("Error writing to file");
    }
}

fn get_chunks(indices: &Vec<usize>, steps: usize) -> Vec<Vec<usize>> {
    let mut start = 0;
    let mut end = steps;
    let len = indices.len();
    let mut return_chunks = vec![];

    while start < len {
        if end > len {
            end = len;
        }

        let chunk: Vec<usize> = (start..end).collect();
        start = end;
        end += steps;
        return_chunks.push(chunk);
    }
    return_chunks
}

/// FracMinHash keeps a hashed k-mer iff `hash < u64::MAX / c`. Because the
/// threshold shrinks as `c` grows, the set of k-mers kept at a large `c` (sparse)
/// is a strict subset of the set kept at a small `c` (dense). So we can turn an
/// already-loaded sketch into a *sparser* one for free by dropping every k-mer
/// whose hash is above the coarser threshold. We can NOT go the other way --
/// making a sketch denser requires re-reading the source sequence.
fn subsample_view(gs: &GenomeSketch, target_c: usize) -> GenomeSketch {
    let mut out = gs.clone();
    if target_c <= gs.c {
        // Already at (or denser than) the requested rate; nothing to drop.
        return out;
    }
    let thresh = u64::MAX / (target_c as u64);
    out.c = target_c;
    out.genome_kmers = gs
        .genome_kmers
        .iter()
        .copied()
        .filter(|h| *h < thresh)
        .collect();
    if let Some(p) = &gs.pseudotax_tracked_nonused_kmers {
        out.pseudotax_tracked_nonused_kmers =
            Some(p.iter().copied().filter(|h| *h < thresh).collect());
    }
    out
}

/// Path of the per-genome dense-sketch cache file for a given source fasta.
fn dense_cache_path(dir: &str, file_name: &str, dense_c: usize, k: usize) -> String {
    let h = fxhash::hash64(&file_name);
    format!(
        "{}/{:016x}.c{}.k{}{}",
        dir, h, dense_c, k, GENOME_SKETCH_SUFFIX
    )
}

/// Obtain a *dense* (`-c = dense_c`) sketch for a genome that passed the screen.
///   * database already as dense / denser -> reuse (subsampling if denser),
///   * database sparser -> (re)sketch the source fasta at `dense_c`,
///     reusing an on-disk cache and an in-memory cache so each fasta is only
///     ever sketched once. This is how a dense database is grown lazily, only
///     for the genomes that actually appear in samples.
fn densify_genome(
    db_sketch: &GenomeSketch,
    dense_c: usize,
    k: usize,
    min_spacing: usize,
    mem_cache: &Mutex<FxHashMap<String, GenomeSketch>>,
    disk_cache: &Option<String>,
) -> Option<GenomeSketch> {
    if db_sketch.c <= dense_c {
        return Some(subsample_view(db_sketch, dense_c));
    }

    let key = db_sketch.file_name.clone();
    if let Some(hit) = mem_cache.lock().unwrap().get(&key) {
        return Some(hit.clone());
    }

    if let Some(dir) = disk_cache {
        let path = dense_cache_path(dir, &key, dense_c, k);
        if Path::new(&path).exists() {
            if let Ok(file) = File::open(&path) {
                let reader = BufReader::with_capacity(10_000_000, file);
                if let Ok(gs) = bincode::deserialize_from::<_, GenomeSketch>(reader) {
                    mem_cache.lock().unwrap().insert(key, gs.clone());
                    return Some(gs);
                }
            }
        }
    }

    let sketched = sketch_genome(dense_c, k, &key, min_spacing, true)?;

    if let Some(dir) = disk_cache {
        let path = dense_cache_path(dir, &key, dense_c, k);
        match File::create(&path) {
            Ok(file) => {
                let mut writer = BufWriter::new(file);
                if let Err(e) = bincode::serialize_into(&mut writer, &sketched) {
                    warn!("Could not write dense-sketch cache {}: {}", path, e);
                }
            }
            Err(e) => warn!("Could not create dense-sketch cache {}: {}", path, e),
        }
    }

    mem_cache.lock().unwrap().insert(key, sketched.clone());
    Some(sketched)
}

/// Two-stage stage 1 + 2: screen `sequence_sketch` against the pooled
/// `screen_index` (a single inverted pass over the sample at the sparse
/// `screen_c`), then return *dense* (`dense_c`) sketches for only the genomes
/// that pass. `plain_genome_sketches` is `Some` only for the non-`.syl2db`
/// (densify) path; a `.syl2db` decodes survivors by index instead. The returned
/// set replaces the full database for the (expensive) dense profiling that follows.
fn compute_dense_survivors(
    args: &ContainArgs,
    screen_index: &ScreenIndex,
    plain_genome_sketches: Option<&Vec<GenomeSketch>>,
    sequence_sketch: &SequencesSketch,
    screen_c: usize,
    dense_c: usize,
    k: usize,
    min_spacing: usize,
    mem_cache: &Mutex<FxHashMap<String, GenomeSketch>>,
    disk_cache: &Option<String>,
    two_stage_db: Option<&TwoStageDb>,
) -> Vec<GenomeSketch> {
    // Stage 1: cheap, permissive screen (query-like settings, no CIs).
    let mut screen_args = args.clone();
    screen_args.pseudotax = false;
    screen_args.minimum_ani = Some(args.screen_ani);
    screen_args.no_ci = true;
    // The screen scores against sketches sub-sampled to `screen_c`, so the dense
    // `-M/--min-number-kmers` floor would over-reject here: a genome/contig has
    // only ~length/screen_c sparse k-mers, so the dense floor M corresponds to a
    // length of M*screen_c (e.g. 50*3000 = 150 kb), silently dropping smaller
    // genomes (viruses, plasmids, short contigs) that single-stage profiling
    // would report. Scale the floor to the screen resolution so a genome that
    // could clear the dense floor also clears the screen; genomes truly below the
    // floor are still rejected at the dense stage.
    screen_args.min_number_kmers = args.min_number_kmers * dense_c as f64 / screen_c as f64;

    // The densify fallback (no .syl2db) re-sketches whole FASTA files keyed by
    // file name, which collapses `--individual-records` databases (many records
    // share a file name) into a single whole-file sketch. Reject that up front;
    // `db-convert` to a .syl2db preserves individual records and works fine.
    if let Some(genome_sketches) = plain_genome_sketches {
        let mut seen = std::collections::HashSet::with_capacity(genome_sketches.len());
        for gs in genome_sketches.iter() {
            if !seen.insert(gs.file_name.as_str()) {
                log::error!(
                    "--two-stage on a raw --individual-records database is unsupported \
                     (densification re-sketches whole files). Convert it with `sylph db-convert` \
                     first (a .syl2db preserves individual records). Exiting."
                );
                std::process::exit(1);
            }
        }
    }

    // Stage 1 = "Path B": one inverted pass over the sample produces, per genome,
    // the same matched-coverage multiset the per-genome `get_stats` loop would
    // collect; feeding it to the same `finalize_stats` + `--screen-min-matches` +
    // `min_number_kmers` checks reproduces the survivor set exactly, in O(sample)
    // rather than O(reference) work. (See experiments/7_mphf_screen_again.)
    let name_of = |g: usize| -> String {
        match two_stage_db {
            Some(db) => db.genome_file_name(g as u32).to_string(),
            None => plain_genome_sketches.unwrap()[g].file_name.clone(),
        }
    };
    let hits: Vec<(u32, Vec<u32>)> = screen_index
        .gather_hits(sequence_sketch)
        .into_iter()
        .collect();
    // Optional per-survivor dump: (genome, matched k-mers, total screen k-mers,
    // naive ANI, adjusted ANI, median coverage).
    let dump: Mutex<Vec<(String, usize, usize, f64, f64, f64)>> = Mutex::new(vec![]);
    let mut survivors: Vec<usize> = hits
        .into_par_iter()
        .filter_map(|(g, covs)| {
            let g = g as usize;
            let n_kmers = screen_index.sparse_count[g] as usize;
            // Mirror get_stats: reject genomes below the (scaled) k-mer floor.
            if (n_kmers as f64) < screen_args.min_number_kmers {
                return None;
            }
            let contain_count = covs.len();
            // finalize_stats applies the screen min-ANI gate (returns None below it).
            let fin = finalize_stats(&screen_args, k, n_kmers, contain_count, covs, None)?;
            // Require at least `--screen-min-matches` matched screen k-mers; this
            // cheaply drops genomes that clear the (permissive) screen ANI on only
            // a few chance-shared k-mers, before they cost a dense decode. Default
            // 1 == current behaviour (get_stats already needs >= 1 match).
            if contain_count < args.screen_min_matches {
                return None;
            }
            if args.screen_dump.is_some() {
                dump.lock().unwrap().push((
                    name_of(g),
                    contain_count,
                    n_kmers,
                    fin.naive_ani * 100.,
                    fin.final_est_ani * 100.,
                    fin.median_cov,
                ));
            }
            Some(g)
        })
        .collect();
    survivors.sort_unstable();
    log::info!(
        "{}: stage-1 screen (c={}, min-ANI {}) kept {} / {} candidate genomes",
        sequence_sketch.file_name,
        screen_c,
        args.screen_ani,
        survivors.len(),
        screen_index.num_genomes()
    );
    if let Some(path) = &args.screen_dump {
        let mut f =
            BufWriter::new(File::create(path).expect("could not create --screen-dump file"));
        writeln!(f, "Genome_file\tscreen_matched_kmers\tscreen_total_kmers\tnaive_ani\tscreen_adjusted_ani\tscreen_median_cov").unwrap();
        for (g, m, t, na, ea, mc) in dump.into_inner().unwrap() {
            writeln!(f, "{}\t{}\t{}\t{:.3}\t{:.3}\t{}", g, m, t, na, ea, mc).unwrap();
        }
        log::info!("Wrote stage-1 screen dump to {}", path);
    }

    // Stage 2. With a .syl2db we decode each survivor's dense block into a
    // short-lived sketch, run the (prefetch-friendly) pass-1 profiling on it, and
    // keep it only if it passes -- so the discarded majority is freed immediately
    // and never cached, keeping peak RAM proportional to the genomes that survive
    // rather than to the (much larger) screen-survivor set. With a plain database
    // the dense sketch is derived/re-sketched by densify_genome and filtered the
    // same way.
    let dense: Mutex<Vec<GenomeSketch>> = Mutex::new(vec![]);
    survivors.par_iter().for_each(|i| {
        let g = match two_stage_db {
            Some(db) => match db.decode_dense(*i as u32) {
                Ok(g) => Some(g),
                Err(e) => {
                    warn!("Could not decode dense block for genome index {}: {}", i, e);
                    None
                }
            },
            None => densify_genome(
                &plain_genome_sketches.unwrap()[*i],
                dense_c,
                k,
                min_spacing,
                mem_cache,
                disk_cache,
            ),
        };
        if let Some(g) = g {
            // Keep only genomes that pass pass-1 profiling; drop (free) the rest.
            if get_stats(args, &g, sequence_sketch, None, false).is_some() {
                dense.lock().unwrap().push(g);
            }
        }
    });
    let dense = dense.into_inner().unwrap();
    log::info!(
        "{}: stage-2 dense profiling (c={}) against {} genomes",
        sequence_sketch.file_name,
        dense_c,
        dense.len()
    );
    dense
}

pub fn contain(mut args: ContainArgs, pseudotax_in: bool) {
    if pseudotax_in {
        args.pseudotax = true;
    }

    let level;
    if args.trace {
        level = log::LevelFilter::Trace;
    } else if args.debug {
        level = log::LevelFilter::Debug;
    } else {
        level = log::LevelFilter::Info;
    }

    simple_logger::SimpleLogger::new()
        .with_level(level)
        .init()
        .unwrap();

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap();

    if args.estimate_read_counts {
        args.estimate_unknown = true;
        log::info!("--estimate-read-counts detected, also enabling -u. Sequence_abundance column will be set to estimated read counts, not abundance. This is still experimental.");
    }

    // Validate --apply-unknown BEFORE the output writer is created: `File::create`
    // below truncates the `-o` path, so any invalid --apply-unknown invocation that
    // errors out afterwards would first destroy an existing output file (and, when
    // `-o` is the input TSV itself, the input). All apply-unknown gating therefore
    // lives here, ahead of output creation.
    if let Some(tsv) = args.apply_unknown.clone() {
        if !pseudotax_in {
            log::error!("--apply-unknown is only supported for `profile` (it rewrites a profile TSV). Exiting.");
            std::process::exit(1);
        }
        if args.merge {
            log::error!("--apply-unknown cannot be combined with --merge: it rescales an existing per-sample profile TSV, so give the same per-sample inputs used to produce it. Exiting.");
            std::process::exit(1);
        }
        if args.estimate_read_counts {
            log::error!("--apply-unknown cannot be combined with --estimate-read-counts (the read-count formula differs). Re-run the original profile with --estimate-read-counts instead. Exiting.");
            std::process::exit(1);
        }
        if args.estimate_unknown {
            log::error!(
                "--apply-unknown already produces the -u profile; do not also pass -u. Exiting."
            );
            std::process::exit(1);
        }
        // `-o` pointing at the input TSV would wipe the input before it is read.
        // Compare canonical paths when both resolve; else a literal comparison.
        if let Some(out) = &args.out_file_name {
            let same = match (std::fs::canonicalize(&tsv), std::fs::canonicalize(out)) {
                (Ok(a), Ok(b)) => a == b,
                _ => Path::new(&tsv) == Path::new(out),
            };
            if same {
                log::error!(
                    "--apply-unknown: the output file (-o `{}`) is the same as the input profile TSV; this would overwrite the input before it is read. Write to a different path. Exiting.",
                    out
                );
                std::process::exit(1);
            }
        }
    }

    let out_writer = match args.out_file_name {
        Some(ref x) => {
            let path = Path::new(&x);
            Box::new(BufWriter::new(File::create(path).unwrap())) as Box<dyn Write + Send>
        }
        None => Box::new(BufWriter::new(io::stdout())) as Box<dyn Write + Send>,
    };

    log::info!("Obtaining sketches...");
    let mut genome_sketch_files = vec![];
    let mut genome_files = vec![];
    let mut two_stage_db_files = vec![];
    let mut read_sketch_files = vec![];
    let mut read_files = vec![];

    let mut all_files = args.files.clone();

    if let Some(ref newline_file) = args.file_list {
        let file = File::open(newline_file).unwrap();
        let reader = BufReader::new(file);
        for line in reader.lines() {
            all_files.push(line.unwrap());
        }
    }

    for file in all_files.iter() {
        let mut genome_sketch_good_suffix = false;
        for suff in QUERY_FILE_SUFFIX_VALID {
            if file.ends_with(suff) {
                genome_sketch_good_suffix = true;
                break;
            }
        }

        let mut sample_sketch_good_suffix = false;
        for suff in SAMPLE_FILE_SUFFIX_VALID {
            if file.ends_with(suff) {
                sample_sketch_good_suffix = true;
                break;
            }
        }
        // reference-delta compressed samples are decoded via --reference in get_seq_sketch
        if file.ends_with(REF_SAMPLE_SUFFIX) {
            sample_sketch_good_suffix = true;
        }

        if file.ends_with(TWO_STAGE_DB_SUFFIX) {
            two_stage_db_files.push(file);
        } else if genome_sketch_good_suffix {
            genome_sketch_files.push(file);
        } else if sample_sketch_good_suffix {
            read_sketch_files.push(file);
        } else if is_fasta(file) {
            genome_files.push(file);
        } else if is_fastq(file) {
            read_files.push(vec![file]);
        } else {
            warn!(
                "{} file extension is not a sketch or a fasta/fastq file.",
                &file
            );
        }
    }

    if args.first_pair.len() != args.second_pair.len() {
        error!("Different number of paired sequences (-1, -2) for sketching. Exiting.");
        std::process::exit(1);
    }

    // zip together the first and second pair files, push them to read_files
    for (first, second) in args.first_pair.iter().zip(args.second_pair.iter()) {
        read_files.push(vec![first, second]);
    }

    for read in args.reads.iter() {
        read_files.push(vec![read]);
    }

    // Interleaved inputs are encoded as a 3-element group (the path repeated) so
    // get_seq_sketch can distinguish them by length from single (1) and paired (2)
    // inputs; only the first element is actually used when sketching.
    for read in args.interleaved.iter() {
        read_files.push(vec![read, read, read]);
    }

    if genome_sketch_files.is_empty() && genome_files.is_empty() && two_stage_db_files.is_empty() {
        log::error!("No genome files found; see sylph query/profile -h for help. Exiting");
        std::process::exit(1);
    }

    if read_sketch_files.is_empty() && read_files.is_empty() {
        log::error!("No read files found; see sylph query/profile -h for help. Exiting");
        std::process::exit(1);
    }

    // Load the reference DB (if given) once, for decoding *.sylspr samples.
    let have_refdelta = read_sketch_files
        .iter()
        .any(|f| f.ends_with(REF_SAMPLE_SUFFIX));
    let ref_db: Option<crate::refdelta::RefIndex> = match &args.reference {
        Some(path) => {
            log::info!(
                "Loading reference database {} for .sylspr decoding...",
                path
            );
            let file = File::open(path)
                .unwrap_or_else(|_| panic!("Could not open reference database {}", path));
            Some(
                crate::refdelta::open_ref_index_file(file).unwrap_or_else(|e| {
                    panic!("{} is not a valid reference database: {}", path, e)
                }),
            )
        }
        None => {
            if have_refdelta {
                log::error!("Reference-delta compressed samples (*.sylspr) were given but --reference was not specified. Exiting.");
                std::process::exit(1);
            }
            None
        }
    };

    // A .syl2db is a self-contained two-stage database: open it (loading only the
    // sparse stage-1 index) and use its per-genome sparse sketches as the screen.
    let two_stage_db: Option<TwoStageDb> = if !two_stage_db_files.is_empty() {
        if two_stage_db_files.len() > 1
            || !genome_sketch_files.is_empty()
            || !genome_files.is_empty()
        {
            log::error!(
                "A two-stage database ({}) must be the only genome input. Exiting",
                TWO_STAGE_DB_SUFFIX
            );
            std::process::exit(1);
        }
        if !args.pseudotax {
            log::error!("Two-stage databases ({}) are only supported for `sylph profile`, not `sylph query`. Exiting", TWO_STAGE_DB_SUFFIX);
            std::process::exit(1);
        }
        log::info!(
            "Opening two-stage database {} (loading stage-1 sparse index)...",
            two_stage_db_files[0]
        );
        let db = crate::twostage_db::open_file(two_stage_db_files[0]).unwrap_or_else(|e| {
            panic!(
                "{} is not a valid two-stage database: {}",
                two_stage_db_files[0], e
            )
        });
        // The database is inherently two-stage; enable the screen-then-densify path.
        args.two_stage = true;
        Some(db)
    } else {
        None
    };

    // A .syl2db screens via its pooled stage-1 index and decodes dense blocks by
    // index for stage 2, so it carries no in-memory `GenomeSketch` list. Only the
    // plain-database path materializes one here.
    let genome_sketches = match &two_stage_db {
        Some(_) => Vec::new(),
        None => get_genome_sketches(&args, &genome_sketch_files, &genome_files),
    };
    log::info!("Finished obtaining genome sketches.");

    match &two_stage_db {
        Some(db) => {
            if db.is_empty() {
                log::error!("Two-stage database contains no genomes. Exiting");
                std::process::exit(1);
            }
        }
        None => {
            if genome_sketches.is_empty() {
                log::error!(
                    "No genome sketches found; see sylph query/profile -h for help. Exiting"
                );
                std::process::exit(1);
            }
            if genome_sketches
                .first()
                .unwrap()
                .pseudotax_tracked_nonused_kmers
                .is_none()
                && args.pseudotax
            {
                log::error!("Attempting profiling, but *.syldb was sketched with the --disable-profiling option. Exiting");
                std::process::exit(1);
            }
        }
    }

    // ---- Two-stage profiling setup ----------------------------------------
    // A .syl2db dictates its own dense rate (`c`) and stage-1 screen rate
    // (`screen_c`); otherwise they come from the database `-c` and the flags.
    let (db_c, db_k, screen_c, dense_c) = match &two_stage_db {
        Some(db) => (db.c, db.k, db.screen_c, db.c),
        None => {
            let db_c = genome_sketches[0].c;
            (
                db_c,
                genome_sketches[0].k,
                args.screen_c.unwrap_or(SCREEN_C_DEFAULT).max(db_c),
                args.dense_c,
            )
        }
    };
    let dense_mem_cache: Mutex<FxHashMap<String, GenomeSketch>> = Mutex::new(FxHashMap::default());
    if args.two_stage {
        if !args.pseudotax {
            log::error!(
                "--two-stage is only supported for `sylph profile`, not `sylph query`. Exiting"
            );
            std::process::exit(1);
        }
        if dense_c > screen_c {
            log::error!("--dense-c ({}) must be <= the screen -c ({}); the dense stage cannot be sparser than the screen. Exiting", dense_c, screen_c);
            std::process::exit(1);
        }
        if let Some(dir) = &args.dense_cache {
            if let Err(e) = std::fs::create_dir_all(dir) {
                log::error!(
                    "Could not create --dense-cache directory {}: {}. Exiting",
                    dir,
                    e
                );
                std::process::exit(1);
            }
        }
        if two_stage_db.is_some() {
            log::info!(
                "Two-stage database profiling: screen at c={} (min-ANI {}), dense profile at c={} by decoding only the screened genomes' compressed blocks.",
                screen_c, args.screen_ani, dense_c
            );
        } else {
            log::info!(
                "Two-stage profiling enabled: screen at c={} (min-ANI {}), dense profile at c={}{}.",
                screen_c, args.screen_ani, dense_c,
                if db_c > dense_c { " by (re)sketching candidate genomes from their source fasta" } else { "" }
            );
        }
    }
    // The sample sketch must not be sparser than the genome rate it is compared
    // against: the full DB rate normally, or the dense rate under two-stage.
    let effective_genome_c = if args.two_stage { dense_c } else { db_c };

    // Stage-1 screen index. A .syl2db carries its own (loaded from file); the
    // plain-database two-stage path builds one in memory once, from each genome
    // subsampled to `screen_c`.
    let inmem_screen_index: Option<ScreenIndex> = if args.two_stage && two_stage_db.is_none() {
        log::info!("Building in-memory stage-1 screen index (c={}).", screen_c);
        let sparse_per_genome: Vec<Vec<u64>> = genome_sketches
            .par_iter()
            .map(|gs| subsample_view(gs, screen_c).genome_kmers)
            .collect();
        Some(ScreenIndex::build(&sparse_per_genome, screen_c, db_k))
    } else {
        None
    };

    let num_raw_read_files = read_files.len();
    let step;
    if let Some(sample_threads) = args.sample_threads {
        if sample_threads > 0 {
            step = sample_threads;
        } else {
            step = 1;
        }
    } else {
        if args.pseudotax {
            step = usize::max(
                args.threads / 3 + 1,
                usize::min(num_raw_read_files, args.threads),
            )
        } else {
            step = usize::max(1, usize::min(num_raw_read_files, args.threads))
        }
    }

    let read_sketch_files_as_vec = read_sketch_files
        .clone()
        .into_iter()
        .map(|x| vec![x])
        .collect::<Vec<Vec<&String>>>();
    read_files.extend(read_sketch_files_as_vec);
    let sequence_index_vec = (0..read_files.len()).collect::<Vec<usize>>();
    let out_writer: Mutex<Box<dyn Write + Send>> = Mutex::new(out_writer);

    // --apply-unknown: rescale an existing (non-`-u`) profile TSV into the profile
    // `-u` would have produced, without re-profiling. Everything `-u` changes is a
    // per-sample scalar applied to two printed columns (see `estimate_true_cov` /
    // `estimate_covered_bases`), so we only need the sample sketches (k-mer identity
    // + read length) and the database (per-genome sizes) -- both required here.
    if let Some(tsv_path) = args.apply_unknown.clone() {
        apply_unknown_from_tsv(
            &args,
            &tsv_path,
            &two_stage_db,
            &genome_sketches,
            &read_files,
            read_sketch_files.len(),
            ref_db.as_ref(),
            effective_genome_c,
            db_k,
            &out_writer,
        );
        log::info!("sylph finished.");
        return;
    }

    let chunks = get_chunks(&sequence_index_vec, step);

    print_header(
        args.pseudotax,
        &mut out_writer.lock().unwrap(),
        args.estimate_unknown,
    );
    // The per-sample profiling body: screen (if two-stage), compute per-genome
    // stats, reassign k-mers for taxonomic profiling, then print. Shared by the
    // normal per-input path and the single merged sample produced by --merge.
    let process_sample = |sequence_sketch: SequencesSketch, first_read_file: &str| {
        {
            let kmer_id_opt;
            if args.seq_id.is_some() {
                kmer_id_opt = Some((args.seq_id.unwrap() / 100.).powf(sequence_sketch.k as f64));
            } else {
                kmer_id_opt = get_kmer_identity(&sequence_sketch, args.estimate_unknown);
                log::debug!(
                    "{} has estimated identity {:.3}.",
                    &first_read_file,
                    kmer_id_opt.unwrap().powf(1. / sequence_sketch.k as f64) * 100.
                );
            }

            // Under two-stage, screen first and replace the full database
            // with dense sketches of only the genomes that pass.
            let dense_local: Vec<GenomeSketch>;
            let active_sketches: &Vec<GenomeSketch> = if args.two_stage {
                let screen_index = match &two_stage_db {
                    Some(db) => &db.screen_index,
                    None => inmem_screen_index.as_ref().unwrap(),
                };
                let plain_genome_sketches = two_stage_db.is_none().then_some(&genome_sketches);
                dense_local = compute_dense_survivors(
                    &args,
                    screen_index,
                    plain_genome_sketches,
                    &sequence_sketch,
                    screen_c,
                    dense_c,
                    db_k,
                    args.min_spacing_kmer,
                    &dense_mem_cache,
                    &args.dense_cache,
                    two_stage_db.as_ref(),
                );
                &dense_local
            } else {
                &genome_sketches
            };
            let active_index_vec = (0..active_sketches.len()).collect::<Vec<usize>>();

            let stats_vec_seq: Mutex<Vec<AniResult>> = Mutex::new(vec![]);
            active_index_vec.par_iter().for_each(|i| {
                let genome_sketch = &active_sketches[*i];
                let res = get_stats(
                    &args,
                    genome_sketch,
                    &sequence_sketch,
                    None,
                    args.log_reassignments,
                );
                if res.is_some() {
                    //res.as_mut().unwrap().genome_sketch_index = *i;
                    stats_vec_seq.lock().unwrap().push(res.unwrap());
                }
            });

            let mut stats_vec_seq = stats_vec_seq.into_inner().unwrap();
            estimate_true_cov(
                &mut stats_vec_seq,
                kmer_id_opt,
                args.estimate_unknown,
                sequence_sketch.mean_read_length,
                sequence_sketch.k,
            );

            if args.pseudotax {
                log::info!(
                    "{} taxonomic profiling; reassigning k-mers for {} genomes...",
                    &first_read_file,
                    stats_vec_seq.len()
                );
                let winner_map = winner_table(&stats_vec_seq, args.log_reassignments);
                let remaining_genomes = stats_vec_seq
                    .iter()
                    .map(|x| x.genome_sketch)
                    .collect::<Vec<&GenomeSketch>>();
                let stats_vec_seq_2 = Mutex::new(vec![]);
                remaining_genomes.into_par_iter().for_each(|genome_sketch| {
                    let res = get_stats(
                        &args,
                        genome_sketch,
                        &sequence_sketch,
                        Some(&winner_map),
                        args.log_reassignments,
                    );
                    if res.is_some() {
                        stats_vec_seq_2.lock().unwrap().push(res.unwrap());
                    }
                });
                stats_vec_seq = derep_if_reassign_threshold(
                    &stats_vec_seq,
                    stats_vec_seq_2.into_inner().unwrap(),
                    args.redundant_ani,
                    sequence_sketch.k,
                );
                //stats_vec_seq = stats_vec_seq_2.into_inner().unwrap();
                estimate_true_cov(
                    &mut stats_vec_seq,
                    kmer_id_opt,
                    args.estimate_unknown,
                    sequence_sketch.mean_read_length,
                    sequence_sketch.k,
                );
                log::info!(
                    "{} has {} genomes passing profiling threshold. ",
                    &first_read_file,
                    stats_vec_seq.len()
                );

                let mut bases_explained = 1.;
                if args.estimate_unknown {
                    bases_explained = estimate_covered_bases(
                        &stats_vec_seq,
                        &sequence_sketch,
                        sequence_sketch.mean_read_length,
                        sequence_sketch.k,
                    );
                    log::info!(
                        "{} has {:.2}% of reads detected in database by profile",
                        &first_read_file,
                        bases_explained * 100.
                    );
                }

                let total_cov = stats_vec_seq.iter().map(|x| x.final_est_cov).sum::<f64>();
                let total_seq_cov = stats_vec_seq
                    .iter()
                    .map(|x| x.final_est_cov * x.genome_sketch.gn_size as f64)
                    .sum::<f64>();
                for thing in stats_vec_seq.iter_mut() {
                    thing.rel_abund = Some(thing.final_est_cov / total_cov * 100.);
                }
                for thing in stats_vec_seq.iter_mut() {
                    if args.estimate_read_counts {
                        thing.seq_abund = Some(
                            (thing.final_est_cov * thing.genome_sketch.gn_size as f64
                                / sequence_sketch.mean_read_length
                                * bases_explained)
                                .round(),
                        );
                    } else {
                        let seq_abund = thing.final_est_cov * thing.genome_sketch.gn_size as f64
                            / total_seq_cov
                            * 100.
                            * bases_explained;
                        thing.seq_abund = Some(seq_abund);
                    }
                }
            }

            if args.pseudotax {
                stats_vec_seq.sort_by(|x, y| {
                    y.rel_abund
                        .unwrap()
                        .partial_cmp(&x.rel_abund.unwrap())
                        .unwrap()
                });
            } else {
                stats_vec_seq
                    .sort_by(|x, y| y.final_est_ani.partial_cmp(&x.final_est_ani).unwrap());
            }

            let mut out_writer = out_writer.lock().unwrap();
            for res in stats_vec_seq {
                print_ani_result(&res, args.pseudotax, &mut out_writer);
            }
        }
    };

    if args.merge {
        // Collapse every read input (raw single/paired/interleaved plus pre-sketched
        // *.sylsp/*.sylspc/*.sylspr samples) into one sketch and profile it once.
        let n_sketch_files = read_sketch_files.len();
        let total = read_files.len();
        let n_raw = total - n_sketch_files;
        // Sketch/load every input concurrently. A sequential loop would sketch each raw group
        // to EOF before opening the next, so raw inputs that are FIFOs fed by one upstream
        // writer (split R1/R2/orphan pipes) could fill their kernel buffers and deadlock the
        // producer -- the same streaming hazard `sketch` guards against. Give every raw input
        // its own thread so all streams stay drained at once; merge order does not change the
        // summed result. Pre-sketched inputs are plain file reads and never block.
        let merge_pool = rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads.max(n_raw).max(1))
            .build()
            .expect("Failed to build merge read-sketching thread pool");
        let collected_opt: Vec<Option<(SequencesSketch, ReadSketchMeta)>> =
            merge_pool.install(|| {
                (0..total)
                    .into_par_iter()
                    .map(|j| {
                        let is_sketch = j >= total - n_sketch_files;
                        get_seq_sketch_with_meta(
                            &args,
                            &read_files[j],
                            is_sketch,
                            ref_db.as_ref(),
                            effective_genome_c,
                            db_k,
                        )
                    })
                    .collect()
            });
        // A `None` means an input could not be sketched or loaded (e.g. incompatible -c/-k,
        // or a pre-sketched sample whose c exceeds the database's). Dropping it would profile
        // a merged sample that silently excludes those reads and misreports
        // containment/coverage, so refuse the whole merge instead of merging a subset.
        let failed: Vec<&str> = collected_opt
            .iter()
            .zip(read_files.iter())
            .filter(|(opt, _)| opt.is_none())
            .map(|(_, rf)| rf[0].as_str())
            .collect();
        if !failed.is_empty() {
            log::error!(
                "--merge: {} read input(s) could not be sketched or loaded ({}); refusing to profile a partial merged sample. Exiting.",
                failed.len(),
                failed.join(", ")
            );
            std::process::exit(1);
        }
        let mut collected: Vec<(SequencesSketch, ReadSketchMeta)> =
            collected_opt.into_iter().flatten().collect();
        if collected.is_empty() {
            log::error!("--merge: no read samples were provided. Exiting.");
            std::process::exit(1);
        }
        let merged_name = args
            .sample_name
            .clone()
            .unwrap_or_else(|| "merged".to_string());
        let n = collected.len();
        let (mut merged, _meta) = if n == 1 {
            collected.pop().unwrap()
        } else {
            merge_sketches(&collected, Some(merged_name.clone()))
        };
        merged.sample_name = Some(merged_name.clone());
        log::info!(
            "Profiling merged sample '{}' built from {} read input(s)...",
            merged_name,
            n
        );
        process_sample(merged, &merged_name);
        log::info!("Finished merged sample {}.", merged_name);
    } else {
        chunks.into_iter().for_each(|chunk| {
            chunk.into_par_iter().for_each(|j| {
                let is_sketch = j >= read_files.len() - read_sketch_files.len();
                let sequence_sketch = get_seq_sketch(
                    &args,
                    &read_files[j],
                    is_sketch,
                    ref_db.as_ref(),
                    effective_genome_c,
                    db_k,
                );
                if let Some(sequence_sketch) = sequence_sketch {
                    process_sample(sequence_sketch, read_files[j][0]);
                }
                if read_files[j].len() > 1 {
                    log::info!("Finished paired sample {}.", &read_files[j][0]);
                } else {
                    log::info!("Finished sample {}.", &read_files[j][0]);
                }
            });
        });
    }

    log::info!("sylph finished.");
}

fn derep_if_reassign_threshold<'a>(
    results_old: &Vec<AniResult>,
    results_new: Vec<AniResult<'a>>,
    ani_thresh: f64,
    k: usize,
) -> Vec<AniResult<'a>> {
    let ani_thresh = ani_thresh / 100.;

    let mut gn_sketch_to_contain = FxHashMap::default();
    for result in results_old.iter() {
        gn_sketch_to_contain.insert(result.genome_sketch, result);
    }

    let threshold = f64::powf(ani_thresh, k as f64);
    let mut return_vec = vec![];
    for result in results_new.into_iter() {
        let old_res = &gn_sketch_to_contain[result.genome_sketch];
        let num_kmer_reassign = (old_res.containment_index.0 - result.containment_index.0) as f64;
        let reass_thresh = threshold * result.containment_index.1 as f64;
        if num_kmer_reassign < reass_thresh {
            return_vec.push(result);
        } else {
            log::debug!(
                "genome {} had num k-mers reassigned = {}, threshold was {}, removing.",
                result.gn_name,
                num_kmer_reassign,
                reass_thresh
            );
        }
    }
    return_vec
}

fn estimate_true_cov(
    results: &mut Vec<AniResult>,
    kmer_id_opt: Option<f64>,
    estimate_unknown: bool,
    read_length: f64,
    k: usize,
) {
    let mut multiplier = 1.;
    if estimate_unknown {
        multiplier = read_length / (read_length - k as f64 + 1.);
    }
    if estimate_unknown && kmer_id_opt.is_some() {
        let id = kmer_id_opt.unwrap();
        for res in results.iter_mut() {
            res.final_est_cov = res.final_est_cov / id * multiplier;
        }
    }
}

/// Per-sample scalars needed to turn a non-`-u` profile into the `-u` profile.
struct UnknownScalars {
    /// `multiplier / id` -- the factor `estimate_true_cov` applies uniformly to
    /// every genome's `final_est_cov` (Eff_cov -> True_cov).
    cov_scale: f64,
    /// `c * (sum of sample k-mer counts) * multiplier` -- the denominator of
    /// `estimate_covered_bases` (independent of which genomes were detected).
    tentative_bases: f64,
}

/// Recompute the per-sample `-u` scalars from a sample sketch, mirroring exactly
/// what `process_sample` does when `estimate_unknown` is set: the k-mer identity
/// (`-I` override or `get_kmer_identity`), the read-length multiplier, and the
/// `estimate_covered_bases` denominator.
fn unknown_scalars_for_sample(args: &ContainArgs, sketch: &SequencesSketch) -> UnknownScalars {
    let id = if let Some(seq_id) = args.seq_id {
        (seq_id / 100.).powf(sketch.k as f64)
    } else {
        // With estimate_unknown effectively on, get_kmer_identity always returns Some.
        get_kmer_identity(sketch, true).unwrap()
    };
    let read_length = sketch.mean_read_length;
    let multiplier = read_length / (read_length - sketch.k as f64 + 1.);
    let total_counts: usize = sketch.kmer_counts.values().map(|c| *c as usize).sum();
    UnknownScalars {
        cov_scale: multiplier / id,
        tentative_bases: sketch.c as f64 * total_counts as f64 * multiplier,
    }
}

/// `profile --apply-unknown`: rewrite an existing (non-`-u`) profile TSV into the
/// profile `-u` would have produced, without re-profiling.
///
/// Everything `-u` changes is a per-sample scalar applied to two printed columns
/// (see [`estimate_true_cov`]/[`estimate_covered_bases`]): `Eff_cov` is scaled by
/// `multiplier/id` (and relabelled `True_cov`), and `Sequence_abundance` is scaled
/// by the estimated fraction of covered bases. `Taxonomic_abundance` and the ANIs
/// are unchanged. We therefore only need the sample sketches (for `id`, read
/// length and total k-mer counts) and the database (for each detected genome's
/// size, `gn_size`) -- both are required on the command line.
///
/// Because the only per-genome input is the already-rounded `Eff_cov` printed in
/// the TSV (not the full-precision internal value a real `-u` run rescales), the
/// output is not guaranteed to be bit-for-bit identical to a real `-u` run:
/// individual `True_cov`/`Sequence_abundance` cells may differ by a unit in the
/// last printed place (the `Sequence_abundance` scalar can drift a little more, as
/// it sums the rounded per-row coverage across all rows).
#[allow(clippy::too_many_arguments)]
fn apply_unknown_from_tsv(
    args: &ContainArgs,
    tsv_path: &str,
    two_stage_db: &Option<TwoStageDb>,
    genome_sketches: &[GenomeSketch],
    read_files: &[Vec<&String>],
    n_sketch_files: usize,
    ref_db: Option<&crate::refdelta::RefIndex>,
    effective_genome_c: usize,
    db_k: usize,
    out_writer: &Mutex<Box<dyn Write + Send>>,
) {
    // Genome_file column -> genome size, from the database (never the TSV).
    // A profile TSV identifies each detected genome only by its `Genome_file`
    // (source fasta path). With `--individual-records` several genome sketches
    // share one fasta and are disambiguated only by `Contig_name`, so keying by
    // `Genome_file` alone would collide and yield wrong sizes for all but one
    // record. Rather than silently mis-scale, refuse such a database: build the
    // map and croak on the first duplicate name.
    let db_genome_sizes: Vec<(String, usize)> = match two_stage_db {
        Some(db) => db
            .genome_sizes()
            .into_iter()
            .map(|(name, size)| (name.to_string(), size))
            .collect(),
        None => genome_sketches
            .iter()
            .map(|gs| (gs.file_name.clone(), gs.gn_size))
            .collect(),
    };
    let mut gn_size_map: FxHashMap<String, usize> = FxHashMap::default();
    for (name, size) in db_genome_sizes {
        if gn_size_map.insert(name.clone(), size).is_some() {
            log::error!(
                "--apply-unknown: the database contains multiple genomes with the same Genome_file `{}` (an --individual-records database). A profile TSV identifies genomes by Genome_file alone, so their sizes cannot be told apart; --apply-unknown does not support such databases. Re-run the original profile with -u instead. Exiting.",
                name
            );
            std::process::exit(1);
        }
    }

    // Per-sample scalars, keyed by the Sample_file value each sketch prints
    // (sample_name if set, else file_name -- matching `AniResult::seq_name`).
    let total = read_files.len();
    let mut scalars: FxHashMap<String, UnknownScalars> = FxHashMap::default();
    for (j, read_file) in read_files.iter().enumerate() {
        let is_sketch = j >= total - n_sketch_files;
        let sketch = get_seq_sketch(args, read_file, is_sketch, ref_db, effective_genome_c, db_k);
        let Some(sketch) = sketch else {
            log::error!(
                "--apply-unknown: sample `{}` could not be loaded/sketched. Exiting.",
                read_file[0]
            );
            std::process::exit(1);
        };
        let seq_name = sketch
            .sample_name
            .clone()
            .unwrap_or_else(|| sketch.file_name.clone());
        // Two inputs printing the same Sample_file (e.g. the same --sample-name, or
        // the same read path from different directories) map to indistinguishable
        // TSV rows, so their scalars would silently overwrite each other and rescale
        // one sample's rows with the other's read-length/count distribution. A real
        // -u run processes each input separately; refuse rather than mis-scale.
        if scalars
            .insert(seq_name.clone(), unknown_scalars_for_sample(args, &sketch))
            .is_some()
        {
            log::error!(
                "--apply-unknown: two sample inputs print the same Sample_file `{}`; their rows in the profile TSV cannot be told apart. Give each a distinct sample so their -u scalars do not collide. Exiting.",
                seq_name
            );
            std::process::exit(1);
        }
    }

    // Read the whole TSV (profiles are small) so we can sum covered bases per
    // sample before rescaling.
    let file = File::open(tsv_path).unwrap_or_else(|e| {
        log::error!(
            "--apply-unknown: could not open profile TSV `{}`: {}. Exiting.",
            tsv_path,
            e
        );
        std::process::exit(1);
    });
    let reader: Box<dyn BufRead> = if tsv_path.ends_with(".gz") {
        Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };
    let mut lines = reader.lines();

    let header = match lines.next() {
        Some(Ok(h)) => h,
        _ => {
            log::error!(
                "--apply-unknown: profile TSV `{}` is empty. Exiting.",
                tsv_path
            );
            std::process::exit(1);
        }
    };
    let cols: Vec<&str> = header.split('\t').collect();
    // Column layout of a pseudotax profile (see print_header/print_ani_result).
    const N_COLS: usize = 15;
    const SEQ_ABUND: usize = 3; // Sequence_abundance
    const COV: usize = 5; // Eff_cov (or True_cov)
    if cols.len() != N_COLS
        || cols[0] != "Sample_file"
        || cols[1] != "Genome_file"
        || cols[SEQ_ABUND] != "Sequence_abundance"
    {
        log::error!(
            "--apply-unknown: `{}` is not a `profile` TSV (expected {} columns starting Sample_file/Genome_file). Exiting.",
            tsv_path,
            N_COLS
        );
        std::process::exit(1);
    }
    if cols[COV] == "True_cov" {
        log::error!(
            "--apply-unknown: `{}` was already produced with -u (True_cov column); it must be a profile made WITHOUT -u. Exiting.",
            tsv_path
        );
        std::process::exit(1);
    }

    // Pass 1: keep rows, and sum gn_size * Eff_cov per sample.
    let rows: Vec<Vec<String>> = lines
        .map(|l| {
            l.unwrap_or_else(|e| {
                log::error!(
                    "--apply-unknown: error reading `{}`: {}. Exiting.",
                    tsv_path,
                    e
                );
                std::process::exit(1);
            })
        })
        .filter(|l| !l.is_empty())
        .map(|l| l.split('\t').map(|s| s.to_string()).collect())
        .collect();

    let mut covered_bases: FxHashMap<String, f64> = FxHashMap::default();
    for row in rows.iter() {
        if row.len() != N_COLS {
            log::error!(
                "--apply-unknown: malformed row in `{}` (expected {} columns): {:?}. Exiting.",
                tsv_path,
                N_COLS,
                row
            );
            std::process::exit(1);
        }
        let sample = &row[0];
        let genome = &row[1];
        if !scalars.contains_key(sample) {
            log::error!(
                "--apply-unknown: TSV references sample `{}`, but no matching sample sketch was given on the command line. Exiting.",
                sample
            );
            std::process::exit(1);
        }
        let gn_size = *gn_size_map.get(genome).unwrap_or_else(|| {
            log::error!(
                "--apply-unknown: genome `{}` from the TSV is not in the given database, so its size is unknown. Exiting.",
                genome
            );
            std::process::exit(1);
        });
        let eff_cov: f64 = row[COV].parse().unwrap_or_else(|_| {
            log::error!(
                "--apply-unknown: could not parse Eff_cov value `{}` in `{}`. Exiting.",
                row[COV],
                tsv_path
            );
            std::process::exit(1);
        });
        *covered_bases.entry(sample.clone()).or_insert(0.) += gn_size as f64 * eff_cov;
    }

    // bases_explained per sample = min( cov_scale * covered_bases / tentative_bases, 1 ).
    let mut bases_explained: FxHashMap<String, f64> = FxHashMap::default();
    for (sample, sc) in scalars.iter() {
        let covered = *covered_bases.get(sample).unwrap_or(&0.);
        let be = if sc.tentative_bases == 0. {
            0.
        } else {
            f64::min(sc.cov_scale * covered / sc.tentative_bases, 1.)
        };
        log::info!(
            "{} has {:.2}% of reads detected in database by profile",
            sample,
            be * 100.
        );
        bases_explained.insert(sample.clone(), be);
    }

    // Pass 2: rescale and print. True_cov = Eff_cov * cov_scale; Sequence_abundance
    // *= bases_explained; every other column is passed through verbatim.
    let mut new_header: Vec<&str> = cols.clone();
    new_header[COV] = "True_cov";
    let mut w = out_writer.lock().unwrap();
    writeln!(w, "{}", new_header.join("\t")).expect("Error writing to file");
    for mut row in rows {
        let sample = row[0].clone();
        let sc = &scalars[&sample];
        let be = bases_explained[&sample];
        let eff_cov: f64 = row[COV].parse().unwrap();
        let seq_abund: f64 = row[SEQ_ABUND].parse().unwrap();
        row[COV] = format!("{:.3}", eff_cov * sc.cov_scale);
        row[SEQ_ABUND] = format!("{:.4}", seq_abund * be);
        writeln!(w, "{}", row.join("\t")).expect("Error writing to file");
    }
}

fn estimate_covered_bases(
    results: &Vec<AniResult>,
    sequence_sketch: &SequencesSketch,
    read_length: f64,
    k: usize,
) -> f64 {
    let multiplier = read_length / (read_length - (k as f64) + 1.);

    let mut num_covered_bases = 0.;
    for res in results.iter() {
        num_covered_bases += (res.genome_sketch.gn_size as f64) * res.final_est_cov
    }
    let mut num_total_counts = 0;
    for count in sequence_sketch.kmer_counts.values() {
        num_total_counts += *count as usize;
    }
    let num_tentative_bases = sequence_sketch.c * num_total_counts;
    let num_tentative_bases = num_tentative_bases as f64 * multiplier;
    if num_tentative_bases == 0. {
        return 0.;
    }
    f64::min(num_covered_bases / num_tentative_bases, 1.)
}

fn winner_table<'a>(
    results: &'a Vec<AniResult>,
    log_reassign: bool,
) -> FxHashMap<Kmer, (f64, &'a GenomeSketch, bool)> {
    let mut kmer_to_genome_map: FxHashMap<_, _> = FxHashMap::default();
    for res in results.iter() {
        //let gn_sketch = &genome_sketches[res.genome_sketch_index];
        let gn_sketch = res.genome_sketch;
        for kmer in gn_sketch.genome_kmers.iter() {
            let v = kmer_to_genome_map.entry(*kmer).or_insert((
                res.final_est_ani,
                res.genome_sketch,
                false,
            ));
            if res.final_est_ani > v.0 {
                *v = (res.final_est_ani, gn_sketch, true);
            }
        }

        if gn_sketch.pseudotax_tracked_nonused_kmers.is_some() {
            for kmer in gn_sketch
                .pseudotax_tracked_nonused_kmers
                .as_ref()
                .unwrap()
                .iter()
            {
                let v = kmer_to_genome_map.entry(*kmer).or_insert((
                    res.final_est_ani,
                    res.genome_sketch,
                    false,
                ));
                if res.final_est_ani > v.0 {
                    *v = (res.final_est_ani, gn_sketch, true);
                }
            }
        }
    }

    //log reassigned kmers
    if log_reassign {
        log::info!("------------- Logging k-mer reassignments -----------------");
        let mut sketch_to_index = FxHashMap::default();
        for (i, res) in results.iter().enumerate() {
            log::info!(
                "Index\t{}\t{}\t{}",
                i,
                res.genome_sketch.file_name,
                res.genome_sketch.first_contig_name
            );
            sketch_to_index.insert(res.genome_sketch, i);
        }
        (0..results.len()).into_par_iter().for_each(|i| {
            let res = &results[i];
            let mut reassign_edge_map = FxHashMap::default();
            for kmer in res.genome_sketch.genome_kmers.iter() {
                let value = kmer_to_genome_map[kmer].1;
                if value != res.genome_sketch {
                    let edge_count = reassign_edge_map
                        .entry((sketch_to_index[value], i))
                        .or_insert(0);
                    *edge_count += 1;
                }
            }
            for (key, val) in reassign_edge_map {
                if val > 10 {
                    log::info!("{}->{}\t{}\tkmers reassigned", key.0, key.1, val);
                }
            }
        });
    }

    kmer_to_genome_map
}

fn print_header(pseudotax: bool, writer: &mut Box<dyn Write + Send>, estimate_unknown: bool) {
    if !pseudotax {
        writeln!(writer,
            //"Sample_file\tQuery_file\tAdjusted_ANI\tNaive_ANI\tANI_5-95_percentile\tEff_cov\tEff_lambda\tLambda_5-95_percentile\tMedian_cov\tMean_cov_geq1\tContainment_ind\tContig_name",
            "Sample_file\tGenome_file\tAdjusted_ANI\tEff_cov\tANI_5-95_percentile\tEff_lambda\tLambda_5-95_percentile\tMedian_cov\tMean_cov_geq1\tContainment_ind\tNaive_ANI\tContig_name",
            ).expect("Error writing to file.");
    } else {
        let cov_head;
        if estimate_unknown {
            cov_head = "True_cov";
        } else {
            cov_head = "Eff_cov";
        }
        writeln!(writer,
            "Sample_file\tGenome_file\tTaxonomic_abundance\tSequence_abundance\tAdjusted_ANI\t{}\tANI_5-95_percentile\tEff_lambda\tLambda_5-95_percentile\tMedian_cov\tMean_cov_geq1\tContainment_ind\tNaive_ANI\tkmers_reassigned\tContig_name", cov_head
            ).expect("Error writing to file.");
    }
}

fn get_genome_sketches(
    args: &ContainArgs,
    genome_sketch_files: &Vec<&String>,
    genome_files: &Vec<&String>,
) -> Vec<GenomeSketch> {
    let mut lowest_genome_c = None;
    let mut current_k = None;

    let genome_sketches = Mutex::new(vec![]);

    for genome_sketch_file in genome_sketch_files {
        let file = File::open(genome_sketch_file).unwrap_or_else(|_| {
            panic!(
                "The sketch `{}` could not be opened. Exiting",
                genome_sketch_file
            )
        });
        let mut genome_reader = BufReader::with_capacity(10_000_000, file);
        let genome_sketches_vec: Vec<GenomeSketch> = if crate::compress::peek_is_compressed(
            &mut genome_reader,
        )
        .unwrap_or(false)
        {
            crate::compress::read_genome_sketches_compressed(&mut genome_reader)
                .unwrap_or_else(|_| panic!("The sketch `{}` is not a valid sketch. Perhaps it is an older, incompatible version ",
                    &genome_sketch_file))
        } else {
            bincode::deserialize_from(&mut genome_reader)
                .unwrap_or_else(|_| panic!("The sketch `{}` is not a valid sketch. Perhaps it is an older, incompatible version ",
                    &genome_sketch_file))
        };
        if genome_sketches_vec.is_empty() {
            continue;
        }
        let c = genome_sketches_vec.first().unwrap().c;
        let k = genome_sketches_vec.first().unwrap().k;
        if lowest_genome_c.is_none() {
            lowest_genome_c = Some(c);
        } else if lowest_genome_c.unwrap() < c {
            lowest_genome_c = Some(c);
        }
        if current_k.is_none() {
            current_k = Some(genome_sketches_vec.first().unwrap().k);
        } else if current_k.unwrap() != k {
            error!("Query sketches have inconsistent -k. Exiting.");
            std::process::exit(1);
        }
        genome_sketches.lock().unwrap().extend(genome_sketches_vec);
    }

    genome_files.into_par_iter().for_each(|genome_file|{
        if lowest_genome_c.is_some() && lowest_genome_c.unwrap() < args.c{
            error!("Value of -c for contain is {} -- greater than the smallest value of -c for a genome sketch {}. Continuing without sketching.", args.c, lowest_genome_c.unwrap());
        }
        else if current_k.is_some() && current_k.unwrap() != args.k{
            error!("-k {} is not equal to -k {} found in sketches. Continuing without sketching.", args.k, current_k.unwrap());
        }
        else {
            if args.individual{
            let indiv_gn_sketches = sketch_genome_individual(args.c, args.k, genome_file, args.min_spacing_kmer, args.pseudotax);
                genome_sketches.lock().unwrap().extend(indiv_gn_sketches);

            }
            else{
                let genome_sketch_opt = sketch_genome(args.c, args.k, genome_file, args.min_spacing_kmer, args.pseudotax);
                if genome_sketch_opt.is_some() {
                    genome_sketches.lock().unwrap().push(genome_sketch_opt.unwrap());
                }
            }
        }
    });

    genome_sketches.into_inner().unwrap()
}

fn get_seq_sketch(
    args: &ContainArgs,
    read_file: &Vec<&String>,
    is_sketch_file: bool,
    ref_db: Option<&crate::refdelta::RefIndex>,
    genome_c: usize,
    genome_k: usize,
) -> Option<SequencesSketch> {
    get_seq_sketch_with_meta(args, read_file, is_sketch_file, ref_db, genome_c, genome_k)
        .map(|(sketch, _)| sketch)
}

/// Like `get_seq_sketch`, but also returns the sample's `ReadSketchMeta` (read count).
/// Used by `--merge`, where the merged `mean_read_length` is weighted by read counts.
fn get_seq_sketch_with_meta(
    args: &ContainArgs,
    read_file: &Vec<&String>,
    is_sketch_file: bool,
    ref_db: Option<&crate::refdelta::RefIndex>,
    genome_c: usize,
    genome_k: usize,
) -> Option<(SequencesSketch, ReadSketchMeta)> {
    if is_sketch_file {
        let read_file = read_file[0];
        let read_sketch_file = read_file;
        let (read_sketch, meta): (SequencesSketch, ReadSketchMeta) = if read_sketch_file
            .ends_with(REF_SAMPLE_SUFFIX)
        {
            let db = ref_db.unwrap_or_else(|| panic!(
                "`{}` is a reference-delta compressed sample (*.sylspr) but no --reference was provided",
                read_sketch_file
            ));
            let file = File::open(read_sketch_file).unwrap_or_else(|_| {
                panic!("The sketch `{}` could not be opened", &read_sketch_file)
            });
            let mut read_reader = BufReader::with_capacity(10_000_000, file);
            crate::refdelta::decompress_seq_with_meta(&mut read_reader, db).unwrap_or_else(|e| {
                panic!(
                    "Could not decode `{}` with the given --reference: {}",
                    read_sketch_file, e
                )
            })
        } else {
            let file = File::open(read_sketch_file).unwrap_or_else(|_| {
                panic!("The sketch `{}` could not be opened", &read_sketch_file)
            });
            let mut read_reader = BufReader::with_capacity(10_000_000, file);
            if crate::compress::peek_is_compressed(&mut read_reader).unwrap_or(false) {
                crate::compress::read_seq_sketch_compressed_with_meta(&mut read_reader)
                    .unwrap_or_else(|e| {
                        panic!("The sketch `{}` could not be read: {}", read_sketch_file, e)
                    })
            } else {
                // Legacy uncompressed *.sylsp samples carry no recorded read count. Merging
                // weights each input's read length by its read count, so a legacy sample would
                // silently contribute a read length of zero and skew the merged sample's
                // read-count/abundance estimates. Refuse it in --merge mode; per-sample
                // profiling (the non-merge path) still accepts it unchanged.
                if args.merge {
                    error!(
                        "`{}` is a legacy uncompressed sample sketch (*.sylsp) with no recorded read count, so it cannot be merged (read length cannot be determined). Re-sketch it with the current version, or exclude it from --merge. Exiting.",
                        read_sketch_file
                    );
                    std::process::exit(1);
                }
                let sketch: SequencesSketch = bincode::deserialize_from(&mut read_reader).unwrap_or_else(|_| panic!("The sketch `{}` is not a valid sketch. Perhaps it is an older incompatible version ", read_sketch_file));
                (sketch, ReadSketchMeta::default())
            }
        };
        if read_sketch.c > genome_c {
            error!("{} value of -c is {}; this is greater than the smallest value of -c = {} for a genome sketch. Exiting.", read_file, read_sketch.c, genome_c);
            return None;
        } else if read_sketch.c < genome_c {
            info!("{} value of -c for reads is {}; this is smaller than the -c for a genome sketch. Using the larger -c {} instead.", read_file, read_sketch.c,  genome_c);
        }

        Some((read_sketch, meta))
    } else {
        if args.c > genome_c {
            info!("{} value of -c for reads is {}; this is smaller than the -c for a genome sketch. Using the larger -c {} instead.", read_file[0], args.c,  genome_c);
        }
        if genome_c < args.c {
            error!("{} error: value of -c for contain = {} -- greater than the smallest value of -c for a genome sketch = {}. Continuing without sketching.", read_file[0], args.c, genome_c);
            None
        } else if genome_k != args.k {
            error!(
                "{} -k {} is not equal to -k {} found in sketches. Continuing without sketching.",
                read_file[0], args.k, genome_k
            );
            None
        } else if read_file.len() == 1 {
            sketch_sequences_needle(read_file[0], args.c, args.k, None, false, false)
        } else if read_file.len() == 2 {
            sketch_pair_sequences(
                read_file[0],
                read_file[1],
                args.c,
                args.k,
                None,
                false,
                DEFAULT_FPR,
                false,
            )
        } else if read_file.len() == 3 {
            sketch_interleaved_sequences(
                read_file[0],
                args.c,
                args.k,
                None,
                false,
                DEFAULT_FPR,
                false,
            )
        } else {
            panic!(
                "Internal Error: read_file has length {}. Something went wrong...",
                read_file.len()
            );
        }
    }
}

fn get_stats<'a>(
    args: &ContainArgs,
    genome_sketch: &'a GenomeSketch,
    sequence_sketch: &SequencesSketch,
    winner_map: Option<&FxHashMap<Kmer, (f64, &GenomeSketch, bool)>>,
    log_reassign: bool,
) -> Option<AniResult<'a>> {
    if genome_sketch.k != sequence_sketch.k {
        log::error!(
            "k parameter for reads {} != k parameter for genome {}",
            sequence_sketch.k,
            genome_sketch.k
        );
        std::process::exit(1);
    }
    if genome_sketch.c < sequence_sketch.c {
        log::error!(
            "c parameter for reads {} > c parameter for genome {}",
            sequence_sketch.c,
            genome_sketch.c
        );
        std::process::exit(1);
    }
    let mut contain_count = 0;
    let mut covs = vec![];
    let gn_kmers = &genome_sketch.genome_kmers;
    if (gn_kmers.len() as f64) < args.min_number_kmers {
        return None;
    }

    let mut kmers_lost_count = 0;
    for kmer in gn_kmers.iter() {
        if sequence_sketch.kmer_counts.contains_key(kmer) {
            if sequence_sketch.kmer_counts[kmer] == 0 {
                continue;
            }
            if winner_map.is_some() {
                let map = &winner_map.unwrap();
                if map[kmer].1 != genome_sketch {
                    kmers_lost_count += 1;
                    continue;
                }
                contain_count += 1;
                covs.push(sequence_sketch.kmer_counts[kmer]);
            } else {
                contain_count += 1;
                covs.push(sequence_sketch.kmer_counts[kmer]);
            }
        }
    }

    let n_kmers = gn_kmers.len();
    let reassign_log = if winner_map.is_some() && log_reassign {
        Some((
            genome_sketch.file_name.as_str(),
            genome_sketch.first_contig_name.as_str(),
            kmers_lost_count,
        ))
    } else {
        None
    };
    let fin = finalize_stats(
        args,
        genome_sketch.k,
        n_kmers,
        contain_count,
        covs,
        reassign_log,
    )?;

    let seq_name = if let Some(sample) = &sequence_sketch.sample_name {
        sample.clone()
    } else {
        sequence_sketch.file_name.clone()
    };
    let kmers_lost = if winner_map.is_some() {
        Some(kmers_lost_count)
    } else {
        None
    };

    Some(AniResult {
        naive_ani: fin.naive_ani,
        final_est_ani: fin.final_est_ani,
        final_est_cov: fin.final_est_cov,
        seq_name,
        gn_name: genome_sketch.file_name.as_str(),
        contig_name: genome_sketch.first_contig_name.as_str(),
        mean_cov: fin.mean_cov,
        median_cov: fin.median_cov,
        containment_index: (contain_count, n_kmers),
        lambda: fin.lambda,
        ani_ci: fin.ani_ci,
        lambda_ci: fin.lambda_ci,
        genome_sketch,
        rel_abund: None,
        seq_abund: None,
        kmers_lost,
    })
}

/// Scalar outputs of the coverage-correction + ANI estimation.
struct Finalized {
    naive_ani: f64,
    final_est_ani: f64,
    final_est_cov: f64,
    median_cov: f64,
    mean_cov: f64,
    lambda: AdjustStatus,
    ani_ci: (Option<f64>, Option<f64>),
    lambda_ci: (Option<f64>, Option<f64>),
}

/// Coverage-correction + ANI math shared by `get_stats` (materialized genome)
/// and the streaming two-stage pass-1 (which never builds a genome `Vec`).
/// Consumes the matched coverage counts `covs` and the genome k-mer count
/// `n_kmers`; applies the minimum-ANI gate (returns None below it). When
/// `reassign_log` is set, logs a drop during pseudotax reassignment.
fn finalize_stats(
    args: &ContainArgs,
    k: usize,
    n_kmers: usize,
    contain_count: usize,
    mut covs: Vec<u32>,
    reassign_log: Option<(&str, &str, usize)>,
) -> Option<Finalized> {
    if covs.is_empty() {
        return None;
    }
    let naive_ani = f64::powf(contain_count as f64 / n_kmers as f64, 1. / k as f64);
    covs.sort();
    let median_cov = covs[covs.len() / 2] as f64;
    let pois = Poisson::new(median_cov).unwrap();
    let mut max_cov = f64::MAX;
    if median_cov < 30. {
        for i in covs.len() / 2..covs.len() {
            let cov = covs[i];
            if pois.cdf(cov.into()) < CUTOFF_PVALUE {
                max_cov = cov as f64;
            } else {
                break;
            }
        }
    }

    let mut full_covs = vec![0; n_kmers - contain_count];
    for cov in covs.iter() {
        if (*cov as f64) <= max_cov {
            full_covs.push(*cov);
        }
    }
    let mean_cov = full_covs.iter().sum::<u32>() as f64 / full_covs.len() as f64;
    let geq1_mean_cov = full_covs.iter().sum::<u32>() as f64 / covs.len() as f64;

    let use_lambda;
    if median_cov > MEDIAN_ANI_THRESHOLD {
        use_lambda = AdjustStatus::High
    } else {
        let test_lambda;
        if args.ratio {
            test_lambda = ratio_lambda(&full_covs, args.min_count_correct)
        } else if args.mme {
            test_lambda = mme_lambda(&full_covs)
        } else if args.nb {
            test_lambda = binary_search_lambda(&full_covs)
        } else if args.mle {
            test_lambda = mle_zip(&full_covs, k as f64)
        } else {
            test_lambda = ratio_lambda(&full_covs, args.min_count_correct)
        };
        if test_lambda.is_none() {
            use_lambda = AdjustStatus::Low
        } else {
            use_lambda = AdjustStatus::Lambda(test_lambda.unwrap());
        }
    }

    let final_est_cov;
    if let AdjustStatus::Lambda(lam) = use_lambda {
        final_est_cov = lam
    } else if median_cov < MAX_MEDIAN_FOR_MEAN_FINAL_EST {
        final_est_cov = geq1_mean_cov;
    } else if args.mean_coverage {
        final_est_cov = geq1_mean_cov;
    } else {
        final_est_cov = median_cov;
    }

    let opt_lambda;
    if use_lambda == AdjustStatus::Low || use_lambda == AdjustStatus::High {
        opt_lambda = None
    } else {
        opt_lambda = Some(final_est_cov)
    };

    let opt_est_ani = ani_from_lambda(opt_lambda, mean_cov, k as f64, &full_covs);

    let final_est_ani;
    if opt_lambda.is_none() || opt_est_ani.is_none() || args.no_adj {
        final_est_ani = naive_ani;
    } else {
        final_est_ani = opt_est_ani.unwrap();
    }

    let min_ani = if args.minimum_ani.is_some() {
        args.minimum_ani.unwrap() / 100.
    } else if args.pseudotax {
        MIN_ANI_P_DEF
    } else {
        MIN_ANI_DEF
    };
    if final_est_ani < min_ani {
        if let Some((gn, ctg, lost)) = reassign_log {
            log::info!(
                "Genome/contig {}/{} has ANI = {} < {} after reassigning {} k-mers ({} contained k-mers after reassign)",
                gn, ctg, final_est_ani * 100., min_ani * 100., lost, contain_count
            );
        }
        return None;
    }

    let (mut low_ani, mut high_ani, mut low_lambda, mut high_lambda) = (None, None, None, None);
    if !args.no_ci && opt_lambda.is_some() {
        let bootstrap = bootstrap_interval(&full_covs, k as f64, args);
        low_ani = bootstrap.0;
        high_ani = bootstrap.1;
        low_lambda = bootstrap.2;
        high_lambda = bootstrap.3;
    }

    Some(Finalized {
        naive_ani,
        final_est_ani,
        final_est_cov,
        median_cov,
        mean_cov: geq1_mean_cov,
        lambda: use_lambda,
        ani_ci: (low_ani, high_ani),
        lambda_ci: (low_lambda, high_lambda),
    })
}

fn ani_from_lambda(lambda: Option<f64>, _mean: f64, k: f64, full_cov: &[u32]) -> Option<f64> {
    lambda?;
    let mut contain_count = 0;
    let mut _zero_count = 0;
    for x in full_cov {
        if *x != 0 {
            contain_count += 1;
        } else {
            _zero_count += 1;
        }
    }

    let lambda = lambda.unwrap();
    let adj_index = contain_count as f64 / (1. - f64::exp(-lambda)) / full_cov.len() as f64;
    let ret_ani;
    //let ani = f64::powf(1. - pi, 1./k);
    let ani = f64::powf(adj_index, 1. / k);
    if ani < 0. || ani.is_nan() {
        ret_ani = None;
    } else {
        if ani > 1. {
            ret_ani = Some(ani)
        } else {
            ret_ani = Some(ani);
        }
    }
    ret_ani
}

fn bootstrap_interval(
    covs_full: &Vec<u32>,
    k: f64,
    args: &ContainArgs,
) -> (Option<f64>, Option<f64>, Option<f64>, Option<f64>) {
    fastrand::seed(7);
    let num_samp = covs_full.len();
    let iters = 100;
    let mut res_ani = vec![];
    let mut res_lambda = vec![];

    for _ in 0..iters {
        let mut rand_vec = vec![];
        rand_vec.reserve(num_samp);
        for _ in 0..num_samp {
            rand_vec.push(covs_full[fastrand::usize(..covs_full.len())]);
        }
        let lambda;
        if args.ratio {
            lambda = ratio_lambda(&rand_vec, args.min_count_correct);
        } else if args.mme {
            lambda = mme_lambda(&rand_vec);
        } else if args.nb {
            lambda = binary_search_lambda(&rand_vec);
        } else if args.mle {
            lambda = mle_zip(&rand_vec, k);
        } else {
            lambda = ratio_lambda(&rand_vec, args.min_count_correct);
        }
        let ani = ani_from_lambda(lambda, mean(&rand_vec).unwrap(), k, &rand_vec);
        if ani.is_some() && lambda.is_some() && !ani.unwrap().is_nan() && !lambda.unwrap().is_nan()
        {
            res_ani.push(ani);
            res_lambda.push(lambda);
        }
    }
    res_ani.sort_by(|x, y| x.partial_cmp(y).unwrap());
    res_lambda.sort_by(|x, y| x.partial_cmp(y).unwrap());
    if res_ani.len() < 50 {
        return (None, None, None, None);
    }
    let suc = res_ani.len();
    let low_ani = res_ani[suc * 5 / 100 - 1];
    let high_ani = res_ani[suc * 95 / 100 - 1];
    let low_lambda = res_lambda[suc * 5 / 100 - 1];
    let high_lambda = res_lambda[suc * 95 / 100 - 1];

    (low_ani, high_ani, low_lambda, high_lambda)
}

fn get_kmer_identity(seq_sketch: &SequencesSketch, estimate_unknown: bool) -> Option<f64> {
    if !estimate_unknown {
        return None;
    }

    let mut median = 0;
    let mut mov_avg_median = 0.;
    let mut n = 1.;
    for count in seq_sketch.kmer_counts.values() {
        if *count > 1 {
            if *count > median {
                median += 1;
            } else {
                median -= 1;
            }
            mov_avg_median += median as f64;
            n += 1.;
        }
    }

    mov_avg_median /= n;
    log::debug!(
        "Estimated continuous median k-mer count for {} is {:.3}",
        &seq_sketch.file_name,
        mov_avg_median
    );

    let mut num_1s = 0;
    let mut num_not1s = 0;
    for count in seq_sketch.kmer_counts.values() {
        if *count == 1 {
            num_1s += 1;
        } else {
            num_not1s += *count;
        }
    }
    //0.1 so no div by 0 error
    let eps = num_not1s as f64 / (num_not1s as f64 + num_1s as f64 + 0.1);
    //dbg!("Automatic id est, 1-to-2 ratio, 2-to-3", eps.powf(1./31.), num_1s as f64 / num_2s as f64, two_to_three);

    if mov_avg_median < MED_KMER_FOR_ID_EST && seq_sketch.mean_read_length < 400. {
        log::info!("{} short-read sample has high diversity compared to sequencing depth (approx. avg depth < 3). Using 99.5% as read accuracy estimate instead of automatic detection for --estimate-unknown.", &seq_sketch.file_name);
        return Some(0.995f64.powf(seq_sketch.k as f64));
    }

    if eps < 1. {
        Some(eps)
    } else {
        Some(1.)
    }
}
