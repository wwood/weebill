use crate::cmdline::MergeArgs;
use crate::constants::*;
use crate::refdelta;
use crate::types::*;
use log::*;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::Path;

fn open_reference(path: &str) -> refdelta::RefIndex {
    let file =
        File::open(path).unwrap_or_else(|_| panic!("Could not open reference database {}", path));
    refdelta::open_ref_index_file(file)
        .unwrap_or_else(|e| panic!("{} is not a valid reference database: {}", path, e))
}

fn read_sample(
    path: &str,
    ref_index: Option<&refdelta::RefIndex>,
) -> (SequencesSketch, ReadSketchMeta) {
    let file =
        File::open(path).unwrap_or_else(|_| panic!("The sketch '{}' could not be opened", path));
    let mut reader = BufReader::with_capacity(10_000_000, file);

    if path.ends_with(REF_SAMPLE_SUFFIX) {
        let idx = ref_index.unwrap_or_else(|| {
            panic!(
                "`{}` is a reference-compressed sample (*.sylspr); pass --reference",
                path
            )
        });
        refdelta::decompress_seq_with_meta(&mut reader, idx).unwrap_or_else(|e| {
            panic!(
                "Could not decode `{}` with the given --reference: {}",
                path, e
            )
        })
    } else if crate::compress::peek_is_compressed(&mut reader).unwrap_or(false) {
        crate::compress::read_seq_sketch_compressed_with_meta(&mut reader)
            .unwrap_or_else(|e| panic!("The compressed sketch `{}` is invalid: {}", path, e))
    } else {
        // Legacy uncompressed *.sylsp samples predate the ReadSketchMeta side-channel and
        // therefore carry no recorded read count. Merging weights each sample's read length
        // by its read count, so a legacy input would silently contribute a read length of
        // zero and corrupt the merged mean_read_length. Refuse it rather than merge garbage;
        // the user can re-sketch to the current format first.
        error!(
            "'{}' is a legacy uncompressed sample sketch (*.sylsp) with no recorded read count, so it cannot be merged (read length cannot be determined). Re-sketch it with the current version, or exclude it. Exiting.",
            path
        );
        std::process::exit(1);
    }
}

fn output_kind(args: &MergeArgs) -> &'static str {
    if args.ref_compress || args.output.ends_with(REF_SAMPLE_SUFFIX) {
        REF_SAMPLE_SUFFIX
    } else if args.compressed || args.output.ends_with(SAMPLE_COMP_FILE_SUFFIX) {
        SAMPLE_COMP_FILE_SUFFIX
    } else {
        SAMPLE_FILE_SUFFIX
    }
}

fn output_path(path: &str, suffix: &str) -> String {
    if path.ends_with(SAMPLE_FILE_SUFFIX)
        || path.ends_with(SAMPLE_COMP_FILE_SUFFIX)
        || path.ends_with(REF_SAMPLE_SUFFIX)
    {
        path.to_string()
    } else {
        format!("{}{}", path, suffix)
    }
}

pub fn merge(args: MergeArgs) {
    simple_logger::SimpleLogger::new()
        .with_level(log::LevelFilter::Info)
        .init()
        .unwrap();

    if args.files.len() < 2 {
        error!("Need at least 2 sketch files to merge. Exiting.");
        std::process::exit(1);
    }

    let needs_reference = args.ref_compress
        || args.output.ends_with(REF_SAMPLE_SUFFIX)
        || args.files.iter().any(|f| f.ends_with(REF_SAMPLE_SUFFIX));
    let ref_index = if needs_reference {
        Some(open_reference(args.reference.as_ref().unwrap_or_else(
            || {
                error!("--reference is required for *.sylspr input or output. Exiting.");
                std::process::exit(1);
            },
        )))
    } else {
        None
    };

    let sketches: Vec<(SequencesSketch, ReadSketchMeta)> = args
        .files
        .iter()
        .map(|f| read_sample(f, ref_index.as_ref()))
        .collect();

    let c = sketches[0].0.c;
    let k = sketches[0].0.k;
    for (i, (sketch, _)) in sketches.iter().enumerate() {
        if sketch.c != c || sketch.k != k {
            error!(
                "Sketch '{}' has c={}, k={} but expected c={}, k={}. All sketches must have the same c and k. Exiting.",
                args.files[i], sketch.c, sketch.k, c, k
            );
            std::process::exit(1);
        }
    }

    let (merged, meta) = merge_sketches(&sketches, args.sample_name.clone());
    let suffix = output_kind(&args);
    let out = output_path(&args.output, suffix);
    if let Some(parent) = Path::new(&out).parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)
                .unwrap_or_else(|_| panic!("Could not create output directory for '{}'", out));
        }
    }
    let mut writer =
        BufWriter::new(File::create(&out).unwrap_or_else(|_| panic!("Could not create '{}'", out)));

    match suffix {
        REF_SAMPLE_SUFFIX => {
            let idx = ref_index
                .as_ref()
                .expect("reference index must be loaded for *.sylspr output");
            refdelta::compress_seq_with_meta(
                &mut writer,
                &merged,
                idx,
                args.reference.as_deref().unwrap_or(""),
                meta,
            )
            .unwrap_or_else(|e| panic!("Could not write reference-compressed output: {}", e));
        }
        SAMPLE_COMP_FILE_SUFFIX => {
            crate::compress::write_seq_sketch_compressed_with_meta(&mut writer, &merged, meta)
                .unwrap_or_else(|e| panic!("Could not write compressed output: {}", e));
        }
        _ => {
            bincode::serialize_into(&mut writer, &merged).unwrap();
        }
    }

    info!(
        "Merged {} sketches into '{}' ({} distinct k-mers, {} reads recorded).",
        sketches.len(),
        out,
        merged.kmer_counts.len(),
        meta.num_reads,
    );
}

pub fn merge_sketches(
    sketches: &[(SequencesSketch, ReadSketchMeta)],
    sample_name: Option<String>,
) -> (SequencesSketch, ReadSketchMeta) {
    assert!(sketches.len() >= 2);

    let c = sketches[0].0.c;
    let k = sketches[0].0.k;
    // Sketches subsampled at different rates (c) or built with different k cannot be
    // summed: their k-mer sets were selected against different thresholds, so a merged
    // count map would mix sampling rates and corrupt containment/coverage. Reject the
    // mismatch here so every caller (merge subcommand, sketch --merge, profile --merge)
    // is protected, not just the ones that pre-check.
    for (sketch, _) in sketches[1..].iter() {
        if sketch.c != c || sketch.k != k {
            error!(
                "Cannot merge sketches with differing parameters: '{}' has c={}, k={} but expected c={}, k={}. All merged inputs must share the same c and k. Exiting.",
                sketch.file_name, sketch.c, sketch.k, c, k
            );
            std::process::exit(1);
        }
    }

    let mut merged_counts = sketches[0].0.kmer_counts.clone();
    let mut total_reads = sketches[0].1.num_reads;
    let mut weighted_length = sketches[0].0.mean_read_length * sketches[0].1.num_reads as f64;
    let mut fallback_length_sum = sketches[0].0.mean_read_length;
    let mut any_paired = sketches[0].0.paired;

    for (sketch, meta) in sketches[1..].iter() {
        for (kmer, count) in sketch.kmer_counts.iter() {
            let slot = merged_counts.entry(*kmer).or_insert(0);
            *slot = slot.checked_add(*count).unwrap_or_else(|| {
                error!(
                    "k-mer count overflowed u32 while merging (a k-mer's summed count exceeds {}). Exiting.",
                    u32::MAX
                );
                std::process::exit(1);
            });
        }
        total_reads += meta.num_reads;
        weighted_length += sketch.mean_read_length * meta.num_reads as f64;
        fallback_length_sum += sketch.mean_read_length;
        any_paired = any_paired || sketch.paired;
    }

    let mean_read_length = if total_reads > 0 {
        weighted_length / total_reads as f64
    } else {
        fallback_length_sum / sketches.len() as f64
    };

    let file_name = sample_name
        .clone()
        .unwrap_or_else(|| sketches[0].0.file_name.clone());
    (
        SequencesSketch {
            kmer_counts: merged_counts,
            c,
            k,
            file_name,
            sample_name,
            paired: any_paired,
            mean_read_length,
        },
        ReadSketchMeta {
            num_reads: total_reads,
        },
    )
}
