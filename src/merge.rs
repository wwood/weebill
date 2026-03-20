use crate::cmdline::MergeArgs;
use crate::constants::*;
use crate::types::*;
use log::*;
use std::fs::File;
use std::io::{BufReader, BufWriter};

pub fn merge(args: MergeArgs) {
    simple_logger::SimpleLogger::new()
        .with_level(log::LevelFilter::Info)
        .init()
        .unwrap();

    if args.files.len() < 2 {
        error!("Need at least 2 sketch files to merge. Exiting.");
        std::process::exit(1);
    }

    let mut sketches: Vec<SequencesSketch> = Vec::new();

    for file in args.files.iter() {
        let mut valid_suffix = false;
        for suff in SAMPLE_FILE_SUFFIX_VALID {
            if file.ends_with(suff) {
                valid_suffix = true;
                break;
            }
        }
        if !valid_suffix {
            error!(
                "'{}' does not have a valid sample sketch suffix (.sylsp). Exiting.",
                file
            );
            std::process::exit(1);
        }

        let f = File::open(file).expect(&format!(
            "The sketch '{}' could not be opened. Exiting.",
            file
        ));
        let reader = BufReader::with_capacity(10_000_000, f);
        let sketch: SequencesSketch = bincode::deserialize_from(reader).expect(&format!(
            "The sketch '{}' is not a valid sketch. Perhaps it is an older, incompatible version.",
            file
        ));
        sketches.push(sketch);
    }

    // Validate all sketches have same c and k
    let c = sketches[0].c;
    let k = sketches[0].k;
    for (i, sk) in sketches.iter().enumerate() {
        if sk.c != c || sk.k != k {
            error!(
                "Sketch '{}' has c={}, k={} but expected c={}, k={}. All sketches must have the same c and k. Exiting.",
                args.files[i], sk.c, sk.k, c, k
            );
            std::process::exit(1);
        }
    }

    let merged = merge_sketches(&sketches, args.sample_name);

    let output_path = if args.output.ends_with(SAMPLE_FILE_SUFFIX) {
        args.output.clone()
    } else {
        format!("{}{}", args.output, SAMPLE_FILE_SUFFIX)
    };

    let mut writer = BufWriter::new(
        File::create(&output_path).expect(&format!("Could not create output file '{}'", output_path)),
    );
    bincode::serialize_into(&mut writer, &merged).unwrap();
    info!(
        "Merged {} sketches into '{}' ({} k-mers, {} total reads).",
        sketches.len(),
        output_path,
        merged.kmer_counts.len(),
        merged.num_reads,
    );
}

pub fn merge_sketches(sketches: &[SequencesSketch], sample_name: Option<String>) -> SequencesSketch {
    assert!(sketches.len() >= 2);

    let mut merged_counts = sketches[0].kmer_counts.clone();
    let c = sketches[0].c;
    let k = sketches[0].k;
    let mut total_reads: u64 = sketches[0].num_reads;
    let mut weighted_length: f64 = sketches[0].mean_read_length * sketches[0].num_reads as f64;
    let mut any_paired = sketches[0].paired;

    for sk in sketches[1..].iter() {
        for (kmer, count) in sk.kmer_counts.iter() {
            *merged_counts.entry(*kmer).or_insert(0) += count;
        }
        total_reads += sk.num_reads;
        weighted_length += sk.mean_read_length * sk.num_reads as f64;
        any_paired = any_paired || sk.paired;
    }

    let mean_read_length = if total_reads > 0 {
        weighted_length / total_reads as f64
    } else {
        0.0
    };

    SequencesSketch {
        kmer_counts: merged_counts,
        c,
        k,
        file_name: sketches[0].file_name.clone(),
        sample_name,
        paired: any_paired,
        mean_read_length,
        num_reads: total_reads,
    }
}
