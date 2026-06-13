use assert_cmd::prelude::*; // Add methods on commands
use sylph::seeding;
use sylph::compress;
use sylph::refdelta;
use sylph::types::{GenomeSketch, SequencesSketch};
use fxhash::FxHashMap;

fn gsketch(file_name: &str, kmers: Vec<u64>) -> GenomeSketch {
    GenomeSketch {
        genome_kmers: kmers,
        pseudotax_tracked_nonused_kmers: None,
        file_name: file_name.to_string(),
        first_contig_name: "c".to_string(),
        c: 200,
        k: 31,
        gn_size: 1000,
        min_spacing: 30,
    }
}

#[test]
fn refdelta_build_two_level_assignment() {
    // species A: rep A_rep + strain A_str; species B: rep B_rep
    let sketches = vec![
        gsketch("A_rep.fa", vec![1, 2, 3, 100, 101, 500]),
        gsketch("A_str.fa", vec![100, 101, 200, 201]),
        gsketch("B_rep.fa", vec![300, 301, 500]),
    ];
    let mut tax: FxHashMap<String, (String, bool)> = FxHashMap::default();
    tax.insert("A_rep.fa".into(), ("A".into(), true));
    tax.insert("A_str.fa".into(), ("A".into(), false));
    tax.insert("B_rep.fa".into(), ("B".into(), true));

    let db = refdelta::build_refdb(&sketches, &tax);
    // ordering: species A (rep first), then strain, then species B
    assert_eq!(db.genomes[0].file_name, "A_rep.fa");
    assert_eq!(db.genomes[1].file_name, "A_str.fa");
    assert_eq!(db.genomes[2].file_name, "B_rep.fa");
    // strains of A are contiguous (ids 0,1)
    assert!(db.genomes[0].is_rep && !db.genomes[1].is_rep);

    // 100,101 are in rep A_rep AND strain A_str -> rep wins (distinctive to A_rep)
    assert_eq!(db.distinctive[0], vec![1, 2, 3, 100, 101]);
    assert_eq!(db.distinctive[1], vec![200, 201]); // strain-only
    assert_eq!(db.distinctive[2], vec![300, 301]);
    // 500 is in two reps -> shared pool
    assert_eq!(db.pool, vec![500]);
}

fn refdelta_roundtrip(sketch: &SequencesSketch, db: &refdelta::RefDb) {
    let lookup = db.build_lookup();
    let mut buf = Vec::new();
    refdelta::compress_seq(&mut buf, sketch, db, &lookup).unwrap();
    let decoded = refdelta::decompress_seq(&buf[..], db).unwrap();
    assert_eq!(*sketch, decoded);
}

#[test]
fn refdelta_compress_decompress_roundtrip() {
    let sketches = vec![
        gsketch("A_rep.fa", vec![1, 2, 3, 100, 101, 500]),
        gsketch("A_str.fa", vec![100, 101, 200, 201]),
        gsketch("B_rep.fa", vec![300, 301, 500]),
    ];
    let mut tax: FxHashMap<String, (String, bool)> = FxHashMap::default();
    tax.insert("A_rep.fa".into(), ("A".into(), true));
    tax.insert("A_str.fa".into(), ("A".into(), false));
    tax.insert("B_rep.fa".into(), ("B".into(), true));
    let db = refdelta::build_refdb(&sketches, &tax);

    // reference DB serialization roundtrip
    let mut dbuf = Vec::new();
    refdelta::write_refdb(&mut dbuf, &db).unwrap();
    let db2 = refdelta::read_refdb(&dbuf[..]).unwrap();
    assert_eq!(db, db2);

    // distinctive hits + pool hit + novel hashes + counts
    let mut counts: FxHashMap<u64, u32> = FxHashMap::default();
    for (h, c) in [(1u64, 5u32), (2, 3), (100, 7), (200, 2), (500, 4), (999, 1), (98765, 9)] {
        counts.insert(h, c);
    }
    let sketch = SequencesSketch {
        kmer_counts: counts,
        c: 200,
        k: 31,
        file_name: "sample.fq".into(),
        sample_name: Some("s".into()),
        paired: true,
        mean_read_length: 149.0,
    };
    refdelta_roundtrip(&sketch, &db);

    // empty sketch and a sketch that is entirely novel also roundtrip
    refdelta_roundtrip(
        &SequencesSketch {
            kmer_counts: FxHashMap::default(),
            c: 200,
            k: 31,
            file_name: "empty.fq".into(),
            sample_name: None,
            paired: false,
            mean_read_length: 0.0,
        },
        &db,
    );
    let mut nov: FxHashMap<u64, u32> = FxHashMap::default();
    nov.insert(7_000_000, 1);
    nov.insert(8_000_000, 2);
    refdelta_roundtrip(
        &SequencesSketch {
            kmer_counts: nov,
            c: 200,
            k: 31,
            file_name: "novel.fq".into(),
            sample_name: None,
            paired: false,
            mean_read_length: 1.0,
        },
        &db,
    );
}

#[test]
fn refdelta_rejects_wrong_reference() {
    // Two DBs with the same shape (genome count, per-array lengths) and identical
    // boundary (first/last) hashes, differing only in an interior hash. A weak
    // fingerprint would collide; the full-content digest must reject decoding a
    // sample compressed against `db_a` using `db_b`.
    let mk = |middle: u64| {
        let sketches = vec![gsketch("g.fa", vec![10, middle, 30])];
        refdelta::build_refdb(&sketches, &FxHashMap::default())
    };
    let db_a = mk(20);
    let db_b = mk(21);
    assert_eq!(db_a.distinctive[0].len(), db_b.distinctive[0].len());
    assert_eq!(db_a.distinctive[0].first(), db_b.distinctive[0].first());
    assert_eq!(db_a.distinctive[0].last(), db_b.distinctive[0].last());
    assert_ne!(db_a.fingerprint, db_b.fingerprint);

    let mut counts: FxHashMap<u64, u32> = FxHashMap::default();
    counts.insert(10, 1);
    counts.insert(20, 2);
    let sketch = SequencesSketch {
        kmer_counts: counts,
        c: 200,
        k: 31,
        file_name: "s.fq".into(),
        sample_name: None,
        paired: false,
        mean_read_length: 1.0,
    };
    let lookup = db_a.build_lookup();
    let mut buf = Vec::new();
    refdelta::compress_seq(&mut buf, &sketch, &db_a, &lookup).unwrap();
    assert!(refdelta::decompress_seq(&buf[..], &db_a).is_ok());
    assert!(refdelta::decompress_seq(&buf[..], &db_b).is_err());
}

#[test]
fn compress_genome_roundtrip() {
    let mut sketches = Vec::new();
    for g in 0..3 {
        let mut genome_kmers: Vec<u64> = (0..5000u64)
            .map(|i| (i.wrapping_mul(2654435761).wrapping_add(g * 7)) % (u64::MAX / 200))
            .collect();
        // intentionally unsorted with a duplicate to exercise gap/dup handling
        genome_kmers.push(genome_kmers[0]);
        sketches.push(GenomeSketch {
            genome_kmers,
            pseudotax_tracked_nonused_kmers: if g % 2 == 0 {
                Some((0..1000u64).map(|i| i * 13 % (u64::MAX / 200)).collect())
            } else {
                None
            },
            file_name: format!("genome_{}.fa", g),
            first_contig_name: format!("contig {} description", g),
            c: 200,
            k: 31,
            gn_size: 4_600_000 + g as usize,
            min_spacing: 30,
        });
    }

    let mut buf = Vec::new();
    compress::write_genome_sketches_compressed(&mut buf, &sketches).unwrap();
    assert_eq!(&buf[..4], b"SYLZ", "compressed output must start with SYLZ magic");
    let decoded = compress::read_genome_sketches_compressed(&buf[..]).unwrap();

    assert_eq!(decoded.len(), sketches.len());
    for (orig, dec) in sketches.iter().zip(decoded.iter()) {
        assert_eq!(orig.file_name, dec.file_name);
        assert_eq!(orig.first_contig_name, dec.first_contig_name);
        assert_eq!(orig.c, dec.c);
        assert_eq!(orig.k, dec.k);
        assert_eq!(orig.gn_size, dec.gn_size);
        assert_eq!(orig.min_spacing, dec.min_spacing);

        let mut a = orig.genome_kmers.clone();
        a.sort_unstable();
        assert_eq!(a, dec.genome_kmers, "genome_kmers must match as a sorted multiset");

        match (&orig.pseudotax_tracked_nonused_kmers, &dec.pseudotax_tracked_nonused_kmers) {
            (Some(oa), Some(da)) => {
                let mut oa = oa.clone();
                oa.sort_unstable();
                assert_eq!(oa, *da);
            }
            (None, None) => {}
            _ => panic!("pseudotax option mismatch"),
        }
    }
}

#[test]
fn compress_seq_roundtrip() {
    let mut kmer_counts: FxHashMap<u64, u32> = FxHashMap::default();
    for i in 0..10000u64 {
        let key = (i.wrapping_mul(11400714819323198485)) % (u64::MAX / 200);
        kmer_counts.insert(key, (i % 4 + 1) as u32);
    }
    let sketch = SequencesSketch {
        kmer_counts: kmer_counts.clone(),
        c: 200,
        k: 31,
        file_name: "sample_1.fastq.gz".to_string(),
        sample_name: Some("my sample".to_string()),
        paired: true,
        mean_read_length: 149.7,
    };

    let mut buf = Vec::new();
    compress::write_seq_sketch_compressed(&mut buf, &sketch).unwrap();
    let decoded = compress::read_seq_sketch_compressed(&buf[..]).unwrap();

    assert_eq!(sketch, decoded, "SequencesSketch must roundtrip exactly");

    // also exercise the None sample_name / unpaired path
    let sketch2 = SequencesSketch {
        kmer_counts,
        c: 200,
        k: 21,
        file_name: "sample_2.fastq".to_string(),
        sample_name: None,
        paired: false,
        mean_read_length: 0.0,
    };
    let mut buf2 = Vec::new();
    compress::write_seq_sketch_compressed(&mut buf2, &sketch2).unwrap();
    assert_eq!(sketch2, compress::read_seq_sketch_compressed(&buf2[..]).unwrap());
}

#[test]
fn compress_detection_does_not_collide_with_legacy() {
    // A real compressed sketch is detected as compressed.
    let sketch = SequencesSketch {
        kmer_counts: FxHashMap::default(),
        c: 200,
        k: 31,
        file_name: "s.fq".to_string(),
        sample_name: None,
        paired: false,
        mean_read_length: 0.0,
    };
    let mut buf = Vec::new();
    compress::write_seq_sketch_compressed(&mut buf, &sketch).unwrap();
    let mut slice: &[u8] = &buf;
    assert!(compress::peek_is_compressed(&mut slice).unwrap());

    // A legacy bincode sketch whose leading length bytes are 1f 8b 08 (the old
    // gzip-style signature) must NOT be misdetected as compressed.
    let legacy_collision: &[u8] = &[0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00];
    let mut slice2 = legacy_collision;
    assert!(!compress::peek_is_compressed(&mut slice2).unwrap());
}

fn test_hash(){

    let key = 19238239812933123;
    println!("{}", format!("{key:b}"));
    let h = seeding::mm_hash64(key);
    println!("{}", format!("{h:b}"));
    let rev = seeding::rev_hash_64(h);
    println!("{}", format!("{rev:b}"));
    assert!(rev == key);

    if is_x86_feature_detected!("avx2"){
        unsafe{
            let key = key as i64;
            println!("{}", format!("{key:b}"));
            use std::arch::x86_64::*;
            use sylph::avx2_seeding::*;
            let mut rolling_kmer_f_marker = _mm256_set_epi64x(0, 0, 0, key);
                let hash_256 = mm_hash256(rolling_kmer_f_marker);
                let v1 = _mm256_extract_epi64(hash_256, 0);
                println!("{}", format!("{v1:b}"));
                assert!(v1 == h as i64);
        }

    }
}
