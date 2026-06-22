use assert_cmd::prelude::*; // Add methods on commands
use fxhash::FxHashMap;
use sylph::compress;
use sylph::refdelta;
use sylph::seeding;
use sylph::types::{GenomeSketch, SequencesSketch};

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

    let db = refdelta::build_refdb_with_pool_min_genomes(&sketches, &tax, 2);
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

#[test]
fn refdelta_pool_min_genomes_assigns_pairs_to_first_owner() {
    let sketches = vec![
        gsketch("A_rep.fa", vec![1, 10, 20]),
        gsketch("B_rep.fa", vec![2, 10, 20]),
        gsketch("C_rep.fa", vec![3, 20]),
    ];
    let mut tax: FxHashMap<String, (String, bool)> = FxHashMap::default();
    tax.insert("A_rep.fa".into(), ("A".into(), true));
    tax.insert("B_rep.fa".into(), ("B".into(), true));
    tax.insert("C_rep.fa".into(), ("C".into(), true));

    let db = refdelta::build_refdb_with_pool_min_genomes(&sketches, &tax, 3);

    // 10 is shared by exactly two reps, so threshold 3 keeps it distinctive and
    // assigns it to the first genome in build order. 20 is shared by three reps
    // and still goes to the shared pool.
    assert_eq!(db.distinctive[0], vec![1, 10]);
    assert_eq!(db.distinctive[1], vec![2]);
    assert_eq!(db.distinctive[2], vec![3]);
    assert_eq!(db.pool, vec![20]);
}

fn open_index(db: &refdelta::RefDb, sparse_div: u64) -> refdelta::RefIndex {
    let mut dbuf = Vec::new();
    refdelta::write_refdb(&mut dbuf, db, sparse_div).unwrap();
    refdelta::open_ref_index(std::io::Cursor::new(dbuf)).unwrap()
}

fn refdelta_roundtrip(sketch: &SequencesSketch, idx: &refdelta::RefIndex) {
    let mut buf = Vec::new();
    refdelta::compress_seq(&mut buf, sketch, idx, "unit-test.sylref").unwrap();
    let decoded = refdelta::decompress_seq(&buf[..], idx).unwrap();
    assert_eq!(*sketch, decoded);
}

fn three_genome_db() -> refdelta::RefDb {
    let sketches = vec![
        gsketch("A_rep.fa", vec![1, 2, 3, 100, 101, 500]),
        gsketch("A_str.fa", vec![100, 101, 200, 201]),
        gsketch("B_rep.fa", vec![300, 301, 500]),
    ];
    let mut tax: FxHashMap<String, (String, bool)> = FxHashMap::default();
    tax.insert("A_rep.fa".into(), ("A".into(), true));
    tax.insert("A_str.fa".into(), ("A".into(), false));
    tax.insert("B_rep.fa".into(), ("B".into(), true));
    refdelta::build_refdb(&sketches, &tax)
}

#[test]
fn refdelta_compress_decompress_roundtrip() {
    let db = three_genome_db();
    // sparse_div = 1 keeps every distinctive k-mer in the stage-1 index so every
    // genome is detectable; the round trip is lossless regardless either way.
    let idx = open_index(&db, 1);

    // distinctive hits + pool hit + novel hashes + counts
    let mut counts: FxHashMap<u64, u32> = FxHashMap::default();
    for (h, c) in [
        (1u64, 5u32),
        (2, 3),
        (100, 7),
        (200, 2),
        (500, 4),
        (999, 1),
        (98765, 9),
    ] {
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
    refdelta_roundtrip(&sketch, &idx);

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
        &idx,
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
        &idx,
    );

    // round trip is also lossless when stage 1 subsamples (some genomes may be
    // missed and fall back to novel coding, but the result must be identical)
    let idx_sparse = open_index(&db, 4);
    let mut counts2: FxHashMap<u64, u32> = FxHashMap::default();
    for (h, c) in [(1u64, 5u32), (100, 7), (300, 2), (500, 4), (424242, 9)] {
        counts2.insert(h, c);
    }
    refdelta_roundtrip(
        &SequencesSketch {
            kmer_counts: counts2,
            c: 200,
            k: 31,
            file_name: "s2.fq".into(),
            sample_name: None,
            paired: false,
            mean_read_length: 100.0,
        },
        &idx_sparse,
    );
}

#[test]
fn refdelta_sparse_hit_detection() {
    let db = three_genome_db();
    let idx = open_index(&db, 1);
    // a sample containing only B_rep's distinctive k-mer (300) should detect
    // exactly genome 2 (build order: A_rep, A_str, B_rep).
    let mut counts: FxHashMap<u64, u32> = FxHashMap::default();
    counts.insert(300, 5);
    let sketch = SequencesSketch {
        kmer_counts: counts,
        c: 200,
        k: 31,
        file_name: "x.fq".into(),
        sample_name: None,
        paired: false,
        mean_read_length: 1.0,
    };
    assert_eq!(idx.hit_genomes(&sketch), vec![2]);
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
    assert_eq!(
        &buf[..4],
        b"SYLZ",
        "compressed output must start with SYLZ magic"
    );
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
        assert_eq!(
            a, dec.genome_kmers,
            "genome_kmers must match as a sorted multiset"
        );

        match (
            &orig.pseudotax_tracked_nonused_kmers,
            &dec.pseudotax_tracked_nonused_kmers,
        ) {
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
    assert_eq!(
        sketch2,
        compress::read_seq_sketch_compressed(&buf2[..]).unwrap()
    );
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

fn test_hash() {
    let key = 19238239812933123;
    println!("{}", format!("{key:b}"));
    let h = seeding::mm_hash64(key);
    println!("{}", format!("{h:b}"));
    let rev = seeding::rev_hash_64(h);
    println!("{}", format!("{rev:b}"));
    assert!(rev == key);

    if is_x86_feature_detected!("avx2") {
        unsafe {
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
