use assert_cmd::prelude::*; // Add methods on commands
use sylph::seeding;
use sylph::compress;
use sylph::types::{GenomeSketch, SequencesSketch};
use fxhash::FxHashMap;

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
    assert_eq!(buf[0], 0x1f, "compressed output must start with gzip magic");
    assert_eq!(buf[1], 0x8b);
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
