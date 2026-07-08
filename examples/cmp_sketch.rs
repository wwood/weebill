// Compare two paired .sylsp (plain bincode SequencesSketch) as (kmer -> count) maps.
use std::fs::File;
use std::io::BufReader;
use weebill::types::SequencesSketch;

fn load(path: &str) -> SequencesSketch {
    let f = BufReader::new(File::open(path).unwrap());
    bincode::deserialize_from(f).unwrap()
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let a = load(&args[1]);
    let b = load(&args[2]);
    println!(
        "A: {} kmers, sum={}, mean_read_len={}",
        a.kmer_counts.len(),
        a.kmer_counts.values().map(|v| *v as u64).sum::<u64>(),
        a.mean_read_length
    );
    println!(
        "B: {} kmers, sum={}, mean_read_len={}",
        b.kmer_counts.len(),
        b.kmer_counts.values().map(|v| *v as u64).sum::<u64>(),
        b.mean_read_length
    );
    let maps_equal = a.kmer_counts == b.kmer_counts;
    println!("kmer_counts maps equal (as sets): {}", maps_equal);
    if !maps_equal {
        let mut only_a = 0usize;
        let mut diff_count = 0usize;
        for (k, va) in a.kmer_counts.iter() {
            match b.kmer_counts.get(k) {
                None => only_a += 1,
                Some(vb) if vb != va => diff_count += 1,
                _ => {}
            }
        }
        let only_b = b
            .kmer_counts
            .keys()
            .filter(|k| !a.kmer_counts.contains_key(k))
            .count();
        println!(
            "keys only in A: {}, only in B: {}, shared-but-count-differs: {}",
            only_a, only_b, diff_count
        );
    }
    println!(
        "meta equal (c,k,paired,mean): {}",
        a.c == b.c
            && a.k == b.k
            && a.paired == b.paired
            && a.mean_read_length == b.mean_read_length
    );
}
