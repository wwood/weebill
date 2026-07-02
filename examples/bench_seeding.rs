// Correctness + throughput benchmark for the k-mer marker extraction kernels.
//
//   cargo run --release --example bench_seeding
//
// Compares the scalar (`fmh_seeds`), AVX2 (`extract_markers_avx2`) and new
// AVX-512 (`extract_markers_avx512`) paths. Verifies the AVX-512 path emits the
// exact same multiset of hashes as the AVX2 path (sketch compatibility), then
// times all three on a large real genome sequence.

use std::time::Instant;
use weebill::seeding::fmh_seeds;
use weebill::types::BYTE_TO_SEQ;

#[cfg(target_arch = "x86_64")]
use weebill::avx2_seeding::extract_markers_avx2;
#[cfg(target_arch = "x86_64")]
use weebill::avx512_seeding::extract_markers_avx512;

fn read_fasta_gz(path: &str, out: &mut Vec<u8>) {
    let mut reader = needletail::parse_fastx_file(path).expect("open fasta");
    while let Some(r) = reader.next() {
        let rec = r.expect("record");
        for &b in rec.seq().iter() {
            // keep only ACGT-ish bytes the sketcher understands
            let up = b.to_ascii_uppercase();
            if matches!(up, b'A' | b'C' | b'G' | b'T') {
                out.push(up);
            }
        }
    }
}

fn avx2(seq: &[u8], c: usize, k: usize) -> Vec<u64> {
    let mut v = Vec::new();
    #[cfg(target_arch = "x86_64")]
    unsafe {
        extract_markers_avx2(seq, &mut v, c, k);
    }
    v
}

fn avx512(seq: &[u8], c: usize, k: usize) -> Vec<u64> {
    let mut v = Vec::new();
    #[cfg(target_arch = "x86_64")]
    unsafe {
        extract_markers_avx512(seq, &mut v, c, k);
    }
    v
}

fn scalar(seq: &[u8], c: usize, k: usize) -> Vec<u64> {
    let mut v = Vec::new();
    fmh_seeds(seq, &mut v, c, k);
    v
}

fn sorted(mut v: Vec<u64>) -> Vec<u64> {
    v.sort_unstable();
    v
}

fn fuzz_correctness() {
    // Deterministic pseudo-random sequences of many lengths; AVX-512 must
    // reproduce the AVX2 multiset exactly for every length/k/c.
    let mut state: u64 = 0x9E3779B97F4A7C15;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    let alphabet = [b'A', b'C', b'G', b'T'];
    let mut failures = 0;
    for len in 0..600usize {
        let mut seq = Vec::with_capacity(len);
        for _ in 0..len {
            seq.push(alphabet[(next() & 3) as usize]);
        }
        for &k in &[21usize, 31] {
            for &c in &[1usize, 3, 200] {
                let a = sorted(avx2(&seq, c, k));
                let b = sorted(avx512(&seq, c, k));
                if a != b {
                    failures += 1;
                    if failures <= 5 {
                        eprintln!(
                            "MISMATCH len={len} k={k} c={c}: avx2={} avx512={}",
                            a.len(),
                            b.len()
                        );
                    }
                }
            }
        }
    }
    if failures == 0 {
        println!("fuzz: AVX-512 multiset == AVX2 multiset for all 600 lengths x {{k=21,31}} x {{c=1,3,200}}  ✓");
    } else {
        println!("fuzz: {failures} MISMATCHES  ✗");
        std::process::exit(1);
    }
}

fn time<F: Fn() -> Vec<u64>>(label: &str, bp: usize, iters: usize, f: F) -> (usize, f64) {
    // warm up
    let n = f().len();
    let t = Instant::now();
    let mut total = 0usize;
    for _ in 0..iters {
        total += f().len();
    }
    let secs = t.elapsed().as_secs_f64();
    let ns_per_bp = secs * 1e9 / (bp as f64 * iters as f64);
    println!(
        "  {label:<8} {ns_per_bp:6.3} ns/bp   ({n} markers, {:.0} Mbp/s)",
        (bp as f64 * iters as f64) / secs / 1e6
    );
    let _ = total;
    (n, ns_per_bp)
}

fn main() {
    println!(
        "AVX2 detected: {}, AVX-512F detected: {}\n",
        is_x86_feature_detected!("avx2"),
        is_x86_feature_detected!("avx512f")
    );

    fuzz_correctness();

    // Build a large realistic sequence (~ tens of Mbp) from the bundled genomes.
    let mut seq = Vec::new();
    for f in [
        "test_files/e.coli-K12.fasta.gz",
        "test_files/e.coli-EC590.fasta.gz",
        "test_files/e.coli-o157.fasta.gz",
    ] {
        read_fasta_gz(f, &mut seq);
    }
    // Repeat to get a stable, cache-cold-ish working set (~40 Mbp).
    let base_len = seq.len();
    while seq.len() < 40_000_000 {
        seq.extend_from_within(0..base_len);
    }
    // sanity: all bytes are ACGT
    debug_assert!(seq.iter().all(|&b| BYTE_TO_SEQ[b as usize] < 4));

    println!("\nBenchmark sequence: {} bp\n", seq.len());
    let iters = 3;
    for &(k, c) in &[(31usize, 200usize), (21, 200), (31, 3), (31, 50)] {
        println!("k={k}, c={c}:");
        let (_, s) = time("scalar", seq.len(), iters, || scalar(&seq, c, k));
        let (_, a) = time("avx2", seq.len(), iters, || avx2(&seq, c, k));
        let (_, v) = time("avx512", seq.len(), iters, || avx512(&seq, c, k));
        println!(
            "  -> AVX-512 speedup vs AVX2: {:.2}x , vs scalar: {:.2}x\n",
            a / v,
            s / v
        );
    }
}
