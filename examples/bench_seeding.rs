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
    let alphabet = *b"ACGT";
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

/// Time the kernels the way read sketching actually calls them: one call per
/// ~150 bp read, into a scratch `Vec` that is reused across reads. This is the
/// regime that matters for `sketch`/`profile` on fastq -- a 150 bp read gives the
/// AVX-512 path only `(150 - k + 1) / 8` iterations of its inner loop, so the
/// per-call priming of the eight lanes is a much larger share of the work than it
/// is on a whole genome.
fn bench_short_reads(seq: &[u8], read_len: usize, c: usize, k: usize) {
    let reads: Vec<&[u8]> = seq.chunks_exact(read_len).take(400_000).collect();
    let bp = reads.len() * read_len;
    println!(
        "short reads: {} x {} bp = {:.1} Mbp, k={k}, c={c}",
        reads.len(),
        read_len,
        bp as f64 / 1e6
    );

    let mut counts = Vec::new();
    for (label, backend) in [
        ("scalar", 0u8),
        ("avx2", 1),
        #[cfg(target_arch = "x86_64")]
        ("avx512", 2),
    ] {
        let mut scratch: Vec<u64> = Vec::new();
        let run = |scratch: &mut Vec<u64>| {
            let mut n = 0usize;
            for r in &reads {
                scratch.clear();
                match backend {
                    0 => fmh_seeds(r, scratch, c, k),
                    #[cfg(target_arch = "x86_64")]
                    1 => unsafe { extract_markers_avx2(r, scratch, c, k) },
                    #[cfg(target_arch = "x86_64")]
                    _ => unsafe { extract_markers_avx512(r, scratch, c, k) },
                    #[cfg(not(target_arch = "x86_64"))]
                    _ => fmh_seeds(r, scratch, c, k),
                }
                n += scratch.len();
            }
            n
        };
        let n = run(&mut scratch); // warm up
        let iters = 3;
        let t = Instant::now();
        for _ in 0..iters {
            run(&mut scratch);
        }
        let secs = t.elapsed().as_secs_f64();
        let ns_per_bp = secs * 1e9 / (bp as f64 * iters as f64);
        println!(
            "  {label:<8} {ns_per_bp:6.3} ns/bp   ({n} markers, {:.0} Mbp/s, {:.0} ns/read)",
            (bp as f64 * iters as f64) / secs / 1e6,
            secs * 1e9 / (reads.len() as f64 * iters as f64)
        );
        counts.push((label, n, ns_per_bp));
    }
    if let (Some(a), Some(v)) = (
        counts.iter().find(|x| x.0 == "avx2"),
        counts.iter().find(|x| x.0 == "avx512"),
    ) {
        assert_eq!(a.1, v.1, "AVX2 and AVX-512 marker counts differ");
        println!("  -> AVX-512 speedup vs AVX2: {:.2}x\n", a.2 / v.2);
    }
}

/// Time the `_positions` kernels, which genome sketching (not read sketching) uses:
/// they emit `(contig, end_position, hash)` tuples, so the AVX-512 path cannot use
/// `vpcompressq` to compact survivors -- it stores its eight hashes and loops over
/// them scalar-wise, exactly as AVX2 does over four -- while still paying twice the
/// per-lane priming.
fn bench_positions(seq: &[u8], c: usize, k: usize) {
    println!("positions kernels ({} bp contig, k={k}, c={c}):", seq.len());
    let mut timings = [0f64; 2];
    let mut counts = [0usize; 2];
    for (i, use512) in [false, true].iter().enumerate() {
        let mut v: Vec<(usize, usize, u64)> = Vec::new();
        let run = |v: &mut Vec<(usize, usize, u64)>| {
            v.clear();
            #[cfg(target_arch = "x86_64")]
            unsafe {
                if *use512 {
                    weebill::avx512_seeding::extract_markers_avx512_positions(seq, v, c, k, 0);
                } else {
                    weebill::avx2_seeding::extract_markers_avx2_positions(seq, v, c, k, 0);
                }
            }
        };
        run(&mut v);
        counts[i] = v.len();
        let iters = 3;
        let t = Instant::now();
        for _ in 0..iters {
            run(&mut v);
        }
        timings[i] = t.elapsed().as_secs_f64() * 1e9 / (seq.len() as f64 * iters as f64);
    }
    assert_eq!(counts[0], counts[1], "positions kernels disagree");
    println!(
        "  avx2   {:6.3} ns/bp\n  avx512 {:6.3} ns/bp\n  -> AVX-512 speedup vs AVX2: {:.2}x\n",
        timings[0],
        timings[1],
        timings[0] / timings[1]
    );
}

fn c_default() -> usize {
    200
}
fn k_default() -> usize {
    31
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
    println!("=== whole-genome regime (one call per contig) ===");
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

    println!("=== genome-sketching regime (the `_positions` kernels) ===");
    for &(k, c) in &[(31usize, 200usize), (31, 50)] {
        bench_positions(&seq[..5_000_000.min(seq.len())], c, k);
    }

    println!("=== crossover sweep: sequence length at which AVX-512 overtakes AVX2 ===");
    println!("  (one call per sequence of the given length, k=31, c=200)");
    for &len in &[
        100usize, 150, 200, 250, 300, 400, 500, 750, 1000, 1500, 2000, 3000, 5000, 10_000, 50_000,
    ] {
        let seqs: Vec<&[u8]> = seq
            .chunks_exact(len)
            .take(200_000.min(20_000_000 / len))
            .collect();
        let bp = seqs.len() * len;
        let mut scratch: Vec<u64> = Vec::new();
        let mut timings = [0f64; 2];
        for (i, use512) in [false, true].iter().enumerate() {
            let run = |scratch: &mut Vec<u64>| {
                for s in &seqs {
                    scratch.clear();
                    #[cfg(target_arch = "x86_64")]
                    unsafe {
                        if *use512 {
                            extract_markers_avx512(s, scratch, c_default(), k_default());
                        } else {
                            extract_markers_avx2(s, scratch, c_default(), k_default());
                        }
                    }
                }
            };
            run(&mut scratch);
            let iters = 3;
            let t = Instant::now();
            for _ in 0..iters {
                run(&mut scratch);
            }
            timings[i] = t.elapsed().as_secs_f64() * 1e9 / (bp as f64 * iters as f64);
        }
        println!(
            "  len {len:>6}: avx2 {:6.3} ns/bp   avx512 {:6.3} ns/bp   -> {:.2}x",
            timings[0],
            timings[1],
            timings[0] / timings[1]
        );
    }
    println!();

    println!("=== read regime (one call per read, reused scratch buffer) ===");
    for &(read_len, k, c) in &[
        (150usize, 31usize, 200usize),
        (150, 31, 100),
        (100, 31, 200),
        (250, 31, 200),
        (10_000, 31, 200), // long reads (nanopore-ish)
    ] {
        bench_short_reads(&seq, read_len, c, k);
    }
}
