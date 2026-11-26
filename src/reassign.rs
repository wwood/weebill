//! Coverage-aware (EM) reassignment of k-mers shared between candidate genomes.
//!
//! sylph's `profile` gives every k-mer shared by several candidate genomes, with its
//! whole read count, to the candidate with the highest ANI ("winner-take-all"). When
//! two related genomes are both in the sample, a shared k-mer's count is really about
//! `λ_X + λ_Y`, so winner-take-all inflates the winner's coverage and strips the
//! loser of k-mers it does contain, dragging its ANI and coverage down (sometimes
//! below the profiling thresholds).
//!
//! The model here: each candidate genome `g` has a k-mer coverage `λ_g` and a
//! presence rate `π_g` (the fraction of its k-mers present in the sampled strain,
//! roughly ANI^k). For a genome k-mer `h` with sample count `c_h`, let `G_h` be the
//! candidates that carry it and `z_g ∈ {0,1}` say whether `g`'s copy is in the sample
//! (`P(z_g = 1) = π_g`). Then `c_h ~ Poisson(Σ_{g ∈ G_h} z_g λ_g)`. For a k-mer
//! carried by one genome this is exactly a zero-inflated Poisson.
//!
//! EM over the latent `z` (enumerating subsets of `G_h`) estimates `(π_g, λ_g)`.
//! The E-step gives each genome an expected histogram of its per-k-mer counts
//! (unique k-mers at their count, shared k-mers spread over the posterior of their
//! share). The M-step reads λ from that histogram with sylph's modal-ratio
//! estimator rather than the Poisson mean, so a thin tail of k-mers whose counts
//! really come from a co-present relative cannot inflate the λ of a low-coverage
//! genome (at λ ≈ 0.04, a few dozen such k-mers would otherwise give λ ≈ 0.3).
//! Above a median count of 2 it switches, as sylph does, to the mean (or median)
//! of the contained counts (see `fit_params`). Counts above a genome's outlier
//! cap are left out of the fit; the caps start from the pass-1 medians and are
//! then refitted from the EM's own histograms until they settle.
//! It starts from the winner-take-all assignment: each genome's initial parameters
//! come from the k-mers it wins, with lost k-mers treated as missing rather than as
//! zeros. That start doubles as the tie-break the model cannot make by itself: a
//! genome whose k-mers are all shared with a better genome (a redundant near-copy
//! that is not really in the sample) starts with `π ≈ 0` and stays there.
//!
//! After convergence, each shared k-mer's count is split among the genomes into
//! integer per-genome counts: a presence pattern `z` is drawn from its posterior, and
//! the count is divided multinomially in proportion to `λ` across the present
//! genomes. The draw is seeded from the k-mer hash, so results are deterministic and
//! independent of thread count (given the candidate order, which the caller fixes). The per-genome count vectors then go through the
//! usual coverage/ANI estimation, just like the winner-take-all ones.

use crate::constants::*;
use crate::types::*;
use fxhash::FxHashMap;
use rayon::prelude::*;
use smallvec::SmallVec;
use statrs::distribution::{DiscreteCDF, Poisson};

/// Maximum EM iterations.
const EM_MAX_ITER: usize = 200;
/// Convergence: stop when no parameter moves by more than this (relative).
const EM_TOL: f64 = 1e-6;
/// A shared k-mer with more carriers than this only considers the most abundant
/// (`π·λ`) ones as possible sources (the posterior enumerates `2^n` subsets).
const EM_MAX_EXACT_MEMBERS: usize = 10;
const EM_MIN_PI: f64 = 1e-6;
const EM_MAX_PI: f64 = 1. - 1e-9;
const EM_MIN_LAMBDA: f64 = 1e-6;
/// Maximum fits while the outlier caps settle (see `em_reassign`).
const EM_MAX_CAP_ROUNDS: usize = 5;
/// Member ids with this bit set carry the k-mer only in
/// `pseudotax_tracked_nonused_kmers`: they can absorb part of its count (as they claim
/// it under winner-take-all) but it is not one of their counted k-mers.
const TRACKED_BIT: u32 = 1 << 31;

/// Pass-1 (independent, unreassigned) estimates for one candidate genome.
pub struct EmCandidate<'a> {
    pub genome: &'a GenomeSketch,
    /// Pass-1 ANI, used only for the winner-take-all starting point.
    pub ani: f64,
    /// Pass-1 effective coverage; fallback λ for a genome that wins nothing.
    pub cov: f64,
    /// Pass-1 median count; sets the outlier cap as in `finalize_stats`.
    pub median_cov: f64,
}

/// Per-genome outcome of the EM reassignment.
pub struct EmAssignment {
    /// Non-zero per-genome counts of the genome's counted k-mers (unique k-mers
    /// keep their full count; shared ones get their share). Ready for
    /// `finalize_stats` as the contained-k-mer coverage vector.
    pub covs: Vec<u32>,
    pub lambda: f64,
    pub pi: f64,
}

struct SharedKmer {
    kmer: Kmer,
    count: u32,
    members: SmallVec<[u32; 4]>,
    /// Count above every carrier's outlier cap (repeat / contaminant): kept out of
    /// the parameter estimates, but still split at the end.
    outlier: bool,
}

#[derive(Clone, Copy)]
struct Params {
    pi: f64,
    lambda: f64,
}

/// Largest count `finalize_stats` keeps for a genome with this median count:
/// the counts whose Poisson(median) CDF stays below `CUTOFF_PVALUE`.
fn outlier_cap(median_cov: f64) -> u32 {
    if median_cov.is_nan() || median_cov >= 30. || median_cov <= 0. {
        return u32::MAX;
    }
    let pois = Poisson::new(median_cov).unwrap();
    let mut c = median_cov as u64;
    while pois.cdf(c) < CUTOFF_PVALUE {
        c += 1;
    }
    // `c` is the first count at/after the median whose CDF reaches the cutoff;
    // finalize_stats keeps everything strictly before it (and at least the median).
    (c.saturating_sub(1).max(median_cov as u64)) as u32
}

/// λ of a zero-truncated Poisson with the given mean of the non-zero counts:
/// solves `λ / (1 - e^-λ) = m`.
fn ztp_lambda(m: f64) -> f64 {
    if m.is_nan() || m <= 1. + 1e-9 {
        return EM_MIN_LAMBDA.max(1e-3);
    }
    let mut lam = m;
    for _ in 0..100 {
        let e = (-lam).exp();
        let f = lam - m * (1. - e);
        let df = 1. - m * e;
        if df.abs() < 1e-12 {
            break;
        }
        let next = lam - f / df;
        if !next.is_finite() || next <= 0. {
            lam /= 2.;
            continue;
        }
        if (next - lam).abs() < 1e-12 {
            lam = next;
            break;
        }
        lam = next;
    }
    lam.max(EM_MIN_LAMBDA)
}

/// Fixed-size E-step chunks (reduced in order: thread-count independent sums).
const EM_CHUNK: usize = 4096;
/// Minimum expected k-mers at the modal count and the one above it for the
/// ratio estimate of λ (sylph's `--min-count-correct` default).
const EM_MIN_RATIO_COUNT: f64 = 3.;

fn hist_add(h: &mut Vec<f64>, j: usize, w: f64) {
    if h.len() <= j {
        h.resize(j + 1, 0.);
    }
    h[j] += w;
}

/// Add `w` × Binomial(c, q) over shares j >= 1 to `h`.
fn add_binomial(h: &mut Vec<f64>, c: u32, q: f64, w: f64) {
    if q >= 1. - 1e-12 {
        hist_add(h, c as usize, w);
        return;
    }
    if q <= 0. || c == 0 {
        return;
    }
    if h.len() <= c as usize {
        h.resize(c as usize + 1, 0.);
    }
    let r = q / (1. - q);
    let mut pmf = (1. - q).powi(c as i32);
    for j in 0..c {
        pmf *= (c - j) as f64 / (j + 1) as f64 * r;
        h[j as usize + 1] += w * pmf;
    }
}

/// Median of a genome's contained counts from its (expected) histogram of
/// non-outlier counts, with `n_outlier` further contained k-mers (all above it)
/// included, picked like `finalize_stats`'s `covs[covs.len() / 2]`.
fn hist_median(h: &[f64], n_outlier: f64) -> Option<f64> {
    let nonzero: f64 = h.iter().skip(1).sum();
    if nonzero <= 0. {
        return None;
    }
    let target = ((nonzero + n_outlier) / 2.).floor();
    let mut acc = 0.;
    for (j, w) in h.iter().enumerate().skip(1) {
        acc += w;
        if acc > target + 1e-9 {
            return Some(j as f64);
        }
    }
    Some(h.len() as f64)
}

/// (π, λ) from an expected histogram of a genome's non-zero, non-outlier counts
/// over `n_eff` k-mers, with `n_outlier` further contained k-mers whose counts were
/// set aside as outliers. λ is the Poisson rate behind the estimator
/// `finalize_stats` uses for the reported coverage, so that the EM splits shared
/// k-mers by the coverages sylph would report:
///
/// * median count `<= MEDIAN_ANI_THRESHOLD`: the ratio of the modal count and the
///   one above it, as in sylph's default estimator (`ratio_lambda`): for Poisson,
///   `P(m+1)/P(m) = λ/(m+1)`. This only looks at the bulk of the distribution, so a
///   thin tail of k-mers whose counts really come from another genome (e.g. ones
///   that the other genome's sketch dropped as duplicates) cannot inflate it the way
///   it inflates a mean -- at low coverage that difference decides the fit. Falls
///   back to the zero-truncated-Poisson mean when the ratio is not supported.
/// * higher median: sylph does not fit a Poisson there (`AdjustStatus::High`) but
///   reports the mean of the contained k-mers' counts, or their median once that
///   reaches `MAX_MEDIAN_FOR_MEAN_FINAL_EST` (unless `mean_coverage`). The modal
///   ratio would be biased low there: neighbouring bins are of similar height, and
///   the mode is the one that happens to be the largest. The mean is over non-zero
///   counts only, so λ is the rate whose zero-truncated mean it is
///   (`λ / (1 - e^-λ) = mean`). Using the mean itself as the rate inflates a
///   genome that shares most k-mers with a more abundant relative several-fold:
///   its share of a shared count is about the count minus the relative's λ, so the
///   2% excess at λ ≈ 4 becomes ~8% (20x K-12 + 5x EC590: 4.04 vs 3.80).
///
/// π is the non-zero fraction over its Poisson expectation `1 - e^-λ` in both
/// cases, not the plain containment sylph reports at high coverage: the model's
/// shared k-mers are present with probability π and then non-zero with probability
/// `1 - e^-λ`, and with the plain fraction each iteration shrinks π by that factor
/// wherever the counts alone cannot tell whether the genome contributes (20x K-12
/// + 5x EC590: π 0.90, ~2000 EC590 k-mers lost).
fn fit_params(h: &[f64], n_eff: f64, n_outlier: f64, mean_coverage: bool) -> Option<Params> {
    let nonzero: f64 = h.iter().skip(1).sum();
    if nonzero <= 0. || n_eff <= 0. {
        return None;
    }
    let sum: f64 = h.iter().enumerate().map(|(j, w)| j as f64 * w).sum();
    let contained = nonzero + n_outlier;
    let median = hist_median(h, n_outlier).unwrap();
    if median > MEDIAN_ANI_THRESHOLD {
        let lambda = if median < MAX_MEDIAN_FOR_MEAN_FINAL_EST || mean_coverage {
            ztp_lambda(sum / contained)
        } else {
            median
        }
        .max(EM_MIN_LAMBDA);
        let pi = (nonzero / n_eff / (1. - (-lambda).exp())).clamp(EM_MIN_PI, EM_MAX_PI);
        return Some(Params { pi, lambda });
    }
    let mut mode = 1;
    for j in 2..h.len() {
        if h[j] > h[mode] {
            mode = j;
        }
    }
    let next = h.get(mode + 1).copied().unwrap_or(0.);
    let lambda = if h[mode] >= EM_MIN_RATIO_COUNT && next >= EM_MIN_RATIO_COUNT {
        (mode + 1) as f64 * next / h[mode]
    } else {
        ztp_lambda(sum / nonzero)
    }
    .max(EM_MIN_LAMBDA);
    let pi = (nonzero / n_eff / (1. - (-lambda).exp())).clamp(EM_MIN_PI, EM_MAX_PI);
    Some(Params { pi, lambda })
}

/// The carriers that may have produced a shared k-mer's count, as indices into
/// `members`: all of them, or the `EM_MAX_EXACT_MEMBERS` with the largest `π·λ`.
fn active_members(members: &[u32], params: &[Params]) -> SmallVec<[usize; 8]> {
    let mut act: SmallVec<[usize; 8]> = (0..members.len()).collect();
    if act.len() > EM_MAX_EXACT_MEMBERS {
        act.sort_by(|&a, &b| {
            let pa = &params[(members[a] & !TRACKED_BIT) as usize];
            let pb = &params[(members[b] & !TRACKED_BIT) as usize];
            (pb.pi * pb.lambda)
                .partial_cmp(&(pa.pi * pa.lambda))
                .unwrap()
                .then(a.cmp(&b))
        });
        act.truncate(EM_MAX_EXACT_MEMBERS);
        act.sort();
    }
    act
}

/// Posterior over presence patterns of a shared k-mer's active carriers, as
/// `(mask, weight)` pairs (bit `i` of `mask` = `act[i]` present), weights summing to 1.
fn posterior(count: u32, members: &[u32], act: &[usize], params: &[Params]) -> Vec<(u32, f64)> {
    let s = act.len();
    let c = count as f64;
    let ps: SmallVec<[Params; 8]> = act
        .iter()
        .map(|&i| params[(members[i] & !TRACKED_BIT) as usize])
        .collect();
    let mut out = Vec::with_capacity((1usize << s) - 1);
    let mut max_lw = f64::NEG_INFINITY;
    for mask in 1u32..(1u32 << s) {
        let mut lw = 0.;
        let mut lam = 0.;
        for (i, p) in ps.iter().enumerate() {
            if mask & (1 << i) != 0 {
                lw += p.pi.ln();
                lam += p.lambda;
            } else {
                lw += (1. - p.pi).ln();
            }
        }
        // Poisson(c; lam) up to the constant -ln(c!).
        lw += c * lam.ln() - lam;
        max_lw = max_lw.max(lw);
        out.push((mask, lw));
    }
    let mut tot = 0.;
    for (_, w) in out.iter_mut() {
        *w = (*w - max_lw).exp();
        tot += *w;
    }
    for (_, w) in out.iter_mut() {
        *w /= tot;
    }
    out
}

/// Winner-take-all start: the first carrier with the highest pass-1 ANI wins
/// (the same rule as `winner_table`); lost k-mers are missing, not zeros.
fn wta_start(
    cands: &[EmCandidate],
    shared: &[SharedKmer],
    uniq_hist: &[Vec<f64>],
    n_kmers: &[f64],
    n_outlier: &[f64],
    mean_coverage: bool,
) -> Vec<Params> {
    let n = cands.len();
    let mut init_hist = uniq_hist.to_vec();
    let mut lost_n = vec![0f64; n];
    for s in shared.iter().filter(|s| !s.outlier) {
        let mut win = s.members[0];
        for &m in s.members.iter().skip(1) {
            if cands[(m & !TRACKED_BIT) as usize].ani > cands[(win & !TRACKED_BIT) as usize].ani {
                win = m;
            }
        }
        for &m in s.members.iter().filter(|&&m| m & TRACKED_BIT == 0) {
            if m == win {
                hist_add(&mut init_hist[m as usize], s.count as usize, 1.);
            } else {
                lost_n[m as usize] += 1.;
            }
        }
    }
    (0..n)
        .map(|g| {
            fit_params(
                &init_hist[g],
                n_kmers[g] - lost_n[g] - n_outlier[g],
                n_outlier[g],
                mean_coverage,
            )
            .unwrap_or(Params {
                pi: EM_MIN_PI,
                lambda: cands[g].cov.max(1e-3),
            })
        })
        .collect()
}

/// EM iterations over the non-outlier shared k-mers, updating `params` in place.
/// Returns the expected count histograms of the last E-step.
fn em_iterate(
    shared: &[SharedKmer],
    uniq_hist: &[Vec<f64>],
    n_kmers: &[f64],
    n_outlier: &[f64],
    mean_coverage: bool,
    params: &mut [Params],
) -> Vec<Vec<f64>> {
    let n = params.len();
    let fit: Vec<&SharedKmer> = shared.iter().filter(|s| !s.outlier).collect();
    let mut hist = uniq_hist.to_vec();
    for iter in 0..EM_MAX_ITER {
        // E-step: each shared k-mer's posterior distribution of every carrier's
        // share, accumulated into expected count histograms. Chunks are fixed-size
        // and reduced in order, so the sums do not depend on the thread count.
        let partial: Vec<FxHashMap<u32, Vec<f64>>> = fit
            .par_chunks(EM_CHUNK)
            .map(|chunk| {
                let mut local: FxHashMap<u32, Vec<f64>> = FxHashMap::default();
                for s in chunk {
                    let act = active_members(&s.members, params);
                    let post = posterior(s.count, &s.members, &act, params);
                    for (i, &ai) in act.iter().enumerate() {
                        let m = s.members[ai];
                        if m & TRACKED_BIT != 0 {
                            continue;
                        }
                        let lam_i = params[m as usize].lambda;
                        let h = local.entry(m).or_default();
                        for &(mask, w) in post.iter() {
                            if mask & (1 << i) == 0 {
                                continue;
                            }
                            let tot: f64 = act
                                .iter()
                                .enumerate()
                                .filter(|(j, _)| mask & (1 << j) != 0)
                                .map(|(_, &aj)| {
                                    params[(s.members[aj] & !TRACKED_BIT) as usize].lambda
                                })
                                .sum();
                            add_binomial(h, s.count, lam_i / tot, w);
                        }
                    }
                }
                local
            })
            .collect();
        hist = uniq_hist.to_vec();
        for local in partial.iter() {
            let mut keys: Vec<&u32> = local.keys().collect();
            keys.sort();
            for m in keys {
                let h = &mut hist[*m as usize];
                for (j, w) in local[m].iter().enumerate() {
                    if *w > 0. {
                        hist_add(h, j, *w);
                    }
                }
            }
        }

        // M-step.
        let mut max_delta: f64 = 0.;
        for g in 0..n {
            let Params { pi, lambda } = params[g];
            let Some(new) = fit_params(
                &hist[g],
                n_kmers[g] - n_outlier[g],
                n_outlier[g],
                mean_coverage,
            ) else {
                continue;
            };
            max_delta = max_delta
                .max((new.pi - pi).abs() / pi.max(1e-12))
                .max((new.lambda - lambda).abs() / lambda.max(1e-12));
            params[g] = new;
        }
        if max_delta < EM_TOL {
            log::debug!("EM reassignment converged after {} iterations", iter + 1);
            break;
        }
    }
    hist
}

/// Run the EM reassignment over the pass-1 candidates of one sample.
/// `mean_coverage` mirrors `profile --mean-coverage` (see `fit_params`).
pub fn em_reassign(
    cands: &[EmCandidate],
    kmer_counts: &FxHashMap<Kmer, u32>,
    mean_coverage: bool,
) -> Vec<EmAssignment> {
    let n = cands.len();
    assert!(n < TRACKED_BIT as usize);

    // Observed k-mer -> carrying candidates.
    let mut carriers: FxHashMap<Kmer, SmallVec<[u32; 4]>> = FxHashMap::default();
    let mut n_obs = vec![0usize; n];
    for (g, cand) in cands.iter().enumerate() {
        for kmer in cand.genome.genome_kmers.iter() {
            if kmer_counts.get(kmer).copied().unwrap_or(0) > 0 {
                n_obs[g] += 1;
                carriers.entry(*kmer).or_default().push(g as u32);
            }
        }
        if let Some(tracked) = &cand.genome.pseudotax_tracked_nonused_kmers {
            for kmer in tracked.iter() {
                if kmer_counts.get(kmer).copied().unwrap_or(0) > 0 {
                    let m = carriers.entry(*kmer).or_default();
                    if !m.iter().any(|&x| x & !TRACKED_BIT == g as u32) {
                        m.push(g as u32 | TRACKED_BIT);
                    }
                }
            }
        }
    }

    // Split into unique (fixed contribution) and shared k-mers.
    let mut uniq_covs: Vec<Vec<u32>> = vec![vec![]; n];
    let mut shared: Vec<SharedKmer> = vec![];
    for (kmer, members) in carriers.into_iter() {
        let count = kmer_counts[&kmer];
        if members.len() == 1 {
            if members[0] & TRACKED_BIT == 0 {
                uniq_covs[members[0] as usize].push(count);
            }
            continue;
        }
        shared.push(SharedKmer {
            kmer,
            count,
            outlier: false,
            members,
        });
    }
    shared.sort_unstable_by_key(|s| s.kmer);
    let n_kmers: Vec<f64> = cands
        .iter()
        .map(|c| c.genome.genome_kmers.len() as f64)
        .collect();

    // Outlier caps start from the pass-1 medians, but a genome sharing most of its
    // k-mers with an abundant relative has a pass-1 median of roughly the sum of
    // both coverages, so its cap is far too loose. Each fit is therefore followed
    // by caps from the genome's own fitted count histogram (the median that
    // `finalize_stats` will see after the split), and refitted until they settle.
    let mut caps: Vec<u32> = cands.iter().map(|c| outlier_cap(c.median_cov)).collect();
    let mut params: Option<Vec<Params>> = None;
    for cap_round in 0..EM_MAX_CAP_ROUNDS {
        for s in shared.iter_mut() {
            let cap = s
                .members
                .iter()
                .map(|&m| caps[(m & !TRACKED_BIT) as usize] as u64)
                .sum::<u64>();
            s.outlier = s.count as u64 > cap;
        }
        // Fixed per-genome statistics: histograms (index = count) of the unique,
        // non-outlier counts.
        let mut uniq_hist: Vec<Vec<f64>> = vec![vec![]; n];
        let mut n_outlier = vec![0f64; n];
        for g in 0..n {
            for &c in uniq_covs[g].iter() {
                if c <= caps[g] {
                    hist_add(&mut uniq_hist[g], c as usize, 1.);
                } else {
                    n_outlier[g] += 1.;
                }
            }
        }
        for s in shared.iter().filter(|s| s.outlier) {
            for &m in s.members.iter().filter(|&&m| m & TRACKED_BIT == 0) {
                n_outlier[m as usize] += 1.;
            }
        }
        let mut p = match params.take() {
            Some(p) => p,
            None => wta_start(
                cands,
                &shared,
                &uniq_hist,
                &n_kmers,
                &n_outlier,
                mean_coverage,
            ),
        };
        let hist = em_iterate(
            &shared,
            &uniq_hist,
            &n_kmers,
            &n_outlier,
            mean_coverage,
            &mut p,
        );
        params = Some(p);
        let new_caps: Vec<u32> = (0..n)
            .map(|g| match hist_median(&hist[g], n_outlier[g]) {
                Some(median) => outlier_cap(median),
                None => caps[g],
            })
            .collect();
        if new_caps == caps {
            break;
        }
        log::debug!(
            "EM outlier caps changed after fit {}; refitting",
            cap_round + 1
        );
        caps = new_caps;
    }
    let params = params.unwrap();

    // Split each shared k-mer's count into integer per-genome shares.
    let draws: Vec<SmallVec<[(u32, u32); 4]>> = shared
        .par_iter()
        .map(|s| {
            let act = active_members(&s.members, &params);
            let post = posterior(s.count, &s.members, &act, &params);
            let mut rng = fastrand::Rng::with_seed(
                s.kmer.wrapping_mul(0x9E37_79B9_7F4A_7C15) ^ DEFAULT_RNG_SEED,
            );
            let u = rng.f64();
            let mut acc = 0.;
            let mut mask = post.last().unwrap().0;
            for &(m, w) in post.iter() {
                acc += w;
                if u < acc {
                    mask = m;
                    break;
                }
            }
            let present: SmallVec<[usize; 8]> = act
                .iter()
                .enumerate()
                .filter(|(i, _)| mask & (1 << i) != 0)
                .map(|(_, &a)| a)
                .collect();
            let lams: SmallVec<[f64; 8]> = present
                .iter()
                .map(|&a| params[(s.members[a] & !TRACKED_BIT) as usize].lambda)
                .collect();
            let tot: f64 = lams.iter().sum();
            let mut shares: SmallVec<[u32; 8]> = SmallVec::from_elem(0, present.len());
            for _ in 0..s.count {
                let mut x = rng.f64() * tot;
                let mut pick = present.len() - 1;
                for (j, l) in lams.iter().enumerate() {
                    if x < *l {
                        pick = j;
                        break;
                    }
                    x -= l;
                }
                shares[pick] += 1;
            }
            present
                .iter()
                .zip(shares.iter())
                .filter(|(&a, &sh)| sh > 0 && s.members[a] & TRACKED_BIT == 0)
                .map(|(&a, &sh)| (s.members[a], sh))
                .collect()
        })
        .collect();

    let mut covs = uniq_covs;
    for v in draws.into_iter() {
        for (m, sh) in v {
            covs[m as usize].push(sh);
        }
    }
    covs.into_iter()
        .zip(params)
        .map(|(covs, p)| EmAssignment {
            covs,
            lambda: p.lambda,
            pi: p.pi,
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn genome(kmers: Vec<Kmer>, name: &str) -> GenomeSketch {
        GenomeSketch {
            genome_kmers: kmers,
            file_name: name.to_string(),
            k: 31,
            c: 1,
            ..Default::default()
        }
    }

    /// Deterministic Poisson draw (Knuth), for synthetic counts.
    fn pois(rng: &mut fastrand::Rng, lam: f64) -> u32 {
        let l = (-lam).exp();
        let (mut k, mut p) = (0, 1.);
        loop {
            p *= rng.f64();
            if p <= l {
                return k;
            }
            k += 1;
        }
    }

    #[test]
    fn test_ztp_lambda() {
        for lam in [0.3, 1., 2.5, 7.] {
            let m = lam / (1. - f64::exp(-lam));
            assert!((ztp_lambda(m) - lam).abs() < 1e-6);
        }
    }

    /// A single genome at high coverage: λ comes from sylph's coverage estimate --
    /// below a median of 15 (or under `--mean-coverage`) the rate whose
    /// zero-truncated mean is the mean of the contained counts, from there the
    /// median -- and π is the contained fraction over `1 - e^-λ`.
    #[test]
    fn test_em_high_coverage_mirrors_sylph() {
        for (lam, seed) in [(6., 5), (25., 6)] {
            let mut rng = fastrand::Rng::with_seed(seed);
            let kmers: Vec<Kmer> = (0..20_000).collect();
            let mut counts = FxHashMap::default();
            for k in kmers.iter() {
                counts.insert(*k, pois(&mut rng, lam));
            }
            let mut nonzero: Vec<u32> = counts.values().copied().filter(|&c| c > 0).collect();
            nonzero.sort();
            let mean = nonzero.iter().sum::<u32>() as f64 / nonzero.len() as f64;
            let median = nonzero[nonzero.len() / 2] as f64;
            assert!(median > MEDIAN_ANI_THRESHOLD);
            let g = genome(kmers, "x");
            let cands = vec![EmCandidate {
                genome: &g,
                ani: 1.,
                cov: lam,
                median_cov: median,
            }];
            let res = em_reassign(&cands, &counts, false);
            let want = if median < MAX_MEDIAN_FOR_MEAN_FINAL_EST {
                ztp_lambda(mean)
            } else {
                median
            };
            assert!(
                (res[0].lambda - want).abs() < 1e-9,
                "{} vs {}",
                res[0].lambda,
                want
            );
            let frac = nonzero.len() as f64 / 20_000. / (1. - (-want).exp());
            assert!(
                (res[0].pi - frac.min(EM_MAX_PI)).abs() < 1e-9,
                "{}",
                res[0].pi
            );
            let res = em_reassign(&cands, &counts, true);
            assert!(
                (res[0].lambda - ztp_lambda(mean)).abs() < 1e-9,
                "{}",
                res[0].lambda
            );
        }
    }

    /// Two co-present high-coverage genomes sharing a third of their k-mers: both λs
    /// are close to the truth (the modal ratio was biased low), and agree with what
    /// sylph estimates from the counts each genome was given.
    #[test]
    fn test_em_high_coverage_split() {
        let mut rng = fastrand::Rng::with_seed(7);
        let (lx, ly) = (6.0, 9.0);
        let shared: Vec<Kmer> = (0..3000).collect();
        let ux: Vec<Kmer> = (10_000..16_000).collect();
        let uy: Vec<Kmer> = (20_000..26_000).collect();
        let mut counts = FxHashMap::default();
        for k in shared.iter() {
            counts.insert(*k, pois(&mut rng, lx) + pois(&mut rng, ly));
        }
        for k in ux.iter() {
            counts.insert(*k, pois(&mut rng, lx));
        }
        for k in uy.iter() {
            counts.insert(*k, pois(&mut rng, ly));
        }
        let gx = genome([shared.clone(), ux].concat(), "x");
        let gy = genome([shared, uy].concat(), "y");
        let cands = vec![
            EmCandidate {
                genome: &gx,
                ani: 0.999,
                cov: 7.,
                median_cov: 7.,
            },
            EmCandidate {
                genome: &gy,
                ani: 0.99,
                cov: 8.,
                median_cov: 8.,
            },
        ];
        let res = em_reassign(&cands, &counts, false);
        for (r, lam) in res.iter().zip([lx, ly]) {
            assert!(
                (r.lambda - lam).abs() < 0.03 * lam,
                "{} vs {}",
                r.lambda,
                lam
            );
            // sylph reports the non-zero mean of the counts each genome was given.
            let got = r.covs.iter().sum::<u32>() as f64 / r.covs.len() as f64;
            let ztp_mean = r.lambda / (1. - f64::exp(-r.lambda));
            assert!(
                (ztp_mean - got).abs() < 0.02 * got,
                "{} vs {}",
                ztp_mean,
                got
            );
        }
    }

    /// A genome whose k-mers are mostly shared with a far more abundant relative
    /// has a pass-1 median of about the sum of both coverages, so outlier caps from
    /// it let through "unique" k-mers that really carry a relative's repeat (one its
    /// sketch dropped as a duplicate). The caps are refitted from the EM's own
    /// estimate, which keeps them out of the low-coverage genome's λ.
    #[test]
    fn test_em_outlier_caps_refit() {
        let mut rng = fastrand::Rng::with_seed(8);
        let (lx, ly) = (15.0, 3.8);
        let shared: Vec<Kmer> = (0..15_000).collect();
        let ux: Vec<Kmer> = (20_000..22_500).collect();
        let uy: Vec<Kmer> = (30_000..32_500).collect();
        let mut counts = FxHashMap::default();
        for k in shared.iter() {
            counts.insert(*k, pois(&mut rng, lx) + pois(&mut rng, ly));
        }
        for k in ux.iter() {
            counts.insert(*k, pois(&mut rng, lx));
        }
        for (i, k) in uy.iter().enumerate() {
            // 40 of y's unique k-mers carry a 2-copy repeat of x.
            let extra = if i < 40 { pois(&mut rng, 2. * lx) } else { 0 };
            counts.insert(*k, pois(&mut rng, ly) + extra);
        }
        let gx = genome([shared.clone(), ux].concat(), "x");
        let gy = genome([shared, uy].concat(), "y");
        let cands = vec![
            EmCandidate {
                genome: &gx,
                ani: 1.,
                cov: 15.,
                median_cov: 17.,
            },
            EmCandidate {
                genome: &gy,
                ani: 0.999,
                cov: 18.,
                median_cov: 18.,
            },
        ];
        let res = em_reassign(&cands, &counts, false);
        assert!(
            (res[1].lambda - ly).abs() < 0.08,
            "lambda_y = {}",
            res[1].lambda
        );
        assert!(
            (res[0].lambda - lx).abs() < 0.2,
            "lambda_x = {}",
            res[0].lambda
        );
    }

    /// Two co-present genomes sharing a third of their k-mers: EM recovers both
    /// coverages and hands each genome back (nearly) all of its k-mers.
    #[test]
    fn test_em_splits_shared_kmers() {
        let mut rng = fastrand::Rng::with_seed(1);
        let (lx, ly) = (2.0, 3.0);
        let shared: Vec<Kmer> = (0..3000).collect();
        let ux: Vec<Kmer> = (10_000..16_000).collect();
        let uy: Vec<Kmer> = (20_000..26_000).collect();
        let mut counts = FxHashMap::default();
        for k in shared.iter() {
            counts.insert(*k, pois(&mut rng, lx) + pois(&mut rng, ly));
        }
        for k in ux.iter() {
            counts.insert(*k, pois(&mut rng, lx));
        }
        for k in uy.iter() {
            counts.insert(*k, pois(&mut rng, ly));
        }
        let gx = genome([shared.clone(), ux].concat(), "x");
        let gy = genome([shared, uy].concat(), "y");
        let cands = vec![
            EmCandidate {
                genome: &gx,
                ani: 0.999,
                cov: 2.5,
                median_cov: 3.,
            },
            EmCandidate {
                genome: &gy,
                ani: 0.99,
                cov: 3.2,
                median_cov: 3.,
            },
        ];
        let res = em_reassign(&cands, &counts, false);
        assert!((res[0].lambda - lx).abs() < 0.1, "{}", res[0].lambda);
        assert!((res[1].lambda - ly).abs() < 0.1, "{}", res[1].lambda);
        assert!(res[0].pi > 0.97, "{}", res[0].pi);
        assert!(res[1].pi > 0.97, "{}", res[1].pi);
        // Each genome contains ~ (1 - e^-λ) of its 9000 k-mers.
        let ex = 9000. * (1. - f64::exp(-lx));
        let ey = 9000. * (1. - f64::exp(-ly));
        assert!(
            (res[0].covs.len() as f64 - ex).abs() < 0.03 * ex,
            "{}",
            res[0].covs.len()
        );
        assert!(
            (res[1].covs.len() as f64 - ey).abs() < 0.03 * ey,
            "{}",
            res[1].covs.len()
        );
        // Every count is conserved.
        let total: u32 = counts.values().sum();
        let assigned: u32 = res.iter().flat_map(|r| r.covs.iter()).sum();
        assert_eq!(total, assigned);
    }

    /// A near-copy of a present genome that is itself absent (its own k-mers have
    /// zero count) must not take a share of the shared k-mers.
    #[test]
    fn test_em_absent_near_copy_gets_nothing() {
        let mut rng = fastrand::Rng::with_seed(2);
        let shared: Vec<Kmer> = (0..9000).collect();
        let ux: Vec<Kmer> = (10_000..10_500).collect();
        let uy: Vec<Kmer> = (20_000..20_500).collect();
        let mut counts = FxHashMap::default();
        for k in shared.iter().chain(ux.iter()) {
            counts.insert(*k, pois(&mut rng, 3.));
        }
        let gx = genome([shared.clone(), ux].concat(), "x");
        let gy = genome([shared, uy].concat(), "y");
        let cands = vec![
            EmCandidate {
                genome: &gx,
                ani: 0.999,
                cov: 3.,
                median_cov: 3.,
            },
            EmCandidate {
                genome: &gy,
                ani: 0.997,
                cov: 3.,
                median_cov: 3.,
            },
        ];
        let res = em_reassign(&cands, &counts, false);
        assert!(res[1].pi < 0.01, "pi_y = {}", res[1].pi);
        assert!(res[1].covs.len() < 50, "{}", res[1].covs.len());
        assert!((res[0].lambda - 3.).abs() < 0.1);
    }

    /// A low-coverage genome next to an abundant relative: a few of its "unique"
    /// k-mers carry the relative's counts (as happens when the relative's sketch
    /// dropped them as duplicates). The modal-ratio M-step must not let that tail
    /// inflate the rare genome's λ.
    #[test]
    fn test_em_rare_genome_robust_to_contaminated_tail() {
        let mut rng = fastrand::Rng::with_seed(4);
        let (lx, ly) = (2.2, 0.035);
        let shared: Vec<Kmer> = (0..3000).collect();
        let ux: Vec<Kmer> = (10_000..27_000).collect();
        let uy: Vec<Kmer> = (30_000..46_000).collect();
        let mut counts = FxHashMap::default();
        for k in shared.iter() {
            counts.insert(*k, pois(&mut rng, lx) + pois(&mut rng, ly));
        }
        for k in ux.iter() {
            counts.insert(*k, pois(&mut rng, lx));
        }
        for (i, k) in uy.iter().enumerate() {
            // 20 contaminated k-mers carry x's coverage of a 4-copy repeat family
            // (x's sketch drops within-genome duplicates, so to the database they
            // look unique to y).
            let lam = if i < 20 { 4. * lx } else { ly };
            counts.insert(*k, pois(&mut rng, lam));
        }
        let gx = genome([shared.clone(), ux].concat(), "x");
        let gy = genome([shared, uy].concat(), "y");
        let cands = vec![
            EmCandidate {
                genome: &gx,
                ani: 1.,
                cov: 2.2,
                median_cov: 2.,
            },
            EmCandidate {
                genome: &gy,
                ani: 0.997,
                cov: 0.035,
                median_cov: 1.,
            },
        ];
        let res = em_reassign(&cands, &counts, false);
        assert!(
            (res[1].lambda - ly).abs() < 0.012,
            "lambda_y = {}",
            res[1].lambda
        );
        assert!(res[1].pi > 0.8, "pi_y = {}", res[1].pi);
        assert!((res[0].lambda - lx).abs() < 0.1);
    }

    /// The result does not depend on rayon's thread count.
    #[test]
    fn test_em_deterministic() {
        let mut rng = fastrand::Rng::with_seed(3);
        let shared: Vec<Kmer> = (0..2000).collect();
        let ux: Vec<Kmer> = (10_000..13_000).collect();
        let uy: Vec<Kmer> = (20_000..23_000).collect();
        let mut counts = FxHashMap::default();
        for k in shared.iter() {
            counts.insert(*k, pois(&mut rng, 4.));
        }
        for k in ux.iter().chain(uy.iter()) {
            counts.insert(*k, pois(&mut rng, 2.));
        }
        let gx = genome([shared.clone(), ux].concat(), "x");
        let gy = genome([shared, uy].concat(), "y");
        let cands = vec![
            EmCandidate {
                genome: &gx,
                ani: 0.99,
                cov: 2.,
                median_cov: 2.,
            },
            EmCandidate {
                genome: &gy,
                ani: 0.98,
                cov: 2.,
                median_cov: 2.,
            },
        ];
        let run = |t| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(t)
                .build()
                .unwrap()
                .install(|| em_reassign(&cands, &counts, false))
        };
        let (a, b) = (run(1), run(4));
        for (x, y) in a.iter().zip(b.iter()) {
            let (mut cx, mut cy) = (x.covs.clone(), y.covs.clone());
            cx.sort();
            cy.sort();
            assert_eq!(cx, cy);
            assert_eq!(x.lambda.to_bits(), y.lambda.to_bits());
        }
    }
}
