pub const EM_ABUND_CUTOFF: f64 = 0.01;
pub const PAIR_REGEX: &str = r"(.+)(_?1|_?2)(\..+)";
pub const CUTOFF_PVALUE: f64 = 0.9999999999;
pub const SAMPLE_SIZE_CUTOFF: usize = 25;
pub const MEDIAN_ANI_THRESHOLD: f64 = 2.;
pub const QUERY_FILE_SUFFIX: &str = ".syldb";
pub const SAMPLE_FILE_SUFFIX: &str = ".sylsp";
pub const QUERY_COMP_FILE_SUFFIX: &str = ".syldbc";
pub const SAMPLE_COMP_FILE_SUFFIX: &str = ".sylspc";
pub const REF_DB_SUFFIX: &str = ".sylref";
pub const REF_SAMPLE_SUFFIX: &str = ".sylspr";
pub const QUERY_FILE_SUFFIX_VALID: [&str; 3] =
    [QUERY_FILE_SUFFIX, ".sylqueries", QUERY_COMP_FILE_SUFFIX];
pub const SAMPLE_FILE_SUFFIX_VALID: [&str; 3] =
    [SAMPLE_FILE_SUFFIX, ".sylsample", SAMPLE_COMP_FILE_SUFFIX];
pub const MIN_ANI_DEF: f64 = 0.9;
pub const MIN_ANI_P_DEF: f64 = 0.95;
pub const MAX_MEDIAN_FOR_MEAN_FINAL_EST: f64 = 15.;
pub const DEREP_PROFILE_ANI: f64 = 0.975;
pub const MAX_DEDUP_COUNT: u32 = 4;
pub const MAX_DEDUP_LEN: usize = 10000000;
pub const DEFAULT_FPR: f64 = 0.0001;
pub const MED_KMER_FOR_ID_EST: f64 = 3.;
pub const DENSE_C_DEFAULT: usize = 50;
pub const SCREEN_C_DEFAULT: usize = 3000;
pub const SCREEN_MIN_ANI_DEFAULT: f64 = 85.;
pub const REF_SPARSE_C_DEFAULT: usize = 3000;
pub const REF_SCREEN_ANI_DEFAULT: f64 = 87.;
pub const MIN_DENSE_KMERS_FOR_ERROR_DEFAULT: usize = 100;
/// Minimum estimated coverage depth (lambda) a hit genome needs before it is scanned for
/// single-substitution error k-mers. Scanning a genome costs one pass over its sequence per
/// chunk of novel k-mers -- proportional to genome length, with no dependence on depth --
/// while the k-mers it yields go as depth x length (see the yield coefficients below). The
/// return per unit of work is therefore a function of depth alone, which is what this
/// thresholds. Shallow genomes cost exactly as much to scan as deep ones and pay back
/// almost nothing.
///
/// The default is deliberately lax: on a human gut sample it halves the genomes scanned for
/// 0.9% more output. Most of the protection against a pathological scan comes from
/// [`MIN_ERROR_SHRINK_DEFAULT`], not from here.
pub const MIN_COVERAGE_FOR_ERROR_DEFAULT: f64 = 0.1;
/// Minimum fraction by which the error-k-mer scan must be predicted to shrink the output
/// before it is run at all. Beyond the per-genome cost gated by the coverage threshold, the
/// scan has a fixed cost that scales with the novel set (indexing every novel k-mer under both
/// of its half-k-mer keys, plus a bloom filter over that index) -- on a soil metagenome with
/// 38M novel k-mers that alone is over a minute. Paying it to shrink the output by 1% is a bad
/// trade; on a human gut sample, where it takes a fifth off, it is an excellent one.
///
/// Across 18 metagenomes (human, marine, soil, bioreactor) this default runs the scan on every
/// sample whose output it shrinks by a tenth or more, and skips the rest -- dropping 95% of the
/// total scan time while keeping a third of the total saving. Lower it to trade time for space.
pub const MIN_ERROR_SHRINK_DEFAULT: f64 = 0.10;
/// How much the output actually shrinks, as a multiple of the share of novel k-mers the scan
/// recodes. A recoded k-mer moves from a Rice-coded novel hash to a (position, base) entry
/// against a genome, and empirically that swap is worth ~2 bytes; the novel section dominates
/// the payload, so the two are proportional. Measured at 0.45 across the same 18 metagenomes
/// (e.g. recoding 48% of the novel k-mers of a human gut sample took 22% off the file).
pub const SHRINK_PER_RECODED_SHARE: f64 = 0.45;
/// Expected recodable k-mers per k-mer *observation* drawn from a genome (i.e. per unit of
/// summed sample count). These are sequencing errors: a fraction 1-(1-e)^k of the reads'
/// k-mers carry a substitution, and each is (almost always) a distinct novel hash one
/// substitution away from the genome. The value implies e ~ 0.003/base.
pub const ERROR_YIELD_PER_OBSERVATION: f64 = 0.09;
/// Overall calibration of the predicted yield against what the scan really recovers, applied
/// to the sum of the sequencing-error and strain-SNP terms.
///
/// The SNP term needs no constant of its own: where the sampled strain differs from the
/// reference at one base, each of the k read k-mers spanning that site sits a single
/// substitution from the reference and so recodes, giving `k * divergence` recodable k-mers
/// per covered genome k-mer -- and the divergence is already measured, per genome, as the
/// stage-1 sparse-screen ANI. That matters because divergence varies hugely by biome: against
/// GTDB, human gut strains sit at ~99% ANI while soil organisms are far further from their
/// nearest representative, and a fixed coefficient fitted on one under-predicts the other by
/// an order of magnitude.
///
/// Calibrated across 18 metagenomes (human, marine, soil, bioreactor) against the k-mers the
/// scan actually recovered; it predicts the right side of the yield gate for 17 of them. It
/// feeds only the scan's cost/benefit gates and its telemetry -- nothing about the encoding
/// depends on it, so a bad estimate can cost time or a missed saving, never correctness.
pub const ERROR_YIELD_CALIBRATION: f64 = 0.6;
pub const GENOME_SKETCH_SUFFIX: &str = ".sylgn";
/// Two-stage seekable database: a small bincoded sparse (screen) index plus
/// per-genome Golomb-Rice compressed dense blocks loaded on demand.
pub const TWO_STAGE_DB_SUFFIX: &str = ".syl2db";
