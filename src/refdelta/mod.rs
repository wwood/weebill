//! Reference-based delta compression of read (sample) sketches.
//!
//! A read sketch from an organism present in a reference database shares most of
//! its FracMinHash hashes with that organism's genome sketch. Instead of storing
//! those hashes as ~56-bit values, we store *which* reference hashes are present.
//! Only hashes explained by no reference (sequencing errors, strain divergence,
//! unknown organisms) are stored explicitly.
//!
//! ## K-mer dereplicated reference database
//!
//! To make "which reference" unambiguous, every k-mer in the database is
//! assigned to at most one owner, using a two-level (species rep / strain)
//! scheme that is robust to strains contaminated with sequence from other
//! species:
//!   * Build order places strains of a species *contiguously*, species
//!     representatives first.
//!   * A k-mer present in exactly one species **representative** is owned by that
//!     rep (representatives are higher quality, so they win over strains — a
//!     contaminant k-mer that really belongs to another species' rep is
//!     attributed there, not to the strain carrying the contamination).
//!   * By default, a k-mer in three or more reps is ambiguous → the shared
//!     **pool**. `ref-build --pool-min-genomes` can change that threshold; with
//!     3, exactly two same-tier genomes keep the k-mer assigned to the first
//!     owner.
//!   * A k-mer in no rep but exactly one strain is owned by that strain.
//!   * A k-mer in no rep and enough strains to meet the threshold → the shared
//!     pool.
//!
//! ## Two-stage reference database layout (`.sylref`)
//!
//! Loading a whole reference into a single hash→location table is expensive for
//! large databases. The `.sylref` file is therefore seekable and read in two
//! stages:
//!   * **Stage 1 (sparse, loaded fully):** a 1/N subsample of every genome's
//!     *distinctive* k-mers, stored uncompressed, mapping each sparse hash to its
//!     owning genome. Cheap to load; querying a sample against it yields the small
//!     set of genomes the sample actually contains ("hit genomes").
//!   * **Stage 2 (dense, loaded on demand):** each genome's distinctive set,
//!     minus the hashes already held in stage 1, is an independently Golomb–Rice-
//!     coded block at a known offset (so each k-mer is stored only once). Only the
//!     hit genomes' blocks are decoded and merged back with their stage-1 hashes
//!     to reconstruct the full distinctive array (cached across samples). The
//!     shared **pool** is the one exception — it is conserved across samples and
//!     is loaded once when the index is opened.
//!
//! Missing a genome in stage 1 only costs compression ratio (its hashes fall back
//! to "novel" full-price coding), never correctness: decompression reads the exact
//! genome ids recorded in the sample and loads precisely those dense blocks.
//!
//! ## Encoding a read sketch
//!
//! Read hashes are partitioned into (a) per-genome distinctive hits, (b) pool
//! hits, (c) novel hashes. Each matched set is encoded as positions into the
//! relevant sorted reference array using whichever of {bitmask, present-Rice,
//! absent-Rice} is smallest (a well-covered genome costs ~0 bytes via an empty
//! absent set). The hit genome ids are themselves delta-coded: because strains of
//! a species are contiguous, a sample's hits cluster into small id ranges.
//! Novel hashes use Golomb–Rice (as for normal sketches); counts are varints.
//! The whole payload is zstd-framed.

mod ref_build;
mod sketch_compress;
mod sketch_decompress;

// Shared constants used by both compress and decompress.
pub(crate) const SKETCH_MAGIC: &[u8; 4] = b"SYLD";
pub(crate) const SKETCH_VERSION: u8 = 4;
pub(crate) const SCHEME_BITMASK: u8 = 0;
pub(crate) const SCHEME_PRESENT_RICE: u8 = 1;
pub(crate) const SCHEME_ABSENT_RICE: u8 = 2;
pub(crate) const ZSTD_LEVEL: i32 = 3;

// Utility functions shared by sketch_compress and sketch_decompress.
// They are pure mathematical operations on 2-bit-packed k-mer data.

/// Forward and reverse-complement 2-bit packings of the k-mer starting at base
/// `start` in a 2-bit-packed contig (`data[i/4] >> (2*(i%4)) & 3`).
/// Matches `seeding::fmh_seeds`: forward k-mer has base 0 in the high bits,
/// reverse complement is built from per-base complements.
#[inline]
pub(crate) fn window_fr(data: &[u8], start: usize, k: usize) -> (u64, u64) {
    let mut f = 0u64;
    let mut r = 0u64;
    for j in 0..k {
        let pos = start + j;
        let nuc_f = ((data[pos / 4] >> (2 * (pos % 4))) & 3) as u64;
        let nuc_r = 3 - nuc_f;
        f = (f << 2) | nuc_f;
        r |= nuc_r << (2 * j);
    }
    (f, r)
}

/// The canonical (min of forward / reverse-complement) FracMinHash hash of the
/// k-mer at `start` with the base at offset `off` replaced by `base`. Mirrors the
/// read sketcher so a reconstructed error k-mer hashes identically.
#[inline]
pub(crate) fn substituted_hash(f: u64, r: u64, k: usize, off: usize, base: u8) -> u64 {
    use crate::seeding::mm_hash64;
    let shift_f = 2 * (k - 1 - off);
    let shift_r = 2 * off;
    let f2 = (f & !(3u64 << shift_f)) | ((base as u64) << shift_f);
    let r2 = (r & !(3u64 << shift_r)) | (((3 - base) as u64) << shift_r);
    mm_hash64(f2.min(r2))
}

/// Initialise the global logger (no-op if already initialised by the test
/// harness or a previous call).
fn init_logger(trace: bool) {
    let level = if trace {
        log::LevelFilter::Trace
    } else {
        log::LevelFilter::Info
    };
    let _ = simple_logger::SimpleLogger::new().with_level(level).init();
}

// Public re-exports — everything callers outside `refdelta` need.

pub use ref_build::{
    build_refdb, build_refdb_with_pool_min_genomes, open_ref_index, open_ref_index_for_compress,
    run_ref_build, write_refdb, GenomeSeq, GenomeSource, ReadSeek, RefDb, RefGenome, RefGenomeMeta,
    RefIndex,
};
pub use sketch_compress::{
    compress_seq, compress_seq_with_meta, compress_seq_with_screen_ani,
    compress_seq_with_screen_ani_and_error_kmers, compress_seq_with_screen_ani_and_telemetry,
    run_ref_compress, RefCompressTelemetry,
};
pub use sketch_decompress::{decompress_seq, decompress_seq_with_meta};
