# Changelog

Weebill is a fork of [sylph](https://github.com/bluenote-1577/sylph); this file records
changes made in weebill.

## Unreleased

### Added

- `profile --apply-unknown <profile.tsv>` converts an existing profile made *without* `-u`/`--estimate-unknown` into the profile `-u` would have produced, without re-profiling. Everything `-u` changes is a per-sample scalar applied to two printed columns (`Eff_cov` is scaled by `read_length/(read_length-k+1) / identity` and relabelled `True_cov`; `Sequence_abundance` is scaled by the estimated covered-bases fraction), while `Taxonomic_abundance` and the ANIs are unchanged. The same pre-sketched sample(s) and the original database are both required on the command line — genome sizes come from the database and the k-mer-count distribution/read length from the sample sketch — so the result matches a real `-u` run to the precision of the input TSV's columns — but is not guaranteed to be bit-for-bit identical, since only the rounded printed columns are available (not the full-precision internal values), so individual `True_cov`/`Sequence_abundance` cells may differ from a real `-u` run by a unit in the last printed place. Works with plain `.syldb`, two-stage `.syl2db` (no dense blocks are decoded) and reference-delta `.sylspr` samples (via `--reference`); the TSV may be gzip-compressed. `--individual-records` databases are rejected, because a profile TSV identifies each genome only by its `Genome_file` and several records then share one, so their sizes cannot be told apart.

### Changed

- Read sketching is now multi-threaded *within* a single input, not just across inputs. A dedicated reader thread handles IO/decompression while rayon workers extract k-mers from batched reads in parallel; only the order-dependent dedup fold stays serial and in file order, so sketches remain byte-for-byte identical to the previous single-threaded output and independent of `-t` (verified with `diff -r` of single-end, paired and interleaved sketches and of `profile --merge` output, at 1 and 4 threads). A single large read file now scales across cores (~2× on 4 threads for a single-end input) where it previously did all k-mer work on one core regardless of `--threads`.

## Version 0.2.0

### Added

- `sketch --reference` now takes `--ref-screen-ani`, `--min-dense-kmers-for-error` and `--no-error-kmer`, the reference-delta compression tunables that were previously only settable on `ref-compress`. Sketching straight to `.sylspr` used the built-in defaults with no way to change them.
- The zstd frames in `.sylspc`/`.syldbc` and `.sylspr` now carry a trailing content checksum, so truncated or bit-rotted sketches are detected on read instead of decoding into a silently corrupt sketch. Readers drain the frame to its end, which is what forces zstd to validate the checksum.
- The seekable databases `.sylref` and `.syl2db` now store an XXH64 of the whole file in their header. They are read a block at a time, so nothing validates them end to end in normal use; `weebill inspect` accepts both and verifies the checksum (reporting `checksum: ok`, or failing with a non-zero exit if the file is corrupt).

- The error-k-mer scan is now gated on what it is predicted to be worth, by two new options:
  - `--min-coverage-for-error` (default 0.1x) decides **which genomes** to scan, on their estimated coverage depth. Scanning a genome costs a pass over its sequence per chunk of novel k-mers -- proportional to genome length, flat in depth -- while the k-mers recovered scale with depth x length. The return per unit of work is therefore set by depth alone, and genome size cancels; `--min-dense-kmers-for-error` conflated the two and could not express this.
  - `--min-error-kmer-shrink` (default 0.10) decides **whether to scan at all**, on the fraction by which the scan is predicted to shrink the output. The scan's fixed cost scales with the sample's novel k-mers, so on a diverse, shallowly covered sample it can spend minutes to take 2% off the file, while on a high-coverage one it takes a fifth off in seconds.

  The prediction combines sequencing errors (which scale with a genome's summed sample counts) and strain SNPs against the reference (which scale with the genome k-mers covered *and* the genome's stage-1 ANI divergence -- soil organisms sit much further from their GTDB representative than human gut strains, and yield far more recodable k-mers per k-mer covered). Calibrated across 18 metagenomes spanning human, marine, soil and bioreactor: the defaults skip 95% of the total scan time to give up a third of the total saving, running the scan on every sample it shrinks by a tenth or more (e.g. human gut, 2.1MB -> 1.6MB in 0.9s) and skipping the rest (e.g. a soil sample that spent 222s to take 1.9% off a 162MB file).
- `ref-compress --telemetry` reports two new per-genome columns, `coverage_depth` and `expected_error_kmers`, so the thresholds above can be recalibrated against what the scan actually recovers.

### Fixed

- The error-k-mer scan no longer does a minute or more of work when no genome is eligible to be scanned. It built its per-chunk novel-k-mer index and bloom filter before checking that it had any genome sequence to scan them against, so a sample with tens of millions of novel k-mers paid the full setup cost to find nothing (78s on a soil metagenome with 38M novel k-mers and zero eligible genomes).

### Changed

- **Breaking:** the `.sylspc`/`.syldbc` format is now version 5 and `.sylspr` is version 5; older files are refused rather than read. Re-sketch (or re-run `ref-compress`) to upgrade.
- **Breaking:** `.sylref` is now version 6 and `.syl2db` version 3 (checksum in the header); older files are refused rather than read. Rebuild with `ref-build` / `db-convert` to upgrade.
- A compressed sample sketch that cannot be read now reports the underlying error (e.g. `Restored data doesn't match checksum`, `incomplete frame`) instead of always blaming an incompatible version.

## Version 0.1.0

### Added

- Renamed the binary to `weebill` (library crate remains `sylph`); README documents the fork, installation, and new sub-commands.
- Added `sketch --merge`/`profile --merge` to combine multiple read inputs (single, paired, interleaved) into one sketch, with strict validation: rejects unusable inputs, directory-like output paths, and legacy `.sylsp` merge targets rather than silently dropping data.
- Added `sketch --tolerate-empty-inputs` so a zero-read input stream (e.g. an empty FIFO) sketches as empty instead of failing; `sketch` now errors if *every* read stream turns out empty, to catch broken upstream producers.
- Paired-end sketching now averages read length across both mates (excluding sub-k mates) instead of using one mate; documented the short-insert overlap limitation.
- Added a compressed sketch format (`--compressed-output`/`--compressed-database`, `.sylspc`/`.syldbc`): sorted, delta-encoded, Golomb–Rice-coded hashes wrapped in a zstd frame — ~30%+ smaller databases, ~55% smaller samples. `query`/`profile`/`inspect` transparently read both legacy and compressed formats.
- Added reference-based delta compression (`weebill ref-build` + `weebill ref-compress`): `ref-build` produces a k-mer-dereplicated, two-stage seekable reference (`.sylref`) from a genome database; `ref-compress` encodes samples against it (`.sylsp` → `.sylspr`), storing only novel/error k-mers explicitly. ~60% smaller than `.sylspc`, ~80% smaller than bincode, losslessly. `ref-build` streams and parallelizes with bounded RAM (`--max-ram`) via on-disk hash-partitioned shards. `query`/`profile` can read `.sylspr` directly given `--reference`.
- Added an AVX-512 k-mer extraction path for FracMinHash sketching (8 k-mers/iteration via `vpcompressq`), selected at runtime with AVX2/scalar fallback. Bit-for-bit identical output to AVX2/scalar; ~1.6–2× faster than AVX2 depending on subsampling rate.
- Pinned `rand` to 0.9 and seed the cuckoo filter RNG for reproducible sketching.
- Sketch read passes (paired/interleaved/single) are now drained concurrently.
- Added CI test workflow, `cargo fmt` check, pre-commit hook, and `ARCHITECTURE.md`.
- Moved to `bird_tool_utils` project conventions.

## Previous
Weebill was forked from sylph commit cf6ee06 (~0.9.0).
