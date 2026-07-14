# Changelog

Weebill is a fork of [sylph](https://github.com/bluenote-1577/sylph); this file records
changes made in weebill.

## Unreleased

### Added

- `sketch --reference` now takes `--ref-screen-ani`, `--min-dense-kmers-for-error` and `--no-error-kmer`, the reference-delta compression tunables that were previously only settable on `ref-compress`. Sketching straight to `.sylspr` used the built-in defaults with no way to change them.
- The zstd frames in `.sylspc`/`.syldbc` and `.sylspr` now carry a trailing content checksum, so truncated or bit-rotted sketches are detected on read instead of decoding into a silently corrupt sketch. Readers drain the frame to its end, which is what forces zstd to validate the checksum.
- The seekable databases `.sylref` and `.syl2db` now store an XXH64 of the whole file in their header. They are read a block at a time, so nothing validates them end to end in normal use; `weebill inspect` accepts both and verifies the checksum (reporting `checksum: ok`, or failing with a non-zero exit if the file is corrupt).

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
