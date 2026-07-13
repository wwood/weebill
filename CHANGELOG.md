# Changelog

Weebill is a fork of [sylph](https://github.com/bluenote-1577/sylph); this file records
changes made in weebill.

## Unreleased

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
