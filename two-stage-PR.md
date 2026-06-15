# Add two-stage profiling (`db-convert` + `profile --two-stage`)

## Summary

This PR adds an optional **two-stage profiling** mode that profiles a sample at
full dense resolution **without ever loading the whole dense database into
memory**. Today, `sylph profile` loads the entire `.syldb` (tens of GB for a
full GTDB) and scores every genome, even though a sample contains only a small
fraction of the reference. Two-stage instead:

1. **screens** the sample against a small in-memory *sparse* index of every
   genome (a coarse `query`), then
2. **densely profiles only the genomes that pass the screen**, decoding their
   per-genome blocks on demand from a seekable on-disk database.

Results match single-stage profiling, while profiling is **3–12× faster** and
uses **7–11× less RAM**, and the on-disk database is **~28% smaller**.

New surface (single-stage profiling is unchanged):

- `sylph db-convert in.syldb -o out` → builds a seekable `out.syl2db`.
- `sylph profile --two-stage out.syl2db sample.sylsp` → screen-then-dense.

See [`docs/two-stage-profiling.md`](docs/two-stage-profiling.md) for the
algorithm, file format, and architecture.

## Benchmarks

GTDB r220-scale reference: **199,923 genomes**, dense `-c 100`. Two samples of
differing complexity, profiled at 1 and 8 threads. Dense `.syl2db` built with
defaults (`--screen-c 3000`); profiled with default `--screen-ani 85`. Reads
were pre-sketched once and reused, so the numbers isolate database loading +
profiling. Peak RAM is `maxRSS` from `/usr/bin/time -v`.

### Wall-clock and peak RAM

| sample | threads | single-stage | two-stage | speedup | single RAM | two-stage RAM | RAM saving |
|---|---|---|---|---|---|---|---|
| marine (diverse, deep-ocean) | 1 | 192 s | 151 s | 1.3× | 53.9 GB | **7.2 GB** | 7.5× |
| marine | 8 | 89 s | **24 s** | **3.7×** | 53.9 GB | **7.2 GB** | 7.5× |
| human gut | 1 | 125 s | 23 s | 5.4× | 51.3 GB | **4.8 GB** | 10.7× |
| human gut | 8 | 75 s | **6 s** | **12.5×** | 51.3 GB | **4.8 GB** | 10.7× |

- The **RAM win is large and constant** (~5–7 GB vs ~52 GB): two-stage holds only
  the small sparse screen index plus the dense blocks of surviving genomes, and
  positional reads keep RSS flat across thread counts.
- The **speedup scales with how few genomes survive the screen.** The human gut
  sample matches a small slice of GTDB (few survivors) → up to 12.5×. The diverse
  marine sample has many survivors, so the win is smaller (3.7× at 8 threads) but
  still substantial, and its RAM win is just as large.

### Database size (storage)

| | size |
|---|---|
| single-stage `.syldb` | 51.7 GB |
| two-stage `.syl2db` | **37.2 GB** |

The `.syl2db` is ~28% smaller because dense blocks are delta + Golomb-Rice coded
(near the entropy floor for FracMinHash hashes) whereas the `.syldb` stores raw
sketches. `db-convert` is a one-time step run once per database from an existing
`.syldb`.

### Accuracy (two-stage vs single-stage, 8 threads)

| sample | single species | two-stage species | shared | only single | only two | max \|Δabund\| | sum \|Δabund\| |
|---|---|---|---|---|---|---|---|
| marine | 710 | 711 | 710 | 0 | 1 | 0.0040 % | 0.0081 % |
| human gut | 169 | 169 | 169 | 0 | 0 | 0.0000 % | 0.0000 % |

Two-stage recovers **every** species single-stage reports (no losses), with
abundance differences at the 4th decimal. The single "only two" genome on the
marine sample is a trace-level genome at the detection boundary.

## How it works (brief)

`db-convert` re-packs a `.syldb` into a `.syl2db`: a small bincoded **stage-1
sparse index** (every genome sub-sampled to `--screen-c`) plus per-genome
**Golomb-Rice-compressed dense blocks** at the database's own `-c`, addressed by
offset. `profile --two-stage` loads only the sparse index, screens the sample to
get candidate genomes (adjusted ANI ≥ `--screen-ani`, default 85), then for each
candidate positionally reads and decodes its dense block, runs the normal dense
scoring, and keeps it only if it passes — dropping the rest immediately. Dense
reads use `pread`, so decoding is lock-free and parallel with low RSS at any
thread count.

Design choices validated by benchmarking (details in the doc):

- **Decode-and-drop** (decode each survivor to a short-lived `Vec`, probe, free
  unless it passes) rather than retaining all survivors — same speed as
  retaining, but peak RAM scales with genomes that survive, not the screen set.
- **`pread`, not `mmap`** — same parallel decode, but no whole-file mapping, so
  RSS stays low at every thread count.
- The defaults `--screen-ani 85` / `--screen-c 3000` were chosen from a sweep as
  the point that keeps full single-stage recovery at minimum cost.

## Implementation notes

- New module `src/twostage_db.rs`: the `.syl2db` format (writer, seekable
  reader, Golomb-Rice codec, `db-convert` handler).
- `src/contain.rs`: stage-1 screen + stage-2 decode-drop (`compute_dense_survivors`),
  with the coverage/ANI math shared with single-stage via `finalize_stats`
  (output validated byte-identical to the previous `get_stats`).
- `src/cmdline.rs`: `db-convert` subcommand and the `profile` two-stage flags.
- **No new dependencies** (hand-rolled bit codec; `std` positional reads).
- Single-stage profiling, the `.syldb`/`.sylsp` formats, and the profile output
  are all unchanged; two-stage is opt-in.

## Usage

```bash
# one-time: convert an existing dense database
sylph db-convert gtdb_r220.syldb -o gtdb_r220 -t 50      # -> gtdb_r220.syl2db

# profile (screen at ANI 85, dense-profile only survivors)
sylph profile --two-stage gtdb_r220.syl2db sample.sylsp -t 8 -o profile.tsv
```
