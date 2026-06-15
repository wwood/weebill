# Two-stage profiling (`.syl2db`)

This document describes the two-stage profiling feature: the new algorithm, the
architectural changes in the codebase, and the new/changed CLI surface.

## Motivation

Single-stage `sylph profile` loads the **entire** dense genome database
(`.syldb`) into memory and scores the sample against every genome. For a
comprehensive reference (e.g. a full GTDB release, hundreds of thousands of
genomes) that database is tens of gigabytes, and both the peak RAM and a large
share of the wall-clock are spent loading/holding genomes that do not occur in
the sample at all. A typical shotgun sample only actually contains a small
fraction of the reference.

Two-stage profiling keeps the same dense scoring math (so results match) but
only pays the dense cost for genomes that could plausibly be present.

## The idea

1. **Stage 1 — screen (cheap, in memory).** Score the sample against a small,
   heavily sub-sampled *sparse* index of every genome and keep the ones whose
   coverage-adjusted ANI clears a permissive threshold (`--screen-ani`,
   default 85). This is essentially a `query` at a coarse `-c`.
2. **Stage 2 — dense (only survivors, from disk).** For each survivor, read and
   decode just that genome's dense block from a seekable on-disk file, run the
   full profiling math on it, and keep it only if it passes the profiling ANI
   threshold. The genomes that fail are dropped immediately and never cached.
   The usual k-mer reassignment and abundance estimation then run on the kept
   genomes.

The full dense database is **never** loaded into memory. Peak RAM is the small
sparse index plus the dense blocks of the (few) genomes that survive — not the
whole reference.

## The `.syl2db` file format

A `.syl2db` is produced from a standard `.syldb` by `sylph db-convert`. Layout:

```
[ header ]  magic "SY2D" (4 B) | version (1 B) | footer_offset (u64 LE)
[ body   ]  per genome, in database order:
              dense genome_kmers block        (Golomb-Rice)
              has-pseudotax flag (1 B)
              dense pseudotax block (optional) (Golomb-Rice)
[ footer ]  bincoded stage-1 index:
              c, k, screen_c, and per genome:
                file_name, first_contig_name, gn_size, min_spacing,
                dense_offset (into body), and the sparse k-mer subsample
```

- **Dense blocks** keep *every* k-mer at the source database's `-c` (so the
  dense stage is exactly as accurate as single-stage), delta + Golomb-Rice
  coded — near the information-theoretic size for random FracMinHash hashes.
- **The footer** is the stage-1 screen index: each genome's k-mers sub-sampled
  to `screen_c` (a strict subset of the dense k-mers; FracMinHash lets you make
  a sketch sparser for free) plus the `dense_offset` needed to find its block.
  Only the footer is read into memory when the database is opened.

## Profiling algorithm (`profile --two-stage`)

```
open(.syl2db):
    read header + footer  ->  in-memory sparse screen index (dense body stays on disk)

for each sample:
    sketch / load the sample at -c (must be <= dense -c)

    # stage 1: screen (parallel over all genomes, in memory)
    survivors = { g : adjusted_ANI(sample, sparse_index[g]) >= screen_ani }

    # stage 2: dense (parallel over survivors only)
    kept = []
    for g in survivors:
        block = pread(dense block of g)          # positional, no shared cursor
        sketch = decode(block)                    # short-lived
        if profiling_passes(sample, sketch):      # full dense get_stats
            kept.push(sketch)                     # keep only winners; drop the rest

    # unchanged downstream
    reassign k-mers (pseudotax) over kept
    estimate coverage / relative abundance
    write profile rows
```

Key properties:

- **Only survivors are decoded.** The dense cost scales with the number of
  genomes that pass the screen, not with the database size.
- **Decode-and-drop.** Each survivor is decoded into a short-lived sketch and
  freed unless it passes profiling, so peak RAM is proportional to genomes that
  actually survive — not to the (much larger) screen-survivor set. (Decoding to
  a contiguous `Vec` and then probing the sample is deliberate: it lets the CPU
  prefetch and hide the sample-hashmap latency; a no-materialization streaming
  variant was measured ~2× slower.)
- **Lock-free parallelism.** Dense blocks are read with positional reads
  (`pread`), so any number of threads decode concurrently with no shared file
  cursor and no lock; only the touched blocks (plus reclaimable OS page cache)
  cost memory, keeping RSS low at every thread count.

## Accuracy / tuning

The dense stage uses the same scoring as single-stage, so a genome that
survives the screen gets an identical profile. The only way two-stage can differ
from single-stage is if the screen drops a genuinely present genome. The screen
is therefore deliberately permissive:

- `--screen-ani` (default **85**): lower is more sensitive (more survivors →
  more dense cost), higher is cheaper but can drop low-abundance genomes. 85 was
  chosen empirically to recover essentially the full single-stage profile.
- `db-convert --screen-c` (default **3000**): the screen index sub-sampling
  rate. Sparser (larger) → smaller/faster screen and smaller `.syl2db`, but a
  noisier screen, which should be paired with a permissive `--screen-ani`.

## Architectural changes

- **`src/twostage_db.rs` (new).** The `.syl2db` format end to end:
  - `BitWriter`/`BitReader` and `write_hashes`/`read_hashes` — delta +
    Golomb-Rice codec for sorted hash sets.
  - `write_two_stage_db` — re-pack a set of `GenomeSketch`es into the seekable
    layout.
  - `TwoStageDb` — an opened database. Holds the in-memory sparse `screen_sketches`
    and reads dense blocks on demand via `decode_dense` (uncached) / `load_dense`
    (cached), using positional reads bounded by adjacent `dense_offset`s.
  - `run_db_convert` — the `db-convert` subcommand handler.
- **`src/contain.rs`.**
  - `compute_dense_survivors` — runs stage 1 (screen) and stage 2 (decode-drop)
    and returns the dense sketches of genomes that pass pass-1 profiling.
  - `get_stats` refactored to share its coverage-correction/ANI math with the
    two-stage path via `finalize_stats` (output unchanged).
  - `contain` opens a `.syl2db` (when the input is one) and switches on the
    two-stage path.
- **`src/cmdline.rs`.** New `DbConvert` subcommand + `DbConvertArgs`; new
  `profile` flags (`--two-stage`, `--dense-c`, `--screen-c`, `--screen-ani`,
  `--dense-cache`).
- **`src/constants.rs`.** `SCREEN_C_DEFAULT = 3000`, `SCREEN_MIN_ANI_DEFAULT = 85`,
  `DENSE_C_DEFAULT = 50`.

No dependencies are added (the format uses a hand-rolled bit codec; dense blocks
are read with `std::os::unix::fs::FileExt::read_at`).

## CLI changes

### New subcommand: `sylph db-convert`

Convert a standard `.syldb` into a seekable two-stage `.syl2db`:

```
sylph db-convert gtdb.syldb -o gtdb [--screen-c 3000] [-t THREADS]
# -> gtdb.syl2db
```

| flag | default | meaning |
|---|---|---|
| `-o, --output` | — | output name (`.syl2db` appended) |
| `--screen-c` | 3000 | sub-sampling rate of the stage-1 screen index (≥ db `-c`) |
| `-t` | 3 | threads |

### `sylph profile` additions

```
# single-stage (unchanged): loads the whole dense database
sylph profile gtdb.syldb sample.sylsp

# two-stage: screen against the sparse index, dense-profile only survivors
sylph profile --two-stage gtdb.syl2db sample.sylsp
```

| flag | default | meaning |
|---|---|---|
| `--two-stage` | off | enable two-stage profiling |
| `--screen-ani` | 85 | min adjusted ANI (0–100) to pass the stage-1 screen |
| `--screen-c` | db's `-c` | stage-1 screen sub-sampling rate (≥ db `-c`) |
| `--dense-c` | from `.syl2db` | dense second-stage `-c` (taken from the `.syl2db`) |
| `--dense-cache` | — | dir of cached per-genome dense sketches, reused across runs |

Output format, columns, and all other profiling options are unchanged.
