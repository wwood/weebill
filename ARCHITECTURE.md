# Architecture

## Overview

Sylph (binary: `weebill`) performs ultrafast genome ANI queries and taxonomic
profiling of metagenomic shotgun samples. The core approach is FracMinHash
sketching: k-mers whose hash falls below `u64::MAX / c` form a uniform
subsample of a sequence, enabling cheap containment estimates between reads
and genomes. Coverage is estimated from the k-mer multiplicity distribution
(Poisson / negative-binomial model in `inference`) and used to correct the
naive containment fraction into an accurate ANI.

## Module Map

Start here when reading the code:

| Module | Responsibility | Key items |
|---|---|---|
| `types.rs` | Shared data types and hash utilities | `SequencesSketch`, `GenomeSketch`, `AniResult`, `Kmer = u64`, `mm_hash` |
| `seeding.rs` | Rolling FracMinHash k-mer extraction | `fmh_seeds`, `mm_hash64`, `rev_mm_hash64`, `rev_hash_64` |
| `avx2_seeding.rs` | AVX2 SIMD path for k-mer hashing (x86_64 only) | `mm_hash256`, `extract_markers_avx2` |
| `sketch.rs` | Sketching pipeline for reads and genomes | `sketch`, `write_read_sketch_file`, `sketch_genome` |
| `contain.rs` | Query and profile modes; two-stage workflow | `contain`, `subsample_view`, `densify_genome` |
| `inference.rs` | Coverage estimation from k-mer multiplicity distributions | `mme_lambda`, `mle_zip`, `binary_search_lambda` |
| `compress.rs` | Compressed sketch I/O (`SYLZ` format, Golomb-Rice + zstd) | `peek_is_compressed`, `write_genome_sketch_compressed`, `read_seq_sketch_compressed` |
| `refdelta/mod.rs` | Module root: shared constants, `window_fr`/`substituted_hash` utilities, public re-exports | `SKETCH_MAGIC`, `SKETCH_VERSION`, `SCHEME_*`, `ZSTD_LEVEL` |
| `refdelta/ref_build.rs` | Building and writing `.sylref`; reading via `RefIndex` | `RefDb`, `RefIndex`, `build_refdb`, `write_refdb`, `open_ref_index`, `run_ref_build` |
| `refdelta/sketch_compress.rs` | Compression of sample sketches into `.sylspr` | `RefCompressTelemetry`, `compress_seq`, `encode_subset`, `find_error_kmers`, `run_ref_compress` |
| `refdelta/sketch_decompress.rs` | Decompression of `.sylspr` back to sample sketches | `decompress_seq`, `decompress_seq_with_meta`, `decode_subset` |
| `twostage_db.rs` | Two-stage seekable genome database (`.syl2db`) | `TwoStageDb`, `run_db_convert` |
| `merge.rs` | Merge multiple sample sketches into one | `merge` |
| `inspect.rs` | YAML metadata inspection of sketch files | `inspect` |
| `cmdline.rs` | All CLI argument structs (clap derive) | `SketchArgs`, `ContainArgs`, `RefBuildArgs`, `RefCompressArgs`, `DbConvertArgs` |
| `constants.rs` | File-suffix constants and algorithm defaults | file suffixes (`.syldb`, `.sylsp`, …), `MIN_ANI_DEF`, `DENSE_C_DEFAULT` |
| `main.rs` | Entry point; dispatches to modules by subcommand | — |

`contain.rs` depends on `inference`, `sketch`, `types`, and both database
formats. `refdelta` (the module) and `twostage_db.rs` depend on `compress` and
`seeding` but are otherwise independent of each other. Within `refdelta`,
`ref_build` owns the database types and I/O; `sketch_compress` and
`sketch_decompress` are siblings that access `ref_build` internals via
`pub(crate)` visibility. Everything else is roughly leaf-level.

## Key Types & Data Flow

**`Kmer = u64`** — a single FracMinHash value. A k-mer is included in a sketch
iff `mm_hash64(canonical_kmer) < u64::MAX / c`. Stored directly; the raw
sequence is not kept.

**`SequencesSketch`** — a read/sample sketch: `kmer_counts: FxHashMap<Kmer,
u32>` (hash → read count), plus `c`, `k`, `mean_read_length`, and file
metadata. Serialized to `.sylsp` (bincode) or `.sylspc` / `.sylspr`
(compressed formats).

**`GenomeSketch`** — a genome sketch: `genome_kmers: Vec<Kmer>` (sorted) plus
an optional `pseudotax_tracked_nonused_kmers` set used to reassign shared
k-mers during profiling. Serialized to `.syldb` (bincode) or `.syldbc` (SYLZ
compressed).

**`AniResult<'a>`** — the output of one genome–sample comparison: naive and
adjusted ANI, coverage (`final_est_cov`), confidence intervals, containment
index, and optional relative abundance. Borrows the source `GenomeSketch`.

**`RefIndex`** (in `refdelta/ref_build.rs`) — the seekable query-time form of a
`.sylref` file. Stage-1 sparse MPHF (boomphf) maps k-mer hashes to genome ids,
loaded fully. Dense per-genome Golomb-Rice blocks are loaded on demand by file
seek.

**`TwoStageDb`** (in `twostage_db`) — the seekable query-time form of a
`.syl2db` file. Footer holds bincoded sparse `GenomeSketch` subsets; body
holds Golomb-Rice dense blocks accessed by stored byte offsets.

### Data flow: build time

```
FASTA → fmh_seeds (seeding) → Vec<Kmer> → GenomeSketch → .syldb / .syldbc
FASTQ → cuckoo-filter dedup → fmh_seeds → FxHashMap<Kmer,u32>
                                         → SequencesSketch → .sylsp / .sylspc / .sylspr
```

### Data flow: query / profile

```
.syldb → Vec<GenomeSketch>        ]
.sylsp → SequencesSketch          ] → intersect k-mers → coverage vec
                                       → inference::mme_lambda → λ → adjusted ANI
                                       → AniResult → TSV
```

### Data flow: two-stage profiling

```
.syl2db or .sylref  (stage-1, sparse, always in RAM)
  + SequencesSketch
  → screen: genomes whose sparse k-mers hit the sample
  → for each hit: seek + decode dense block → full GenomeSketch (cached)
  → standard dense profiling pass
```

### Data flow: ref-delta compression

```
ref-build:  .syldb files → k-mer dereplication (distinctive sets + shared pool)
            → MPHF over sparse subset → seekable .sylref

ref-compress: .sylsp + .sylref stage-1 → find hit genomes → load dense blocks
              → partition hashes: (per-genome distinctive | pool | novel)
              → delta-encode as positions into sorted reference arrays
              → zstd-frame → .sylspr
```

## Design Decisions & Invariants

**FracMinHash (threshold) not MinHash (top-k).** Because `hash < u64::MAX / c`
defines membership, making a sketch coarser is a strict subset operation:
`subsample_view` (in `contain.rs`) drops hashes above the new threshold without
re-reading any sequence. The reverse is impossible—densifying requires
re-sketching the source FASTA (`densify_genome`).

**Hashes, not k-mers, are stored.** Containment only needs membership testing;
the actual sequence is irrelevant. The hash value is the only thing written to
disk and compared at query time.

**`mm_hash64` has a known first-step bug** (uses `+` where the minimap2
algorithm intends `-(key)-1`). This is preserved intentionally for backwards
sketch compatibility; the comment says "TODO: fix after release". `rev_mm_hash64`
exactly inverts *this* buggy function; `rev_hash_64` inverts the *intended*
function and therefore does NOT invert `mm_hash64`. Using the wrong inverse
would corrupt ref-delta decoding.

**Golomb-Rice on sorted hash deltas** (in `compress.rs`). Uniform hashes have
geometrically-distributed gaps — exactly the distribution Golomb-Rice is
optimal for. This reaches near-entropy compaction without a secondary
compression pass, which matters because DEFLATE/gzip adds latency with no
gain on already-random bitstreams. zstd wrapping then handles the repetitive
metadata cheaply.

**FxHashMap serialized as a sequence, not a map** (custom serde in `types.rs`).
The comment records a "magnitude" speedup vs. the default map encoding. The
visitor reconstructs the map from `(Kmer, u32)` pairs.

**Two independent two-stage formats** serve different use cases. `.sylref` +
`.sylspr` dereplicates k-mers across genomes (each k-mer owned by one genome
or the shared pool) — maximal compression, requires a curated taxonomy-aware
build. `.syl2db` keeps per-genome complete k-mer sets with no dereplication —
simpler build, usable with any existing `.syldb`, trades some disk efficiency
for build simplicity.

**Strains are stored contiguously in `.sylref`**, species representatives
first. This means a sample's hit genome ids cluster into small ranges, making
varint-delta encoding of genome ids efficient in `.sylspr`.

## Entry Points

| Task | Start at |
|---|---|
| Building a genome database | `sketch.rs::sketch` → `sketch_genome` |
| Sketching reads | `sketch.rs::sketch` → `write_read_sketch_file` |
| Querying / profiling | `contain.rs::contain` |
| Two-stage profiling path | `contain.rs::contain` → `TwoStageDb` or `RefIndex` screen then dense decode |
| Coverage/ANI statistics | `inference.rs::mme_lambda`, `estimate_lambda` |
| k-mer extraction kernel | `seeding.rs::fmh_seeds` (scalar) or `avx2_seeding.rs::extract_markers_avx2` (SIMD) |
| Compressed I/O (SYLZ format) | `compress.rs` |
| Reference-delta build | `refdelta/ref_build.rs::run_ref_build` |
| Reference-delta compress/decompress | `refdelta/sketch_compress.rs::run_ref_compress` |
| Two-stage DB conversion | `twostage_db.rs::run_db_convert` |
| Adding a new output format | Follow the pattern in `compress.rs` (magic, version byte, zstd frame) |

## Out of Scope

This file should be updated when modules are added, removed, or restructured,
but not for routine function-level changes — it documents shape, not implementation detail.
