<p align="center">
  <img src="https://raw.githubusercontent.com/wwood/weebill/main/weebill.png" alt="weebill logo" width="50%" />
</p>

**Weebill** is a fork of [sylph](https://github.com/bluenote-1577/sylph), the fast and precise
species-level metagenomic profiler (ANI querying + taxonomic profiling). See the
[sylph repository](https://github.com/bluenote-1577/sylph) and
[sylph documentation](https://sylph-docs.github.io/) for the underlying method. Weebill maintains compatibility with sylph's command-line interface and sketch formats where possible. 

Weebill is currently in development and is experimental. Efforts are made to contribute improvements back upstream to sylph, but some features are beyond its scope.

## Contents

- [Installation](#installation)
- [Usage examples](#usage-examples)
  - [From reads to a profile (two-stage)](#from-reads-to-a-profile-two-stage)
  - [Pooling samples with `profile --merge`](#pooling-samples-with-profile---merge)
  - [Reference-delta samples with `profile --reference`](#reference-delta-samples-with-profile---reference)
- [Changes in the weebill fork](#changes-in-the-weebill-fork)
- [Minor changes since the fork](#minor-changes-since-the-fork)
- [Development](#development)
- [Citation](#citation)

## Installation

Prebuilt static Linux binaries (x86_64) are attached to each release. Download one from the
[**releases page**](https://github.com/wwood/weebill/releases), then place the `weebill`
binary somewhere on your `PATH` (e.g. `~/.cargo/bin` or `~/.local/bin`) and make it executable:

```sh
chmod +x weebill
```

Alternatively, build from source with the [Rust toolchain](https://www.rust-lang.org/tools/install)
(`cargo`). Install the `weebill` binary straight from GitHub:

```sh
cargo install weebill
```

This builds and installs `weebill` into `~/.cargo/bin` (make sure it is on your `PATH`). To build
from a local clone instead:

```sh
git clone https://github.com/wwood/weebill
cd weebill
cargo install --path .
```

## Usage examples

Weebill has a two-step model: **sketch** sequences into compact indexes, then **profile** (or
**query**) those indexes. Reads become *sample* sketches; genomes become a *database* sketch.
The examples below assume `weebill` is on your `PATH`; add `-t <threads>` to any command to
parallelise.

### From reads to a profile (two-stage)

This is the recommended everyday workflow. Sketch the genomes once, convert the database into a
two-stage seekable database (`.syl2db`), sketch your reads into compressed sample sketches, then
profile with `--two-stage` — which is much faster and uses far less RAM than a standard profile
(see [Changes in the weebill fork](#changes-in-the-weebill-fork)).

```sh
# 1. Sketch a genome database (standard .syldb; db-convert reads this format)
weebill sketch -g genomes/*.fa -o gtdb

# 2. Convert it into a two-stage seekable database (-> gtdb.syl2db)
weebill db-convert gtdb.syldb -o gtdb

# 3. Sketch metagenome reads into compressed sample sketches (-> sketches/*.sylspc)
#    single-end (one .sylspc per input file):
weebill sketch -r sampleA.fastq.gz -r sampleB.fastq.gz --compressed-database sketches/
#    paired-end:
weebill sketch -1 sampleA_1.fq.gz -2 sampleA_2.fq.gz --compressed-database sketches/

# 4. Two-stage taxonomic profile: relative abundance + ANI per detected species (TSV to stdout)
weebill profile --two-stage gtdb.syl2db sketches/*.sylspc > profile.tsv
```

Compressed sample sketches are smaller on disk and are read transparently by `profile`,
`query`, and `inspect`. `profile` also accepts raw fastq/fasta directly and will sketch them on
the fly, e.g. `weebill profile --two-stage gtdb.syl2db -r sampleA.fastq.gz`. Swap `profile` for
`query` (against `gtdb.syldb`) to get nearest-neighbour containment ANI instead of a profile.

### Pooling samples with `profile --merge`

`--merge` combines several read inputs — any mix of pre-sketched samples and raw reads — into a
single sample (summing per-genome k-mer counts) and profiles that one pooled sample instead of
each input separately. Handy for co-assemblies or for pooling multiple sequencing runs of the
same biological sample.

```sh
weebill profile --two-stage gtdb.syl2db --merge -S patient1_pooled \
    sketches/run1.sylspc sketches/run2.sylspc sketches/run3.sylspc \
    -o patient1_pooled.profile.tsv
```

`-S` names the merged sample in the `Sample_file` column (default: `merged`). The same
`--merge`/`-S` flags exist on `weebill sketch` if you want to write the pooled sketch to disk
rather than profile it immediately, and on the standalone `weebill merge` sub-command.

### Reference-delta samples with `profile --reference`

Reference-delta compression stores each sample as only the reference hashes it contains plus any
novel hashes, producing very small `.sylspr` files (see
[Changes in the weebill fork](#changes-in-the-weebill-fork)). Build the reference once, compress
your samples against it, then profile the `.sylspr` samples by passing the same reference so they
can be decoded.

```sh
# 1. Build a reference database once from the genome database sketch (-> gtdb.sylref)
weebill ref-build gtdb.syldb -o gtdb

# 2. Sketch reads straight to reference-delta samples (-> sketches_ref/*.sylspr)
weebill sketch -r sampleA.fastq.gz --reference gtdb.sylref -d sketches_ref/
#    (or compress the sample sketches made above: weebill ref-compress -r gtdb.sylref sketches/*.sylspc -d sketches_ref/)

# 3. Profile the .sylspr samples — --reference is REQUIRED to decode them
weebill profile gtdb.syldb sketches_ref/*.sylspr --reference gtdb.sylref > profile.tsv
```

`--reference` writes `.sylspr` directly and so cannot be combined with `--compressed-database`;
give the output directory with `-d`. The profile produced is identical to profiling the
uncompressed samples — reference-delta compression is lossless.

## Changes in the weebill fork

The binary is installed as `weebill`. Weebill changes:

- **Lighter weight profiling** - the `profile` command can use a "2 stage" profiling approach which is
~7–27× faster and uses ~6× less RAM (a flat ~4.4 GB vs the whole ~26 GB database) than standard sylph profiling (when input is sketches - FASTA/FASTQ inputs are also faster but more modestly). The on-disk database is also ~22% smaller. The profiles produced are effectively identical to standard sylph profiling, and the speed/RAM boost means that choosing smaller `c` values is more computationally feasible. To use this mode, see `weebill db-convert` and `weebill profile --two-stage`.
- **Compressed sketches** — `weebill sketch --compressed-output`/`--compressed-database` write
  `.sylspc` samples and `.syldbc` databases (~55% smaller samples, ~30%+ smaller databases). Hashes
  are sorted, delta-encoded and Golomb–Rice coded, then wrapped in a zstd frame. `query`, `profile`,
  and `inspect` read the legacy and compressed formats interchangeably (detected by content, not
  extension).
- **Reference-delta compression** — `weebill ref-build` dereplicates database sketches into a
  two-stage seekable reference (`.sylref`), and `weebill ref-compress` encodes a sample
  (`.sylsp` → `.sylspr`) by storing only which reference hashes are present plus any novel hashes
  (~60% smaller than `.sylspc`, ~80% smaller than bincode, losslessly). `query`/`profile` can read
  `.sylspr` samples directly via `--reference <ref.sylref>`. `ref-build` is streaming and parallel
  with RAM bounded by `--max-ram`, and the two-stage `.sylref` loads only the genome blocks a sample
  needs. `ref-compress` also supports `--decompress`, `inspect`, and `verify` modes.
- **Error-k-mer encoding** — `ref-build --store-genomes` additionally stores each species
  representative's nucleotide sequence (2-bit packed, ~0.25 byte/bp) in the `.sylref`. `ref-compress`
  then recognises sample hashes that are a *single-base substitution* of a reference k-mer — the
  dominant kind of sequencing-error hash — and stores them as compact `(genome position, k-mer
  offset, replacement base)` triples (positions sorted and Golomb–Rice delta-coded, the three fields
  in separate arrays for the zstd wrapper) instead of full-price "novel" hashes. The genome sequence
  itself acts as the perfect-hash fingerprint, so matches are reconstructed exactly and the round trip
  stays lossless. On 100× Illumina-like reads at 0.5% error this roughly halves the `.sylspr`
  (~91% of non-reference hashes are single-error k-mers); the cost is the one-time per-reference
  genome storage, amortised across all samples compressed against it.
- **`weebill merge`** — a new sub-command that merges several sample sketches into a single sketch
  (summing per-genome k-mer counts). It reads any mix of compressed (`.sylspc`) and
  reference-compressed (`.sylspr`, via `--reference`) inputs, and writes the result as compressed
  `.sylspc` by default (or `.sylspr` with `--ref-compress`). Legacy uncompressed `.sylsp` samples
  record no read count and so cannot be merged — re-sketch them with `--compressed-database` first.

## Minor changes since the fork

Smaller improvements and fixes beyond the headline features above:

- **Multi-threaded read sketching** — sketching a read input now uses threads at two levels. Across
  inputs, single-end, paired-end and interleaved passes are drained concurrently rather than one at a
  time. Within a single input, a dedicated reader thread does the IO/decompression while rayon workers
  extract k-mers from batches of reads in parallel; only the order-dependent dedup fold stays serial,
  so the sketch is byte-for-byte identical to the single-threaded result (and independent of `-t`).
  This lets one large read file scale across cores (~2× on 4 threads for a single-end input) where it
  previously ran the k-mer work on one core regardless of `-t`.
- **`profile --apply-unknown`** — converts a profile produced *without* `-u`/`--estimate-unknown`
  into the profile `-u` would have produced, without re-profiling. `-u`'s effect is a per-sample
  rescale of just two columns (`Eff_cov` → `True_cov`, and `Sequence_abundance` by the estimated
  covered-bases fraction; `Taxonomic_abundance` and the ANIs are unchanged), and those scalars are
  recoverable from the sample sketch (k-mer-count distribution + read length) plus the database
  (genome sizes). Pass the original TSV, the same pre-sketched sample(s) and the original database
  (both required): `weebill profile --apply-unknown profile.tsv gtdb.syl2db sample.sylspr --reference gtdb.sylref`.
  The result matches a real `-u` run to the precision of the input TSV, but is **not guaranteed to be
  bit-for-bit identical**: `--apply-unknown` can only read the already-rounded columns printed in the
  TSV (e.g. `Eff_cov` at 3 decimal places), whereas a real `-u` run rescales the full-precision internal
  values, so individual `True_cov`/`Sequence_abundance` cells may differ by a unit in the last printed
  place (the `Sequence_abundance` scalar can drift slightly more, as it sums the rounded per-row coverage).
- **`sketch`/`profile --merge`** — a `--merge` flag on `sketch` and `profile` (distinct from the
  `merge` sub-command) combines all read inputs into one sample sketch; `-S` names it. Used in the
  [pooling example](#pooling-samples-with-profile---merge) above.
- **`--tolerate-empty-inputs`** — treat a read input containing zero reads (an empty file or an SRA
  FIFO with no unpaired reads) as a valid zero-read sketch instead of an error, so a single empty
  stream does not abort a `--merge`. A run where *every* input is empty is still an error.
- **AVX-512 k-mer extraction** — an 8-lane AVX-512 seeding path (with `vpcompressq`) alongside the
  existing AVX2 path, for faster k-mer extraction on capable CPUs.
- **Reproducible sketching** — `rand` is pinned and the read-deduplication cuckoo filter RNG is
  seeded, so repeated sketching of the same input is deterministic.
- **More accurate paired read lengths** — the read-length estimate for paired-end data uses the
  mean of both mates and excludes sub-*k* mates, improving `--estimate-unknown` coverage scaling.
- **Corruption detection in the new sketch/database formats** — the compressed sketches are
  checksummed by their zstd frame and the seekable databases carry a whole-file checksum, both
  verified on read (e.g. via `weebill inspect`).

## Development

Formatting and tests are checked in CI (`cargo fmt --check` and `cargo test`) on every push and
pull request. Clippy also runs, but is advisory: the code inherited from sylph has pre-existing
findings, so its output is reported without failing the build.

A versioned pre-commit hook is provided to reject commits with unformatted Rust code. After cloning, enable it once with:

```sh
cd weebill
git config core.hooksPath .githooks
cargo test
```

The hook runs `cargo fmt --all -- --check`; run `cargo fmt` to fix formatting before committing.

## Citation

Jim Shaw and Yun William Yu. Rapid species-level metagenome profiling and containment estimation with sylph (2024). Nature Biotechnology.

A manuscript describing the weebill modifications to sylph might appear in future.
