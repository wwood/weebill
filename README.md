<p align="center">
  <img src="https://raw.githubusercontent.com/wwood/weebill/main/weebill.png" alt="weebill logo" width="50%" />
</p>

**Weebill** is a fork of [sylph](https://github.com/bluenote-1577/sylph), the fast and precise
species-level metagenomic profiler (ANI querying + taxonomic profiling). See the
[sylph repository](https://github.com/bluenote-1577/sylph) and
[sylph documentation](https://sylph-docs.github.io/) for the underlying method. Weebill maintains compatibility with sylph's command-line interface and sketch formats where possible. 

Weebill is currently in development and is experimental. Efforts are made to contribute improvements back upstream to sylph, but some features are beyond its scope.

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

## Installation

The quickest way is the installer script attached to each
[release](https://github.com/wwood/weebill/releases), which drops a prebuilt static Linux binary
(x86_64) into `~/.cargo/bin`:

```sh
curl --proto '=https' --tlsv1.2 -LsSf https://github.com/wwood/weebill/releases/latest/download/weebill-installer.sh | sh
```

Otherwise, build from source. This requires the [Rust toolchain](https://www.rust-lang.org/tools/install)
(`cargo`); install the `weebill` binary straight from GitHub:

```sh
cargo install --git https://github.com/wwood/weebill
```

This builds and installs `weebill` into `~/.cargo/bin` (make sure it is on your `PATH`). To build
from a local clone instead:

```sh
git clone https://github.com/wwood/weebill
cd weebill
cargo install --path .
```

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

### Releasing

Write user-facing notes under `## Unreleased` in [CHANGELOG.md](CHANGELOG.md) as you go, then:

```sh
python scripts/release.py --version X.Y.Z
```

This rolls the `Unreleased` notes into a `## Version X.Y.Z` section, bumps the version in
`Cargo.toml`, runs the tests and `dist plan`, commits, tags `vX.Y.Z`, publishes to crates.io and
pushes. The tag triggers `.github/workflows/release.yml` (generated by
[cargo-dist](https://opensource.axo.dev/cargo-dist/) — regenerate it with `dist generate` after
changing the `[workspace.metadata.dist]` config rather than editing it by hand), which builds the
binaries and attaches them plus the installer to the GitHub release.

## Citation

Jim Shaw and Yun William Yu. Rapid species-level metagenome profiling and containment estimation with sylph (2024). Nature Biotechnology.

A manuscript describing the weebill modifications to sylph might appear in future.
