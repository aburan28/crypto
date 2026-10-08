# Releases

Every merge to `main` publishes a GitHub release. There is no tagging step
and no manual trigger to remember: `.github/workflows/release.yml` runs on
`push: main` (and on `workflow_dispatch` if you need to re-cut one).

## What gets built

| Artifact | Contents / build | Runner |
| --- | --- | --- |
| `crypto-<tag>-x86_64-unknown-linux-gnu.tar.gz` | `crypto`, `ic`, `ecbench`, `curve_cover_check`, `hyperelliptic-cover`; `cargo build --release --bins` | ubuntu |
| `crypto-<tag>-x86_64-unknown-linux-musl.tar.gz` | same five executables, static via `musl-tools` | ubuntu |
| `crypto-<tag>-aarch64-apple-darwin.tar.gz` | same five executables | macos |
| `crypto-<tag>-x86_64-pc-windows-msvc.zip` | same five executables, with `.exe` suffixes | windows |
| `ecc2k130-cpu-<tag>-x86_64-linux.tar.gz` | `make -C ecc2k130 cpu` | ubuntu |
| `ecc2k130-cuda-<tag>-x86_64-linux.tar.gz` | `make -C ecc2k130 gpu` | ubuntu + nvcc |

Each release also carries `SHA256SUMS`; verify with `sha256sum -c SHA256SUMS`.

## Using the packaged cover tools

Download the archive for your platform from the
[latest release](https://github.com/aburan28/crypto/releases/latest), verify
it against `SHA256SUMS`, and unpack it. The cover tools are ordinary
executables; neither command requires `cargo run` or a Rust toolchain:

```sh
./hyperelliptic-cover --help
./curve_cover_check --help
```

Both catalog tools require caller-supplied data. `curve_cover_check` reads and
writes catalog files supplied with `--registry PATH` and `--output PATH`, and
`hyperelliptic-cover catalog-export` requires `--registry PATH` (with optional
`--output PATH`). Registries and generated catalogs are not bundled in the
executable archive. The other `hyperelliptic-cover` subcommands accept an
explicit supported model shape and perform no registry lookup. Consult each
command's `--help` for the interface in that release. The release job starts
both tools on every native target before it archives them.

## Versioning

The tag is `v<Cargo.toml version>-main.<run number>`, e.g. `v0.1.0-main.7`.
The run number makes every merge unique without needing a version bump, and
because these are ordinary (not pre-) releases, `/releases/latest` always
resolves to the newest `main` build. Bumping `version` in `Cargo.toml`
changes the prefix; nothing else needs editing.

## Two things worth knowing

**The native client is built `-march=x86-64-v3`, not `native`.**
`ecc2k130/Makefile` defaults `MARCH` to `native`, which bakes the build
host's CPU into the binary and SIGILLs on anything older. The release job
overrides it to the AVX2 + BMI2 + FMA baseline (Haswell / Zen 1 onwards).
For a machine you control, building locally with the default `MARCH=native`
will be faster — the release binary is the portable one.

**The CUDA binary is compiled, not tested.** GitHub runners have no GPU, so
the `cuda` job proves the device code compiles for every architecture in the
Makefile's gencode list (sm_80, sm_86, sm_89, sm_90, sm_120) and no more. It
needs CUDA >= 12.8 because `compute_120` does; the job installs `cuda-nvcc`
and `cuda-cudart-dev` only, about 330 MB rather than the full toolkit.
That job pins `runs-on: ubuntu-24.04` rather than `ubuntu-latest`: the
NVIDIA apt repository is addressed by distro and the CUDA version is
pinned, so a floating runner image would eventually 404 the install and
block every release. Bump the two together or not at all.

## Coverage this added

The `native` job also runs `make test` in `gpu/ecc`, `gpu/ecc2k`,
`gpu/semaev`, `gpu/macaulay` and `gpu/btcpuzzle`, which build the `.cu` and
`.cuh` sources as host C++ and run their self-checks. No workflow compiled
those before.

## If a release fails

Pull requests touching `release.yml`, `Cargo.toml`, Rust sources under `src/`
or `ecc2k130/Makefile` build every artifact as a dry run and stop short of
publishing, so the macOS and Windows legs are proven before a merge can turn a
release red.

The `release` job requires all four build jobs, so a failure publishes
nothing rather than a partial set. The matrix is `fail-fast: false`, so one
platform breaking still tells you about the other three in the same run.
