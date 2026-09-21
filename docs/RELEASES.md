# Releases

Every merge to `main` publishes a GitHub release. There is no tagging step
and no manual trigger to remember: `.github/workflows/release.yml` runs on
`push: main` (and on `workflow_dispatch` if you need to re-cut one).

## What gets built

| Artifact | Built by | Runner |
| --- | --- | --- |
| `crypto-<tag>-x86_64-unknown-linux-gnu.tar.gz` | `cargo build --release --bins` | ubuntu |
| `crypto-<tag>-x86_64-unknown-linux-musl.tar.gz` | same, static via `musl-tools` | ubuntu |
| `crypto-<tag>-aarch64-apple-darwin.tar.gz` | same | macos |
| `crypto-<tag>-x86_64-pc-windows-msvc.zip` | same | windows |
| `ecc2k130-cpu-<tag>-x86_64-linux.tar.gz` | `make -C ecc2k130 cpu` | ubuntu |
| `ecc2k130-cuda-<tag>-x86_64-linux.tar.gz` | `make -C ecc2k130 gpu` | ubuntu + nvcc |

Each release also carries `SHA256SUMS`; verify with `sha256sum -c SHA256SUMS`.

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

Pull requests touching `release.yml`, `Cargo.toml` or `ecc2k130/Makefile`
build every artifact as a dry run and stop short of publishing, so the
macOS and Windows legs are proven before a merge can turn a release red.

The `release` job requires all four build jobs, so a failure publishes
nothing rather than a partial set. The matrix is `fail-fast: false`, so one
platform breaking still tells you about the other three in the same run.
