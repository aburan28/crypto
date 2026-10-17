# Releases

Every push to `main` starts the release workflow. A GitHub release is
published only after every required build job succeeds. The workflow can also
be started with `workflow_dispatch`; publication still requires a `main` ref.
This page describes archives produced by the workflow revision in the same
commit. It is not evidence that an older or currently latest release already
contains newly configured files or executables; bind package claims to the
asset's source commit and successful Actions run.

## What gets built

| Artifact | Contents / build | Runner |
| --- | --- | --- |
| `crypto-<tag>-x86_64-unknown-linux-gnu.tar.gz` | `crypto`, `ic`, `ecbench`, `curve_cover_check`, `curve_standards`, `hyperelliptic-cover`, `ecdlp-audit`; selected `cargo build --locked --release --bin …` targets | ubuntu |
| `crypto-<tag>-aarch64-apple-darwin.tar.gz` | same seven executables | macos |
| `crypto-<tag>-x86_64-pc-windows-msvc.zip` | same seven executables, with `.exe` suffixes | windows |
| `isogeny-<tag>-x86_64-unknown-linux-gnu.tar.gz` | dedicated `isogeny` executable, README, and build receipt | ubuntu |
| `isogeny-<tag>-x86_64-unknown-linux-musl.tar.gz` | same dedicated executable and records, statically linked | ubuntu |
| `isogeny-<tag>-aarch64-apple-darwin.tar.gz` | same dedicated executable and records | macos |
| `isogeny-<tag>-x86_64-pc-windows-msvc.zip` | same dedicated executable and records, with `.exe` suffix | windows |
| `ecc2k130-cpu-<tag>-x86_64-linux.tar.gz` | `make -C ecc2k130 cpu` | ubuntu |
| `ecc2k130-cuda-<tag>-x86_64-linux.tar.gz` | `make -C ecc2k130 gpu` | ubuntu + nvcc |

A release produced from this workflow revision also carries `SHA256SUMS`;
verify with `sha256sum -c SHA256SUMS`. Its three Rust archives contain a staged
`docs/curves` catalog snapshot,
the standalone hyperelliptic command manual and schemas, the catalog-wide
transfer-assessment ledger, the fixed-F8 order certificate and its
implementation-separated native replay command, and non-authoritative
`docs/ecdlp-audit` example records. The main Rust bundle does not currently
have a musl artifact; the dedicated
`isogeny` tool retains its separate musl build.

## Using the packaged cover tools

Once a matching workflow revision has been published, download the archive for
your platform from its release page (the
[latest release](https://github.com/aburan28/crypto/releases/latest) may still
refer to an earlier commit), verify it against `SHA256SUMS`, and unpack it.
The cover tools are ordinary
executables; neither command requires `cargo run` or a Rust toolchain:

```sh
./hyperelliptic-cover --help
./curve_cover_check --help
```

The archive carries the release's `docs/curves` registry, cover report, linkage
graph, standards registry, hyperelliptic command contract, evidence schemas,
and transfer-assessment ledger. Catalog commands still take paths:
`curve_cover_check` accepts `--registry PATH` and `--output PATH`, and
`hyperelliptic-cover catalog-export` requires `--registry PATH` (with an
optional `--output PATH`). Run from the extracted archive root to use relative
`docs/curves/...` paths, or pass explicit paths. The other
`hyperelliptic-cover` subcommands accept an explicit supported model shape and
perform no registry lookup. See
`docs/curves/HYPERELLIPTIC_INFRASTRUCTURE.md` for the exact request contract,
infinity-model boundary, and structured unsupported cases.

Before compression, the release job enters the staged archive tree, requires
byte-identical native regeneration of the saved cover catalog and link graph,
and runs one construction, verification, binary-Jacobian identity operation,
and full catalog export against the staged registry on every native target.
It also replays the fixed-F8 order certificate with the implementation-separated
native verifier embedded in `ecdlp-audit`. The cover CI runs that verifier's
full 45-test suite, including the fixed-F8 verifier's 12 replay-and-fail-closed
mutation tests, then executes both the repository-bound plan replay and the
fixed-F8 certificate replay.
That certificate establishes only `#Jac(H)(GF(8)) = 486`; it does not certify
the norm map or its arithmetic. These are installed-layout smoke tests, not
general mathematical validation. The assessment ledger retains
repository-relative provenance paths; replaying all of those hashes and
reports still requires the matching full checkout.

## Using the packaged ECDLP audit CLI

Releases produced from this workflow revision are configured to contain
`ecdlp-audit` and a copy of
`docs/ecdlp-audit` whose records remain explicitly non-authoritative. The
portable audit checks include process startup and the bounded fixed-F8 order
certificate replay:

```sh
./ecdlp-audit --help
./ecdlp-audit verify-fixed-f8-order \
  --input docs/curves/fixed-f8-genus3-order-certificate.json --repo-root .
```

The archived planning mock is an example, not a self-contained validation
fixture. `ecdlp-audit verify-plan` checks repository-bound source, report, and
deliverable paths and pinned Git objects. Run that command from the matching
full repository checkout with `--repo-root`; missing archive paths must not be
interpreted as a mathematical, implementation, security, or completeness
finding. Workflow configuration alone does not establish that a cross-platform
release run or published package has passed; use the Actions run and release
assets for the relevant commit as evidence.

## Versioning

The workflow reads `major.minor.patch` from `Cargo.toml`, adds the GitHub run
number to `patch`, and tags the built commit as `v<major>.<minor>.<sum>`. For
example, crate version `0.1.0` in run 100 produces `v0.1.100`. Because these
are ordinary (not pre-) releases, `/releases/latest` resolves to the newest
successfully published `main` build. Bumping `version` changes that base.

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

Pull requests touching `release.yml`, its scope classifier, `Cargo.toml`, Rust
sources under `src/`, `ecc2k130/Makefile`, or the staged `docs/curves` and
`docs/ecdlp-audit` payloads, READMEs, or licenses build the affected artifacts
as a dry run and stop short of publishing. Merge-queue runs use the same
local-diff classifier, so records-only package changes cannot silently skip
the corresponding archive lane.

The `release` job requires all four build jobs, so a failure publishes
nothing rather than a partial set. The matrix is `fail-fast: false`, so one
platform breaking still tells you about the other three in the same run.
