# Local Cargo entry point validation, 2026-10-09

The package adds build targets and instructions for the existing research sources.
It changes no construction algorithm, frozen input, recorded experiment, canonical
curve record, or verification routine. The underlying source revision is
`d30db4ebe3ea6aeb379b763f15cbbee6164207f5` on
`codex/p192-p224-large-degree-isogenies-20261008`.

Host: Apple M4 Pro, macOS 26.6, Rust/Cargo 1.93.1. Builds and tests used release
optimization. The package lockfile contains only this package and the repository's
`isogeny_algos` path dependency.

| Check | Observed result |
| --- | --- |
| Locked Cargo build of `isogeny-kernel` and `isogeny-probe` | PASS |
| `cargo test --release --locked --manifest-path tools/isogeny-runner/Cargo.toml --target-dir target/isogeny-local` | 4 extension tests and 2 supervisor tests passed |
| Supervised P224 degree 1471, seed 1, 1800 s limit | PASS; one map; no unavailable memory samples |
| Supervised P192 degree 10453, seed 1, 1800 s limit | PASS; one map; no unavailable memory samples |
| Byte comparison of both newly constructed map files against the sealed study files | Identical |
| Independent replay of the newly constructed P224 degree-1471 map | PASS; exact map identity, kernel, codomain, public subgroup, and 20 scalar-transport checks |
| Original study manifest after adding this package | All 209 bound files verified |
| Required root release-library gate | Still blocked by the separately recorded 643 baseline compiler errors; no root-library files changed here |

Construction stdout hashes, reproduced exactly:

- P224/1471: `ecf1109dbff2a72000fd84d9939127db68e2287e55e4ec10d143229874c230c8`.
- P192/10453: `aff763af2e9adb95765a7f29223aed7dc46c81f488301a3ba06db10de5595f29`.

The fresh checks ran in `/private/tmp/isogeny-local-validation.W2kf58`. The identical,
durable map files are in
[`kernel-p224-1471/p224/ell-1471.json`](../../research/p192_p224_large_degree_isogenies_20261009/kernel-p224-1471/p224/ell-1471.json)
and
[`kernel-p192-10453/p192/ell-10453.json`](../../research/p192_p224_large_degree_isogenies_20261009/kernel-p192-10453/p192/ell-10453.json).
The fresh independent replay used the current `replay-one` source built with its own
Cargo manifest into the existing `/private/tmp/p192-p224-continuation.f9sck4/replay`
build cache. The README uses a checkout-local target directory for convenience;
the target directory does not change the verifier source or inputs.

The P192 map's independent certificate remains the sealed study's
[`kernel-p192-10453/replay.json`](../../research/p192_p224_large_degree_isogenies_20261009/kernel-p192-10453/replay.json).
Its exact verification was not repeated for this packaging change; the new construction
output is byte-identical to that certificate's input. Coverage remains one of two
Frobenius eigenlines for each fixed case.

The existing standalone algorithm suite previously passed 114 tests, and the independent
polynomial/kernel routines previously passed 12 tests. Their sources remain unchanged.
The root gate's retained failure log is
[`validation/root-release-library.log`](../../research/p192_p224_large_degree_isogenies_20261009/validation/root-release-library.log).
These checks establish reproducibility of the new build entry points; no new timing
comparison or degree-enumeration claim is made.
