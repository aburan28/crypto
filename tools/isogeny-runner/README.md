# Local P192/P224 isogeny runs

This standalone Cargo package builds the existing construction and process-supervision
sources from the [2026-10-09 study](../../research/p192_p224_large_degree_isogenies_20261009/).
The frozen sources and recorded evidence are unchanged. It does not build the repository's
root library. The independent replay package also builds separately from that library.

The kernel constructor currently supports exactly two public, fixed-seed cases:

| Source | Prime degree | Extension used for kernel construction | Output coverage |
| --- | ---: | --- | --- |
| P192 | 10453 | Fp^3 | One Frobenius eigenline of two |
| P224 | 1471 | Fp^5 | One Frobenius eigenline of two |

Each output carries `status: PASS` for its checks and `degree_coverage: PARTIAL`.
The unresolved second eigenline remains outside these construction results. Arbitrary
degrees are supported by the separate `isogeny-algos isogenies` command, through its
modular-polynomial route; they are not supported by this fixed-case kernel constructor.

## Build

Run these commands from the repository checkout containing this branch. Rust 1.93.1
was used for validation. The supervisor supports macOS and Linux.

```sh
cargo build --release --locked \
  --manifest-path tools/isogeny-runner/Cargo.toml \
  --target-dir target/isogeny-local

cargo build --release --locked \
  --manifest-path research/p192_p224_large_degree_isogenies_20261009/Cargo.toml \
  --target-dir target/isogeny-local --bin replay-one
```

The replay package has ordinary Cargo dependencies; its first build needs either access
to the Cargo registry or a populated local Cargo cache. The constructor depends only
on the repository's dependency-free `isogeny_algos` crate.

## Construct and independently verify P224, degree 1471

Use a fresh output directory for each run. This example creates one in the system's
temporary directory and prints its path so the JSON records can be inspected afterward.
Move it to durable storage if you want to retain the run.

```sh
ISOGENY_RUN_OUT=$(mktemp -d "${TMPDIR:-/tmp}/isogeny-p224-1471.XXXXXX")
printf '%s\n' "$ISOGENY_RUN_OUT"

target/isogeny-local/release/isogeny-probe \
  target/isogeny-local/release/isogeny-kernel \
  "$ISOGENY_RUN_OUT" p224 1471 1800 "$(git rev-parse HEAD)" --single-map

target/isogeny-local/release/replay-one \
  "$ISOGENY_RUN_OUT" docs/curves/registry.json
```

For P192, create another output directory and replace `p224 1471` with `p192 10453`
in the supervisor command. Keep `--single-map`: it records the constructor's actual
one-map coverage. Independent verification of degree 10453 took approximately 11 minutes
on the study's Apple M4 Pro host; other hosts can differ.

The supervisor's positional arguments are:

```text
isogeny-probe CONSTRUCTOR OUTPUT_DIR CURVE ELL TIMEOUT_SECONDS SOURCE_COMMIT [--single-map]
```

It supplies the published source order and seed 1 automatically, enforces the construction
time limit and an 8 GiB resident-memory cap, and records the executable hash, commands,
output hashes, exit status, elapsed time, and sampled peak memory. Independent replay
runs afterward and is outside this construction timeout. The supervisor can finish
normally after recording a timeout or failed construction: inspect the receipt's `status`;
its own process exit status is not a construction certificate. Replay requires a passing
construction receipt and writes its certificate only after all checks pass.

Files in a P224 run are:

- `search.json`: search summary and protocol hash.
- `p224/ell-1471.json`: explicit kernel, codomain and rational x-map.
- `p224/ell-1471.stderr.txt`: construction stage messages.
- `p224/ell-1471.receipt.json`: bounded construction outcome and hashes.
- `replay.json`: independent verification certificate.
- `curves.json`: independently checked destination curve records.

The replay checks kernel validity, the exact rational-map curve identity, codomain
consistency, 20 fresh scalar-transport cases, and preservation of the public subgroup.
Replay outputs do not automatically update the canonical registry.

## Other prime degrees using the general CLI

Build the standalone CLI into the same output directory:

```sh
cargo build --release --locked --manifest-path isogeny_algos/Cargo.toml \
  --target-dir target/isogeny-local --bin isogeny-algos

target/isogeny-local/release/isogeny-algos version
```

For example, run a bounded P192 degree-509 construction attempt:

```sh
ISOGENY_RUN_OUT=$(mktemp -d "${TMPDIR:-/tmp}/isogeny-p192-509.XXXXXX")
target/isogeny-local/release/isogeny-probe \
  target/isogeny-local/release/isogeny-algos \
  "$ISOGENY_RUN_OUT" p192 509 600 "$(git rev-parse HEAD)"
```

Omit `--single-map` for this CLI: the supervisor expects the complete pair of split
eigenline maps. Choose a split prime from the study's screening records. The general
modular-polynomial route timed out at this degree in the recorded 600-second run;
this command reproduces a bounded attempt rather than a guaranteed successful map.
See [the CLI guide](../../isogeny_algos/docs/USAGE.md) for direct commands, arbitrary
curve parameters, order certification, and JSON status semantics.

## Checks

```sh
cargo test --release --locked --manifest-path tools/isogeny-runner/Cargo.toml \
  --target-dir target/isogeny-local

cargo test --release --locked --manifest-path isogeny_algos/Cargo.toml \
  --target-dir target/isogeny-local --lib

cargo test --release --locked \
  --manifest-path research/p192_p224_large_degree_isogenies_20261009/Cargo.toml \
  --target-dir target/isogeny-local \
  --bin verification-poly-checks --bin verification-kernel-checks
```

These test the existing extension arithmetic, supervisor invariants, isogeny algorithms,
and independent polynomial and kernel verification routines. They do not weaken the
repository's separate root-library validation gate for publishing the branch.
