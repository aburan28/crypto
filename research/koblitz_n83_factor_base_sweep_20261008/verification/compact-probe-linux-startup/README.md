# Compact-probe Linux startup check, 2026-10-09

This is build and startup evidence for source commit
`2742a09899298fc1a11802a3dbc1e83faa01df20`. It does not import a
retained factor base, build the compact index, search for a relation, measure
rank or supply an index-calculus runtime. The one-hour retained-base pilot
remains exhausted. The Linux cross-build took 3m36s of separate development
compute; the two empty/invalid-panel container checks took 447ms and 383ms
of supervisor process wall. Neither interval is a factor-base timing.

The build used the cached Zig 0.16.0 toolchain and the repo's pinned stable
Rust compiler. The exact command was:

```sh
env TMPDIR=/private/tmp CARGO_BUILD_JOBS=4 CARGO_PROFILE_RELEASE_CODEGEN_UNITS=256 CARGO_TARGET_DIR=/private/tmp/n83-linux-target RUSTC=/Users/adamburan/.rustup/toolchains/stable-aarch64-apple-darwin/bin/rustc CARGO_ZIGBUILD_CACHE_DIR=/Volumes/SSD990/llm/tmp/n83-cargo-zigbuild-cache ZIG_LOCAL_CACHE_DIR=/Volumes/SSD990/llm/tmp/n83-zig-cache ZIG_GLOBAL_CACHE_DIR=/Volumes/SSD990/llm/tmp/n83-zig-cache cargo zigbuild --release --no-default-features --target x86_64-unknown-linux-musl --example koblitz_n83_factor_base_export --locked --offline
```

`build.log` ends with a successful release build. The executable was a
statically linked x86-64 Linux ELF, 4.1 MiB on this host, with SHA-256
`2faa775afabdce988dbbf179a721ee746aadf92a3b6ecbede56c96dba9a502d8`.
The locally cached container image resolved to
`sha256:0dd364ba7e10242f07755449e3a3d0e35f9efd987952737b90def6709ab0c5ce`;
Docker did not pull an image. Both supervisor calls used 30-second worker
wall caps, 128 MiB memory, zero swap, one CPU, no network and the exact
unordered K=64 state cap of 172,640. The source snapshot contained only
the exporter, compact index and primary adapter Rust files plus the
attestation. `config.json` retains their SHA-256 digests.

The first `empty-panel/` run reached a generic missing-file error. That
error by itself did not prove which startup stage failed, so the second
`invalid-panel/` run supplied two public files containing `{}`. Its stderr
reports `panel or replay receipt is not complete and bound`. In the worker,
that check follows cgroup verification and compiled-source comparison; the
specific error therefore establishes that both gates returned successfully
before the deliberately invalid panel was rejected. Each outer receipt
correctly records `PRODUCER_FAILURE_worker_exit`, exit code 1, no worker
receipt, unchanged frozen inputs, and null rank and total-runtime fields.
These are expected negative controls, not successful probes.

Run `python3 research/koblitz_n83_factor_base_sweep_20261008/verification/compact-probe-linux-startup/verify.py`
from the checkout to replay the stored receipt hashes, exact synthetic
inputs, source-file hashes from the frozen Git commit and the retained
release-library, exporter-example, study-Python and boundary-Python test
logs. Those suites passed 2,249, 18, 15 and 16 tests respectively. The verifier
returns `PASS_receipts` for the two expected failures. It cannot certify a
future retained-base index build, relation yield, memory peak or wall time;
those require a separately authorized run with its own immutable receipts.
