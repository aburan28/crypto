# Exact execution and reproduction

Native analysis and orchestration source: `examples/conductor_gap_ir.rs`.
The profiling executable is the archived ecbench, not a rebuilt solver.

```sh
export CARGO_HOME=/workspace/toolchains/cargo
export RUSTUP_HOME=/workspace/toolchains/rustup
cargo build --locked --release --example conductor_gap_ir -j3
cargo test --locked --release --example conductor_gap_ir -j3
cargo test --locked --release --lib -j3 -- --test-threads=3
```

Copy the original ecbench before any build overwrites it. Its SHA-256 must be
`429de171c022b71427759b47d411d8f8d6dd52e9f2daa819de23f4a67536194a`.
The recorded source commit is `7df956f97dad09c96dab12772b2b4cb244146f`.
Valgrind was extracted locally from Debian's `valgrind_3.24.0-3_amd64.deb` at
`https://deb.debian.org/debian/pool/main/v/valgrind/`; the package digest is in
`tool-hashes.sha256`. No system package or privileged host configuration changed.

```sh
study=research/ecbench_conductor_gap_ir_20261014
archive=research/ecbench_conductor_gap_20261012
runner=target/release/examples/conductor_gap_ir
ecbench=/workspace/toolchains/volcano-ir/ecbench
valgrind=/workspace/toolchains/volcano-ir/valgrind/usr/bin/valgrind
export VALGRIND_LIB=/workspace/toolchains/volcano-ir/valgrind/usr/libexec/valgrind
$runner freeze "$study" "$archive" "$ecbench"
# Commit protocol, native source, exact jobs and input hashes before profiling.
$runner run "$study" "$archive" "$ecbench" "$valgrind"
$runner analyze "$study" "$archive"
node "$study/render.mjs"
```

The initial preparation helper used the already-built dependency artifacts:

```sh
rustc --edition=2021 -C opt-level=2 -L dependency=target/release/deps \
  --extern crypto_lib=target/release/deps/libcrypto_lib-f2b277864209ea01.rlib \
  --extern serde_json=target/release/deps/libserde_json-d00c0fe8816b70e5.rlib \
  examples/conductor_gap_ir.rs \
  -o /workspace/toolchains/volcano-ir/conductor_gap_ir_prepare
```

Preparation generated all 100 frozen input files without executing a solve.
The full library-test compiler and queued Cargo example build are paused for
the profiles and resumed afterward. The same native helper can execute the
census while those processes are stopped. Required Cargo library and example
checks still complete before publication; no gate is removed.

`freeze` refuses an existing inputs directory. `run` refuses an existing raw
directory. Each job invokes `timeout --signal=TERM --kill-after=15 1800` around
Valgrind, preserving exit status, input, child output, stderr, complete raw
parts, parsed Ir receipt and SHA-256 receipt. Jobs are executed sequentially;
no builds or tests overlap profiling. A failure remains in analysis.json and
makes inference inadmissible. `analyze` verifies every raw checksum, re-parses
all parts, checks original point/scalar identities and deterministic counters.

For archive-only replay, extract the lossless raw archive to this study's
`raw/` directory in a fresh checkout, then run the native `analyze` command.
Regenerated analysis must compare byte-identically to the committed file.
Presentation code only lays out recorded native values; it does no inference.
