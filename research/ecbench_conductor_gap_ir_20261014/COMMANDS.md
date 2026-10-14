# Exact execution and reproduction

Native analysis and orchestration source: `examples/conductor_gap_ir.rs`.
Version 1 used the archived ecbench. The complete version 2 census uses the
newly built, hash-pinned executable described below; the archived executable
remains the source of the exact frozen inputs.

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
The first library-test compilation and queued example build were paused,
resumed, then terminated when the profiling fix required a new source snapshot.
Their logs are retained and supply no completed test result. Version 2 rebuilt
the profiler and native driver before starting the full census. Required final
library and example checks run after profiling; no validation gate is removed.

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

## Version 2 execution

The measurement-helper change is committed at
`a92902a93b1caa467156725f90401a9ac53ce2b4`, on base
`275beb93f3bf33ade6c19f19e8b336027a00b297`. Cargo's default release build
produced the profiler; its SHA-256 is
`ee7b82cf016f2433e99edd096f6a5dbfb3b3f9c4447f1583fae22864252a12e1`.
The configuration was committed before the census at the head recorded in
`preregistration-v2-commit.txt`. The native driver SHA-256 is
`841210a2d4a280f405f8b470e356c70d03b8072bef6d0a84430dc36a5e420bca`.
The driver sets `ECBENCH_CALLGRIND_SOLVE=1`,
`ECBENCH_CALLGRIND_PHASE_DUMPS=0`, and `RAYON_NUM_THREADS=1` for every child.

```sh
CARGO_HOME=/workspace/toolchains/cargo RUSTUP_HOME=/workspace/toolchains/rustup \
  /workspace/toolchains/cargo/bin/cargo build --locked --release --bin ecbench -j3
CARGO_HOME=/workspace/toolchains/cargo RUSTUP_HOME=/workspace/toolchains/rustup \
  /workspace/toolchains/cargo/bin/cargo build --locked --release --example conductor_gap_ir -j3
CARGO_HOME=/workspace/toolchains/cargo RUSTUP_HOME=/workspace/toolchains/rustup \
  /workspace/toolchains/cargo/bin/cargo test --locked --release --example conductor_gap_ir -j3
VALGRIND_LIB=/workspace/toolchains/volcano-ir/valgrind/usr/libexec/valgrind \
  target/release/examples/conductor_gap_ir run \
  research/ecbench_conductor_gap_ir_20261014 \
  research/ecbench_conductor_gap_20261012 \
  /workspace/toolchains/volcano-ir/ecbench-v2 \
  /workspace/toolchains/volcano-ir/valgrind/usr/bin/valgrind \
  > research/ecbench_conductor_gap_ir_20261014/census-v2.log 2>&1
target/release/examples/conductor_gap_ir analyze \
  research/ecbench_conductor_gap_ir_20261014 \
  research/ecbench_conductor_gap_20261012
```

The v1 commands above are historical, not instructions to pool or overwrite
either interrupted attempt. Extract each attempt's `raw-callgrind.tar.xz` in
its own directory. Preserve the raw bytes and missing receipts as recorded.

## Final library validation configuration

After the complete census, `finish-local.sh` runs native analysis and the full
`cargo test --locked --release --lib -j3 -- --test-threads=3` command. The test
build uses `CARGO_PROFILE_RELEASE_CODEGEN_UNITS=256` and a separate target
directory, `/workspace/toolchains/volcano-ir/test-target`, retaining release
optimization level 3 and all default features and tests. This bounds the
monolithic compiler work without filtering the suite. The test executable
supplies no measurement values. The heavy test build takes the repository's
`isolated_bench.py busy` lock, after the profiling process exits.
See `finish-local.status` and `library-tests-final.log` for the actual outcome.

## Archive and publication checks

The primary lossless archive is `raw-callgrind.tar.xz`, extracted with
`tar -xJf raw-callgrind.tar.xz` from the study directory. All three raw archives
passed `tar --compare`; the primary archive was extracted into a fresh scratch
study and native `analyze` regenerated both analysis files byte-identically.
The redundant interrupted ZIP containers remain locally outside the public
payload; every raw byte is published in the verified tar.xz containers.

```sh
node research/ecbench_conductor_gap_ir_20261014/render.mjs
node research/ecbench_conductor_gap_ir_20261014/render_interrupted.mjs
/workspace/toolchains/volcano-report/bin/python \
  research/ecbench_conductor_gap_ir_20261014/render_documents.py
# Run once to append the native observations to the canonical dashboard.
node research/ecbench_conductor_gap_ir_20261014/publish_views.mjs
python3 scripts/build_ic_leaderboard.py
python3 scripts/build_lab_browser.py
python3 scripts/site/panel_index.py
python3 scripts/site/test_build.py
python3 scripts/rho_status/test_rho_status.py
```

Browser generation requires all 98 committed session directories, including
the 39 initially excluded by this checkout. Only those missing session inputs
were materialized before regenerating; the resulting counts match the base.
Curve-name lint also requires the tracked Markdown/HTML and their hash-pin
files to be materialized. The rho-status suite needs network-enabled sandbox
execution for its local HTTP mock server; no remote campaign action is taken.
