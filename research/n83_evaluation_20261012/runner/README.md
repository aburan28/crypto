# Native degree-83 diagnostic runner

This Rust package orchestrates native producers sequentially. It does not perform curve arithmetic and does not establish an admitted speedup. The frozen `protocol.json` specifies the exact public-point a=1/K600 fixture and estimator before timed solves.

Build from the repository root with a local Rust toolchain. The runner requires Linux and `/usr/bin/sha256sum`:

```sh
CARGO_TARGET_DIR=target/n83-runner-reproduction cargo build --release --locked --manifest-path research/n83_evaluation_20261012/runner/Cargo.toml
CARGO_TARGET_DIR=target/n83-runner-reproduction cargo test --release --locked --manifest-path research/n83_evaluation_20261012/runner/Cargo.toml
```

The committed `protocol.json` is the exact execution record for the original run, including its absolute binary/input paths. Never edit or overwrite that frozen record to reproduce it. Before launching a reproduction, create a new protocol identity in a new directory, record the frozen protocol's SHA256 as its parent, and update the copy's `cwd`, `target_file`, and every binary/base/target path in `runs[].argv` for the local checkout. Record the chosen CPU and resource envelope before timing. Preserve the frozen seeds, arm order, caps and analysis definition when reproducing this design; describe any intentional changes as a separate experiment. Build the declared producer revisions and record their binary/input hashes before launching.

For example, choose a new absolute destination and copy the execution record there before editing the new protocol:

```sh
N83_REPRO_DIR=/absolute/path/to/new/n83-reproduction
mkdir -p "$N83_REPRO_DIR"
cp research/n83_evaluation_20261012/runner/protocol.json "$N83_REPRO_DIR/protocol.json"
sha256sum research/n83_evaluation_20261012/runner/protocol.json
```

After the new protocol has its own reproduction identity, parent digest and local paths, hash it and run this payload **inside the native isolated_bench reservation/monitor**, with no concurrent builds, tests or heavy jobs. The `panel` output directory must not already exist:

```sh
sha256sum "$N83_REPRO_DIR/protocol.json"
target/n83-runner-reproduction/release/n83-evaluation-runner panel --protocol "$N83_REPRO_DIR/protocol.json" --out-dir "$N83_REPRO_DIR/panel"
```

There are 22 children: one excluded warmup for each method, five same-seed strong/strong A/A pairs, and five alternating strong/IC pairs. Each child receives exactly one point, one CPU, a 12 GiB RLIMIT_AS limit, and a 1800-second process-group timeout. The child environment is cleared and rebuilt from a whitelist; all thread counts are one. A target scalar is never supplied through the measured child environment.

Every arm has create-new raw stdout/stderr and start/final receipts. Completed files become read-only. Receipts retain command, input and binary SHA256 digests, clean environment, exact resource envelope, monotonic whole-process wall, user/system CPU, Linux wait4 RSS (KiB multiplied by 1024), exit/signal/timeout disposition, and producer timing and semantic rows. RLIMIT_AS allocation failure is not inferred to be OOM from an exit code alone. Producer online time remains distinct from launch/load/preparation-inclusive whole-process time.

IC succeeds only with one scalar-bearing row on the requested point, producer group verification, all five finite nonnegative online stage fields, and exactly one stdout summary with rank and column count 600, target counts 1/1/0, and the declared retained-base hash. The difference between the five-stage sum and online interval is retained; a nonzero residual is explicitly marked as failing exact five-stage conformance. Independent generic arithmetic replay and independent-host admission remain separate obligations.

After the panel completes:

```sh
target/n83-runner-reproduction/release/n83-evaluation-runner analyze --pairs "$N83_REPRO_DIR/panel/pairs.jsonl" --out "$N83_REPRO_DIR/panel/statistics.json" --draws 20000 --seed 2026101283
```

The primary estimator is the geometric mean of within-pair `strong_online_ms / ic_online_ms`, computed as the exponentiated mean log ratio. The 95% interval uses 20,000 paired percentile bootstrap draws. The same resampling indexes produce a separate whole-process ratio interval. Failed, unverified, mismatched-target, mismatched-envelope, missing-timing, and inconsistent-scalar pairs remain explicit exclusions. One retained pair has no confidence interval. Five repeats concern variation of one point on one host; they do not represent independent targets or hardware classes.

`panel --resume true` skips finalized receipts only when the supplied protocol exactly matches its frozen copy. It refuses to reuse an incomplete arm directory. Preserve any incomplete attempt and use a new preregistered run identity to rerun it.

The narrow validation suite covers truncated JSON, missing rank, exit-zero failed targets, incorrect base digests, duplicate summaries, mismatched public points, and insufficient replication. Test fixtures use a task-specific directory under the system temporary directory with unique process/time filenames. The original native process timeout smoke receipt is retained under `/workspace/n83_evaluation/runner/smoke/timeout_process/`; it is a runner check, not a cryptanalytic timing row.

The repository delivery source includes formatting and a portable test-fixture path. The [original compiled scratch source](../evidence/source/runner_main_compiled.rs) is preserved for the final timing manifest. These delivery-only changes do not alter the frozen runtime binary, protocol or panel records; a reproduction records its rebuilt source and binary identities separately.
