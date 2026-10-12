# Native degree-83 diagnostic runner

This Rust package orchestrates native producers sequentially. It does not perform curve arithmetic and does not establish an admitted speedup. The frozen `protocol.json` specifies the exact public-point a=1/K600 fixture and estimator before timed solves.

Build with the session's cached toolchain:

```sh
PATH=/workspace/n83_evaluation/toolchain/rustup/toolchains/stable-x86_64-unknown-linux-gnu/bin:/usr/bin:/bin CARGO_HOME=/workspace/n83_evaluation/toolchain/cargo RUSTUP_HOME=/workspace/n83_evaluation/toolchain/rustup CARGO_TARGET_DIR=/workspace/n83_evaluation/runner/build cargo build --release --offline --manifest-path /workspace/n83_evaluation/runner/Cargo.toml
```

Run the panel **inside the native isolated_bench reservation/monitor**, with no concurrent builds, tests, or heavy jobs:

```sh
/workspace/n83_evaluation/runner/build/release/n83-evaluation-runner panel --protocol /workspace/n83_evaluation/runner/protocol.json --out-dir /workspace/n83_evaluation/panel
```

There are 22 children: one excluded warmup for each method, five same-seed strong/strong A/A pairs, and five alternating strong/IC pairs. Each child receives exactly one point, one CPU, a 12 GiB RLIMIT_AS limit, and a 1800-second process-group timeout. The child environment is cleared and rebuilt from a whitelist; all thread counts are one. A target scalar is never supplied through the measured child environment.

Every arm has create-new raw stdout/stderr and start/final receipts. Completed files become read-only. Receipts retain command, input and binary SHA256 digests, clean environment, exact resource envelope, monotonic whole-process wall, user/system CPU, Linux wait4 RSS (KiB multiplied by 1024), exit/signal/timeout disposition, and producer timing and semantic rows. RLIMIT_AS allocation failure is not inferred to be OOM from an exit code alone. Producer online time remains distinct from launch/load/preparation-inclusive whole-process time.

IC succeeds only with one scalar-bearing row on the requested point, producer group verification, all five finite nonnegative online stage fields, and exactly one stdout summary with rank and column count 600, target counts 1/1/0, and the declared retained-base hash. The difference between the five-stage sum and online interval is retained; a nonzero residual is explicitly marked as failing exact five-stage conformance. Independent generic arithmetic replay and independent-host admission remain separate obligations.

After the panel completes:

```sh
/workspace/n83_evaluation/runner/build/release/n83-evaluation-runner analyze --pairs /workspace/n83_evaluation/panel/pairs.jsonl --out /workspace/n83_evaluation/panel/statistics.json --draws 20000 --seed 2026101283
```

The primary estimator is the geometric mean of within-pair `strong_online_ms / ic_online_ms`, computed as the exponentiated mean log ratio. The 95% interval uses 20,000 paired percentile bootstrap draws. The same resampling indexes produce a separate whole-process ratio interval. Failed, unverified, mismatched-target, mismatched-envelope, missing-timing, and inconsistent-scalar pairs remain explicit exclusions. One retained pair has no confidence interval. Five repeats concern variation of one point on one host; they do not represent independent targets or hardware classes.

`panel --resume true` skips finalized receipts only when the supplied protocol exactly matches its frozen copy. It refuses to reuse an incomplete arm directory. Preserve any incomplete attempt and use a new preregistered run identity to rerun it.

The narrow validation suite covers truncated JSON, missing rank, exit-zero failed targets, incorrect base digests, duplicate summaries, mismatched public points, and insufficient replication. A native process timeout smoke receipt is retained under `smoke/timeout_process/`; it is a runner check, not a cryptanalytic timing row.
