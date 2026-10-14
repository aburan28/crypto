# Executed commands and accounting

Compiler/build environment for the library and worker:

```sh
RUSTUP_HOME=/workspace/n83_evaluation/toolchain/rustup
CARGO_HOME=/workspace/n83_evaluation/toolchain/cargo
CARGO_TARGET_DIR=/workspace/n83_evaluation/build
CARGO_BUILD_JOBS=4
PATH=/workspace/n83_evaluation/toolchain/rustup/toolchains/stable-x86_64-unknown-linux-gnu/bin:/usr/bin:/bin
cargo test --release --locked --lib
cargo test --release --locked --example n83_rho_speed_eval
cargo build --release --locked --example n83_rho_speed_eval
```

These variables were passed to each cargo command. Full library:3992 passed,
112 ignored,0 failed,305.50 s. Worker:7 passed. Exact worker commands/attempts
are in checks/worker_validation.md. Runner commands and source hashes are in
native_runner/README.md and checks/native_runner;22 tests passed. Native C++
generator/checker and Lean commands are in verification/README.md and audit.

The following stages ran sequentially, each in its own native reservation:

```sh
/workspace/n83_evaluation/frozen_binaries/isolated_bench reserve --cpus 3 --out /workspace/n83_speed/isolation/aa.jsonl --label N83-speed-AA -- /workspace/n83_speed/frozen_binaries/runner aa --protocol /workspace/n83_speed/protocol.json --out-dir /workspace/n83_speed/panel/aa
/workspace/n83_evaluation/frozen_binaries/isolated_bench reserve --cpus 3 --out /workspace/n83_speed/isolation/arithmetic.jsonl --label N83-speed-arithmetic -- /workspace/n83_speed/frozen_binaries/runner arithmetic --protocol /workspace/n83_speed/protocol.json --out-dir /workspace/n83_speed/panel/arithmetic
/workspace/n83_evaluation/frozen_binaries/isolated_bench reserve --cpus 3 --out /workspace/n83_speed/isolation/policy.jsonl --label N83-speed-policy -- /workspace/n83_speed/frozen_binaries/runner policy --protocol /workspace/n83_speed/protocol.json --out-dir /workspace/n83_speed/panel/policy --arithmetic-selection /workspace/n83_speed/panel/arithmetic/selection.json
/workspace/n83_evaluation/frozen_binaries/isolated_bench reserve --cpus 3 --out /workspace/n83_speed/isolation/holdout.jsonl --label N83-speed-holdout -- /workspace/n83_speed/frozen_binaries/runner holdout --protocol /workspace/n83_speed/protocol.json --out-dir /workspace/n83_speed/panel/holdout --arithmetic-selection /workspace/n83_speed/panel/arithmetic/selection.json --policy-selection /workspace/n83_speed/panel/policy/selection.json
/workspace/n83_speed/frozen_binaries/runner analyze --holdout-dir /workspace/n83_speed/panel/holdout --out /workspace/n83_speed/panel/statistics.json
```

Stage stdout/stderr were captured to /workspace/n83_speed/{STAGE}.stdout.log
and {STAGE}.stderr.log, copied unchanged under records. Native stage manifests
contain every exact solver argv, clean environment, input/config/output path,
resource cap and source/binary hash. Campaign run ID:run-20261012T1543Z.

Independent replay attempt1 wrapped the checker in
`/usr/bin/time -v -o /workspace/n83_speed/verification/final_replay_time.log`.
That utility was absent, so no checker child or receipt was created. The error
and empty stdout remain in checks/replay_attempt1.*. Attempt2 used the native
isolator's own wait4 wall/CPU/RSS accounting:

```sh
/workspace/n83_evaluation/frozen_binaries/isolated_bench run --cpus 3 --out /workspace/n83_speed/isolation/verification_attempt2.jsonl --label N83-speed-independent-replay-native -- /workspace/n83_speed/verification/check_speed /workspace/repos/crypto/experiments/koblitz-single-target-n83-20261006 /workspace/n83_optimization/fixtures /workspace/n83_speed/fixtures /workspace/n83_speed/panel 139 /workspace/n83_speed/verification/final_receipt.json
```

It completed139 independent scalar replays and all fixture/count gates, exit0,
no observed contention. The separate file auditor command/attempts are retained
in verification/README.md and provenance logs. Rendering reads native records;
it computes no research statistics. Earlier evidence and PDFs are not rewritten.
