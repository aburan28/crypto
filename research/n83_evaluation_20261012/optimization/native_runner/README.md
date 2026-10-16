# Native finite optimization campaign

This standalone Rust package freezes and runs the declared public synthetic n83 training and holdout stages. It does not edit the earlier degree-83 evidence. The human study protocol lives at `../PROTOCOL.md`; a root-authored machine protocol provides local input/binary paths and resource conditions.

Build/test from the repository root before any timed child:

```sh
CARGO_BUILD_JOBS=2 CARGO_TARGET_DIR=target/n83-optimization-runner cargo test --release --locked --manifest-path research/n83_evaluation_20261012/optimization/native_runner/Cargo.toml
CARGO_BUILD_JOBS=2 CARGO_TARGET_DIR=target/n83-optimization-runner cargo build --release --locked --manifest-path research/n83_evaluation_20261012/optimization/native_runner/Cargo.toml
```

For the active session the executable is `/workspace/n83_optimization/runner/build/release/n83-optimization-runner`. Each stage command is a payload to run under the native `isolated_bench` reservation/monitor. Stop all builds, tests and heavy jobs during timing. Create only the common output parent in advance; each stage output directory must be new.

```sh
n83-optimization-runner kernel --protocol /absolute/path/protocol.json --out-dir /absolute/path/panel/kernel
n83-optimization-runner policy --protocol /absolute/path/protocol.json --kernel-selection /absolute/path/panel/kernel/selection.json --out-dir /absolute/path/panel/policy
n83-optimization-runner holdout --protocol /absolute/path/protocol.json --kernel-selection /absolute/path/panel/kernel/selection.json --policy-selection /absolute/path/panel/policy/selection.json --out-dir /absolute/path/panel/holdout
n83-optimization-runner ic-batch --protocol /absolute/path/protocol.json --out-dir /absolute/path/panel/ic
```

The machine schema is `n83.optimization-campaign/v1`, with `protocol_id`, `cwd`, `rho_binary`, `ic_binary`, `base_file`, `expected_base_hash`, `fixtures_dir`, `validation_sidecar`, `generator_receipt`, optional `host_context_file`, and `resource` fields `cpu:3`, `memory_bytes:12884901888`, `rho_timeout_s:60`, `ic_timeout_s:1800`. If `human_protocol_sha256` is provided, its `human_protocol` file is checked before any launch. Additional descriptive metadata is retained in the frozen protocol.

Every stage creates read-only configuration files and a manifest containing the complete ordered block, source/executable/input hashes, fresh host metadata and prerequisite selection hashes before its first child. Kernel order alternates across nine target/seed cases; policy candidate order rotates by case; holdout default/selected order alternates across fifteen cases. There are exactly 18 kernel, 108 policy and 30 held-out rho children. Training seeds are 202610128311–319, three per target; held-out seeds are 202610128341–355, three per target. No partial-block winner or fastest-seed choice exists.

The native monitor adapts the earlier process launcher in a separate source file. It pins each child to CPU3, clears and whitelists a one-thread environment, sets RLIMIT_AS12GiB, records wait4 RSS in bytes and monotonic process wall/CPU, enforces process-group termination, and preserves raw stdout/stderr and semantic receipts. Its 5 ms exit poll is common to all new arms; wall may include nominal 0–5 ms exit-observation delay plus scheduler and polling-body delays, with no strict bound under arbitrary scheduling. The scalar sidecar is opened only by the orchestrator; measured worker arguments contain public points, seeds and policy files.

Kernel/policy selection requires nine verified, uncapped completions. Wrong answers reject a candidate. Fold additionally requires matching verified outcomes, walk counters and reported charges versus Wide for all nine same-input kernel cases. No step-by-step trace is exported. The primary selection cost is the arithmetic mean capped whole-process wall; verified runs charge `min(actual_wall,60)`, failures/caps charge60 and remain explicit. Only the all-nine-complete pool can win. The applicable default wins when within2% of the minimum. No eligible candidate freezes `FAILED_NO_COMPLETE_CANDIDATE` and blocks dependent stages. Immutable `selection.json` files bind the decision to the stage's records and protocol.

The stock IC command uses rank seed29 and the retained base on the same five held-out points. Its 1800 s process budget differs from rho's 60 s cap. IC success requires rank600, expected base hash, counts5/5/0, exactly five matching verified target rows and all five phase fields per point, plus the parent's known-answer check. Preparation and process costs are retained separately; this shared setup is not five cold single-target IC costs. The extractor scans its finite index once and has no20-million-probe cap. Nonzero timer residuals remain raw.

After timing and independent native replay, analyze immutable held-out records:

```sh
n83-optimization-runner analyze --holdout-dir /absolute/path/panel/holdout --ic-dir /absolute/path/panel/ic --out /absolute/path/new_holdout_statistics.json
```

The primary held-out estimator is the geometric mean default/selected ratio of three-seed arithmetic process-wall means within each of five targets, with 20,000 target-block percentile bootstrap draws and seed202610128399. Online default/selected is secondary. Prepared IC/selected-rho online uses one IC interval and the mean of three selected-rho intervals per point, with the same target-block method and seed202610128400. These directions, arithmetic expected-cost descriptions, caps and missing observations are explicit. Full-workload performance ratios require all relevant verified uncensored draws; partial or failed batches cannot supply an unqualified comparison. No zero imputation, new host population claim or calibrated common operation unit is introduced.

All outputs use create-new files. Interruptions preserve a partial block and never select it. The earlier frozen protocol/source/binaries remain unchanged. Reproduction uses a new machine protocol identity and new output directories, with local paths and source/binary hashes frozen before launching.
