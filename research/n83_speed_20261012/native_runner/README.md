# Native further-speed campaign

This add-only runner copies the previous finite campaign's native monitor and
statistical primitives. It launches one fresh single-target process at a time
on CPU 3, with one OS thread, a 12 GiB address-space limit and a 60 second cap.
Every attempt, failure and raw producer receipt remains in its new stage.
The validation scalar sidecar is read only by the parent and never passed to
measured children.

The frozen plan has 139 launches: ten A/A controls, eighteen arithmetic-kernel
controls, 81 policy runs, and thirty holdout runs. All use the Fold field kernel.
The current control is reference scalar arithmetic, 64 lockstep lanes, DP6,
escape off. The arithmetic stage checks nine matched verified exported
outcomes, counters and charges before Barrett may be selected; it does not
export a complete walk trace. Each policy must finish and verify all nine
training cases. Whole-process arithmetic mean selects the winner, with the
stage default retained when within two percent of the minimum. Every candidate
receives all nine cases. A/A runs never enter selection or holdout statistics.

Build with the repository's pinned Rust environment:

```sh
cargo test --release --locked --manifest-path research/n83_speed_20261012/native_runner/Cargo.toml
cargo build --release --locked --manifest-path research/n83_speed_20261012/native_runner/Cargo.toml
```

Each output directory and JSON output must be new. Source, executable, input,
protocol and selection digests and the complete ordered plan are written before
the first child in each stage. The human protocol freezes both later branches
before any measurement. No resume, overwrite or adaptive pruning is supported.

```sh
n83-speed-runner aa --protocol PROTOCOL.json --out-dir PANEL/aa
n83-speed-runner arithmetic --protocol PROTOCOL.json --out-dir PANEL/arithmetic
n83-speed-runner policy --protocol PROTOCOL.json --out-dir PANEL/policy --arithmetic-selection PANEL/arithmetic/selection.json
n83-speed-runner holdout --protocol PROTOCOL.json --out-dir PANEL/holdout --arithmetic-selection PANEL/arithmetic/selection.json --policy-selection PANEL/policy/selection.json
n83-speed-runner analyze --holdout-dir PANEL/holdout --out STATISTICS.json
```

The protocol schema is `n83.speed-campaign/v1`. Required fields are
`protocol_id`, `cwd`, `rho_binary`, `fixtures_dir`, `validation_sidecar`,
`generator_receipt`, and resource fields `cpu:3`, `memory_bytes:12884901888`,
`rho_timeout_s:60`. The generator receipt must declare PASS, three training
targets and five holdout targets. The sidecar uses the prior generator's
`entries` shape, with IDs `training_1` through `training_3` and `holdout_1`
through `holdout_5`; each point file contains exactly one coordinate pair.

Optional predeclared digest checks are `human_protocol` with
`human_protocol_sha256`, `expected_rho_binary_sha256`,
`expected_generator_receipt_sha256`, `expected_validation_sidecar_sha256`,
and `expected_fixture_file_sha256` keyed by fixture ID. Each is checked before
any child. `host_context_file` attaches earlier host context to fresh CPU/host
evidence. Neither a known scalar nor the sidecar path enters child argv or its
clean environment.

Training walk seeds are 202610126211 through 202610126219. Holdout seeds are
202610126241 through 202610126255. Three consecutive seeds belong to each point.
The holdout primary is the geometric mean of five current/selected ratios of
three-seed arithmetic wall means. Its native 20,000-replicate percentile
bootstrap resamples the five paired target clusters, using seed 202610126299.
The same five targets remain the statistical units. All thirty arms must
complete and verify for a primary ratio; failed attempts cost 60 seconds only
in the separately labeled capped expected-cost description. No missing or
failed time is replaced by zero.

This campaign compares two rho configurations. It makes no new IC timing,
cross-method denominator, calibrated-operation or global-optimality claim.
The nominal monitor exit-poll delay is zero to five milliseconds plus
scheduling and poll-body overhead, applied to every arm. Independent native
scalar replay is a separate correctness gate outside measured process time.
