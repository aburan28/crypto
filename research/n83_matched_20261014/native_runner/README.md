# Native matched workload controller

This new Rust controller preserves the prospective protocol in `../PROTOCOL.md`
and never edits prior campaigns. It generates and freezes all32 process launches
before timing: five identical rho A/A pairs, three genuine32-target shared-setup
IC/rho cohort pairs, and eight actual fresh-process singleton pairs. There are218
planned target answers including10 A/A answers. Each cohort solves one independent
point at a time; no target-dependent rho collision state is shared.

The new `examples/n83_rho_matched_eval.rs` fixes the already selected
Fold/reference-u64/256-lane/DP6/escape-off policy. It prepares the exact legacy
curve and normal basis once per process, and all declared seed-specific G-only
jump tables before target computation. Each `solve` owns a fresh distinguished-
point table and walk state. Both one-point and32-point modes use the same worker.
IC uses the existing wide `koblitz_orbit_dlp_fast` point-only interface, rank
seed29, rank threads1 and target threads1. Preparation and output are charged
to measured whole-process costs for both methods.

Build and test before any reserved timing (one compiler/profile for workers):

```sh
CARGO_BUILD_JOBS=1 cargo test --release --locked --manifest-path research/n83_matched_20261014/native_runner/Cargo.toml
CARGO_BUILD_JOBS=1 cargo build --release --locked --manifest-path research/n83_matched_20261014/native_runner/Cargo.toml
CARGO_BUILD_JOBS=1 cargo test --release --locked --example n83_rho_matched_eval
CARGO_BUILD_JOBS=1 cargo build --release --locked --example n83_rho_matched_eval --example koblitz_orbit_dlp_fast
```

The machine schema is `n83.matched-campaign/v1` with these required fields:

```json
{
  "schema": "n83.matched-campaign/v1",
  "protocol_id": "a-new-reproduction-identity",
  "cwd": "/absolute/repository",
  "rho_binary": "/absolute/frozen/new-rho-worker",
  "expected_rho_binary_sha256": "...",
  "ic_binary": "/absolute/frozen/current-IC-worker",
  "expected_ic_binary_sha256": "...",
  "base_file": "/absolute/retained-K600.jsonl",
  "expected_base_file_sha256": "...",
  "expected_base_hash": "d4a4da26d6e7cb09367de922d9806631fe49989f70d378cae76f347f8d76a9f7",
  "fixtures_dir": "/absolute/frozen/fixtures",
  "validation_sidecar": "/absolute/frozen/fixtures/validation_sidecar.json",
  "expected_validation_sidecar_sha256": "...",
  "generator_receipt": "/absolute/frozen/fixtures/generator_receipt.json",
  "expected_generator_receipt_sha256": "...",
  "expected_fixture_file_sha256": {"target_01": "...", "target_32": "..."},
  "human_protocol": "/absolute/study/PROTOCOL.md",
  "human_protocol_sha256": "...",
  "resource": {"cpu": 3, "memory_bytes": 12884901888, "timeout_s": 2100}
}
```

All32 fixture digests are mandatory. Additional descriptive fields are retained
in the frozen protocol. Point files and the parent-only validation scalar sidecar
must agree exactly. The machine generator receipt must be `PASS` for32 targets.
Every binary, input and human-protocol hash is checked before a child. The child
environment is cleared and reconstructed with a one-thread whitelist.

The separate prepare step supports committing the matrix before execution:

```sh
n83-matched-runner prepare --protocol /absolute/protocol.json --out-dir /absolute/NEW-campaign
n83-matched-runner run --protocol /absolute/protocol.json --out-dir /absolute/campaign
```

Run the second command under the native `isolated_bench` reservation/monitor.
No compilation, tests or independent replay may overlap reserved timing. The
controller's wait4 monitor pins each child to CPU3, sets12GiB RLIMIT_AS and the
same2100-second process budget for both arms, records user/system CPU and exact
child peak RSS in bytes, and terminates a timed-out process group. Its nominal
5ms exit-observation poll plus scheduling/body delay is recorded. Batch and
singleton inputs, per-seed jump preparation and all internal capped-walk
restarts remain charged; an internal lane cap is not a process failure.

Each process directory contains `started.json`, raw `stdout.jsonl`, `stderr.txt`,
`heartbeat.jsonl`, producer `targets.jsonl`, wait4 `execution.json` and normalized
`observation.json`. Existing complete receipt sets are verified against raw
hashes, target/binary/envelope bindings and a reconstruction of the observation
before a run is skipped. An interrupted directory without a complete observation
is rejected; it is preserved and cannot be silently rerun or overwritten.
Failed complete receipts are retained and the remaining predeclared runs proceed.

Collect normalized native-analysis input only after the whole matrix is present:

```sh
n83-matched-runner collect --protocol /absolute/protocol.json --out-dir /absolute/campaign --out /absolute/NEW-statistics_input.json --independent-replay-receipt /absolute/replay.json
```

The independent native receipt must be `PASS`, bind `campaign_manifest_sha256`,
cover `checked_target_rows:218`, and state `invalid_rows:0`.
It must also bind `protocol_sha256` and pass all five exported-work A/A pairs
with `checked_aa_semantic_pairs:5` and `aa_semantic_gate:true` (no full trace).
Without this receipt,
collection is permitted only as an explicitly unverified snapshot; all normalized
verification flags are false and successful-workload statistics must be withheld.
The receipt checker independently replays every scalar and IC four-point witness.
Research audit CPU/wall/memory is metered separately from production solver costs.

Whole-process cohort and singleton observations are primary comparable costs.
Stock IC producer online excludes target parsing/validation that rho online
includes, and its five stage timers have a retained outer-clock residual. Those
fields remain raw secondary diagnostics, with no silent clock correction or
claim of exact five-phase conformance. Operation units are not calibrated across
methods. Physical-host exclusivity, common S, ecbench L2 admission and a second
host remain unset.
