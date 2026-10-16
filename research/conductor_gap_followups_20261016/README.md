# Conductor-gap follow-up experiments — 16 October 2026

This study executes the five experiments requested after the [original conductor-gap study](../conductor_gap_cm_orbits_20261015/REPORT.md). Results and completion status are tracked independently for each requirement. The existing PR is [#1888](https://github.com/aburan28/crypto/pull/1888).

| Requirement | Requested primary scope | Current status | Evidence |
|---|---|---|---|
| Interleaved crossover | 20 gaps × 4 target blocks × 3 arms = 240 attempts | Running in serial campaign | prepared/crossover.json; native/followups.rs |
| Relation-yield distribution | 20 gaps × 10 targets × 2 endpoints = 400 attempts | Queued in active serial campaign | prepared/yield.json; native/folded_events.cpp |
| Orbit geometry | 100 gaps × 2 endpoints × 4 lengths = 800 arithmetic jobs, then 40 selected/held-out solves | Order certificates verified; point/table campaign queued | results/orbit-order-certificate.jsonl |
| Conductor-size trend | 8 classes × 4 strata × 5 targets × 2 methods = 320 attempts | 32 models independently certified; timing queued | results/dose-independent-certificate.jsonl |
| Rho partition controls | 20 gaps × 5 targets × 2 endpoints × 2 salts = 400 attempts; 16 coordinate controls | Primary and coordinate panels implemented and queued | prepared/rho.json |

[Pre-run protocol](protocols/PROTOCOL.md) freezes the workload, caps, selection rules, timing boundary, clustering and censor handling. Native sources prepare public unknown-log targets, run the collector, and independently replay orbit points, relation events, rank progression and recovery. An unchanged baseline collector remains available for instrumentation checks. Python is not used for arithmetic, orchestration, verification or analysis.

The first construction script had a GP multiline syntax error; its raw error output is preserved. The first small control correctly rejected a query budget larger than its subgroup order; its raw input/error and corrected successful run are preserved. Neither is counted as a completed primary experiment. PARI subsequently verified all 32 models against every conductor-divisor class polynomial and certified all 100 scalar-orbit orders.

Before launch, run native controls, example tests, strict Clippy and required repository checks. Freeze inputs/executable hashes and commit the protocol. The campaign will run one measured process at a time, pinned to CPU 0; full process CPU includes native work and Rust replay. There is no physical-host isolation claim.

Rebuild with `cargo +1.98 build --release --example conductor_gap_followups`; use the binary's `prepare`, `controls`, `certificate-inputs`, `freeze`, `campaign` and `replay DIR` commands. Never overwrite a sealed attempt. The native scheduler checkpoints immutable receipts and retains failures. Its `results/execution.jsonl` is the live journal once launched. Campaign PID 17400 launched after the input/source freeze. See results/execution.jsonl and immutable per-attempt receipts for actual status; queued work is not completed work.

The local coordination tools are unavailable, so no cairn task/lease or shared reservation is claimed. The five backlog cards remain blocked for distributed handoff while local candidate-code execution is tracked here. Outstanding work is not represented as completed.

## Prelaunch validation

The full release library suite passed with network access enabled: 4,119 passed, 0 failed, 117 ignored. The initial restricted-sandbox run had nine local-network permission failures, all preserved; the affected module separately passed 37 tests with network access. Four native follow-up tests and strict Clippy passed. Sixty paired small instrumentation trials preserved exact base, row, query and recovery outputs; three workflow-style fresh native controls replayed successfully.

The repository-wide workflow budget checker reports six inherited violations in two untouched weak-curve workflows; their bytes match origin/main. See results/workflow-budget-base-proof.json and the raw stderr. The maintained conductor workflow retains a 25-minute budget. No baseline workflows or frozen historical pins were patched to hide these failures.

## Initial execution checkpoint

The first four crossover attempts recovered their targets and passed independent replay. A large-input logging calibration paused only the scheduler, let its current child finish, and ran the unchanged collector alone before resuming. Exact bases, relations, query counts and recovery matched; collection CPU was 117.877176301 seconds with events versus 125.112412321 without. This single calibration does not establish a logging-overhead ratio, and historical per-query timing equivalence is not assumed. See calibration/production-floor-noevents/COMPARISON.md. No primary attempt was discarded or rerun.

## Censor-aware analysis and timing deviations

The frozen benchmark source and executable remain unchanged. A separate `conductor_gap_followup_posthoc` example implements common-prefix checks when wall-censored same-crater arms stop at different query counts; its three regression tests and strict Clippy passed. Use the frozen driver for `audit-all`, then `target/release/examples/conductor_gap_followup_posthoc analyze` for final statistics. The original frozen `analyze` command assumes equal complete streams and must not be used for unequal wall censors.

Artifact timestamps flag a possible roughly 0.1-second overlap between production-calibration startup and crossover-0003’s Rust verification/archive tail. Its native collector had already exited. The raw attempt is retained and full-process CPU remains descriptive L0; results/timing-deviations.json records the evidence and required whole-gap sensitivity check. The later postprocessor build began after crossover-0012’s complete wrapper finished. Neither diagnostic pause stopped a measured child.

## Checkpoint: first 182 attempts

The first 182 crossover receipts are published: 178 verified recoveries and four wall-censored observations (crossover-0175 through crossover-0178). All four censored point/event streams passed replay; none recovered a target. Crater-CM has 58 recoveries and two censors in 60 attempts, crater-recurrence has 59 recoveries and two censors in 61 attempts, and floor-recurrence has 61 recoveries in 61 attempts. These incomplete, unequal arm counts are execution status, not a comparison or a hardness result. All remaining jobs retain their original scope and caps.

The explicitly dispatched native-controls workflow passed on source head fcba0c7a, including conditional-proofs and native-controls jobs: [run 38846171917](https://github.com/aburan28/crypto/actions/runs/38846171917). This does not establish full standard CI. GitHub now reports a merge conflict with the base branch; conflict resolution and validation remain required. The PDF/figures, all five final reports, conductor-stratum analysis, whole-gap CPU sensitivity, final independent audit, 32-model canonical registration and affected catalog/dashboard/index updates are pending, not omitted or completed.
