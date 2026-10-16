# Conductor-gap follow-up experiments — 16 October 2026

This study executes the five experiments requested after the [original conductor-gap study](../conductor_gap_cm_orbits_20261015/REPORT.md). Results and completion status are tracked independently for each requirement. The existing PR is [#1888](https://github.com/aburan28/crypto/pull/1888).

| Requirement | Requested primary scope | Current status | Evidence |
|---|---|---|---|
| Interleaved crossover | 20 gaps × 4 target blocks × 3 arms = 240 attempts | Implemented; campaign not launched | prepared/crossover.json; native/followups.rs |
| Relation-yield distribution | 20 gaps × 10 targets × 2 endpoints = 400 attempts | Implemented; campaign not launched | prepared/yield.json; native/folded_events.cpp |
| Orbit geometry | 100 gaps × 2 endpoints × 4 lengths = 800 arithmetic jobs, then 40 selected/held-out solves | Order certificates verified; point/table campaign not launched | results/orbit-order-certificate.jsonl |
| Conductor-size trend | 8 classes × 4 strata × 5 targets × 2 methods = 320 attempts | 32 models constructed and independently certified; timing not launched | results/dose-independent-certificate.jsonl |
| Rho partition controls | 20 gaps × 5 targets × 2 endpoints × 2 salts = 400 attempts; 16 coordinate controls | Primary driver implemented; coordinate panel in progress; not launched | prepared/rho.json |

[Pre-run protocol](protocols/PROTOCOL.md) freezes the workload, caps, selection rules, timing boundary, clustering and censor handling. Native sources prepare public unknown-log targets, run the collector, and independently replay orbit points, relation events, rank progression and recovery. An unchanged baseline collector remains available for instrumentation checks. Python is not used for arithmetic, orchestration, verification or analysis.

The first construction script had a GP multiline syntax error; its raw error output is preserved. The first small control correctly rejected a query budget larger than its subgroup order; its raw input/error and corrected successful run are preserved. Neither is counted as a completed primary experiment. PARI subsequently verified all 32 models against every conductor-divisor class polynomial and certified all 100 scalar-orbit orders.

Before launch, run native controls, example tests, strict Clippy and required repository checks. Freeze inputs/executable hashes and commit the protocol. The campaign will run one measured process at a time, pinned to CPU 0; full process CPU includes native work and Rust replay. There is no physical-host isolation claim.

Rebuild with `cargo +1.98 build --release --example conductor_gap_followups`; use the binary's `prepare`, `controls`, `certificate-inputs`, `freeze`, `campaign` and `replay DIR` commands. Never overwrite a sealed attempt. The native scheduler checkpoints immutable receipts and retains failures. Its `results/execution.jsonl` is the live journal once launched. No process is currently claimed to be running by this file.

The local coordination tools are unavailable, so no cairn task/lease or shared reservation is claimed. The five backlog cards remain blocked for distributed handoff while local candidate-code execution is tracked here. Outstanding work is not represented as completed.

## Prelaunch validation

The full release library suite passed with network access enabled: 4,119 passed, 0 failed, 117 ignored. The initial restricted-sandbox run had nine local-network permission failures, all preserved; the affected module separately passed 37 tests with network access. Four native follow-up tests and strict Clippy passed. Sixty paired small instrumentation trials preserved exact base, row, query and recovery outputs; three workflow-style fresh native controls replayed successfully.

The repository-wide workflow budget checker reports six inherited violations in two untouched weak-curve workflows; their bytes match origin/main. See results/workflow-budget-base-proof.json and the raw stderr. The maintained conductor workflow retains a 25-minute budget. No baseline workflows or frozen historical pins were patched to hide these failures.
