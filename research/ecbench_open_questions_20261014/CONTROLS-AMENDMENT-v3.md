# Extended structural controls and failure replay, protocol 3

2026-10-14, before executing the supplemental replay or extended controls.
All nine measured specs, targets, seeds, caps and sealed sessions are unchanged.
The v2 receipt controls.json and first failed control log remain immutable.

The standard ecbench audit selects verified deterministic records. Supplement
it by replaying all 69 exhausted deterministic records and the one applicability
error with the exact measured binary, using its native hidden exec command.
Inputs come from the sealed specs/plan and must reproduce workload/method IDs.
Compare recovered/exhausted status, group-operation and native counters, phase
GAE, total GAE, factor-base identity, field-operation counts, automorphisms and
pricing/determinism flags. Exclude clocks and solver wall statistics.
Preserve every input/output outside the sealed session; fail on any stable
difference. This supplements, never replaces, the standard audit.

Add functional public-point recovery on the already declared degree7 E0
subgroup of order29: derive the point through ecbench's public target law,
pull it back to the checked genus3 Jacobian, invoke existing generic Jacobian
BSGS, and verify the returned scalar on both groups. Supply no planted log.
This is a correctness prerequisite, not a measured IC method or speed claim.
Generic square-root search at the full N83 order is not a scalable Jacobian IC.

Correct the general-norm request's identity field width to n. The API explicitly
returns unsupported; v2's identity with width3 did not test a correctly typed
input. Keep that receipt and record the correctly typed request in
controls-v3.json. Add exact ord_n(2) arithmetic for the declared degrees,
including131 beyond the existing u128 polynomial helper's range.

Native code is in the ecbench_open_questions controls-extended and
ecbench_open_questions_replay_failures examples. All original proof, full-cost,
N83, independence and hardware gates stay unresolved.
