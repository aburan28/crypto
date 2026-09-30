# Source-bound F5 and rho arms on the frozen n17a1 target

The public target is `[52411,72106]` from the immutable
[target allocation](../target-panel.json). It was derived by the registered
SHA-256 x-lift/cofactor law without supplying a logarithm. Both arms receive
that exact point and its fixture seed as provenance. The target is checked
against the n17a1 subgroup independently before the first job.

The F5 arm reuses the binary and full source/dependency/build receipt from
the independently audited [complete disclosed-point pilot](../../generic-backend-recovery-pilot/RESULTS.md).
This is the exact `765c3c5f19032bd852163805f257c56babef2040` worker,
not a new build or an unreviewed replacement. Its standard-subspace dimension
six base has 63 geometric lifts, 62 distinct usable subgroup images, and 29
sign/Frobenius matrix columns. The F5 solver uses three summands, degree three,
4096 nodes per PDP attempt, dense final relation LA, eight-trial collection
batches, and at most 256 ordinary trials. Its StdRng08 seed `2026092955`
matches the static SAT arm's ordinary query seed. Every failed attempt,
dependency, rank step and target attempt must remain in the worker report.
The worker's target-dependent interval begins after reusable relation logs
are ready and includes its scalar replay.

The rho reference uses the same source-bound worker and public target, one
Rayon thread, signed/Frobenius walk, seed `2026092958`, one requested walk,
and at most 65536 iterations per restart under the source's 64-restart
policy. It must independently recover and replay the scalar. Both processes
use the same physical local host, one worker thread, and no hard outer wall
or RSS limit. The OS child-process high-water RSS is retained via `wait4`
on macOS; it is diagnostic, not a hard memory guard. Process wall includes
launch, whereas online wall is the worker's exclusive target interval.
No run is retried or replaced after a
failure. Each arm has its own candidate/reference identity, workload and run
ID. The paired point is controlled even though workload IDs differ by query
seed and method policy.

Registration seals the worker archive digest, build receipt, exact jobs,
candidate/workload manifests, runner source and this resource envelope before
either job runs. The independent `generic_admission` checks source binding,
actual base, query law, dispatch, matrix/rank/logs, phase closure and scalar
certificate for F5; `admit_rho` checks the same public-point rho dispatch,
exclusive clocks and scalar replay. A bounded incomplete report, crashed
process or failed audit stays a row with no completed online cost. These arms
alone never produce a speedup or global-tournament promotion. The frozen
four-arm allocation requires the SAT and qualified incumbent rows, same-point
timing conditions, and an order/host-contamination review before a headline
comparison.
