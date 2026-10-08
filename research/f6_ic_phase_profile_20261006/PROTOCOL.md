# F6-IC complete-target stage attribution after two failed pair-index pilots

Registered before profiling code and measurement. PR #1457's attempt-local
pair index cut counted geometric additions by 21.77% on frozen n17 target
T7 but missed its online wall gate. The dependent PR #1458 shared the
index across attempts and cut additions by 91.49%, yet again missed the
online wall gate. This diagnostic asks what fraction of the complete
F6-IC target PDP interval the existing Macaulay stage and exact geometric
node oracle actually occupy. Counted point additions alone do not answer
that question.

Use the exact prepared n17 `icv1-f2m17-tm101-00378d4e` curve, 62 usable
base points, 29 folded columns, dimension-6 standard source subspace,
certified prepared logs, seed `20261004039`, one Rayon thread, 32 maximum
target trials, degree-three inherited F4, and the prior archived T1 and
T7 public targets. T1 is the one-attempt control; T7 uses 11 attempts.
Compare inherited F4, original F6-IC, and the PR #1457 attempt-local pair
arm, all from one new binary. These are previously seen inputs used for
mechanistic diagnosis, not an independent speedup validation. Freeze new
IC1 candidate IDs, same workload IDs, exact source, binary, and input
hashes before any measured worker process. Run two fresh-process
repetitions per arm per target in alternating order. Preserve failures,
timeouts, raw outputs, status, timestamps, memory peak when available,
and five exclusive online phases.

Reset the existing `F4Profile` immediately before the target online
interval and read it immediately after. Its build, reduction, and readback
nanoseconds are descriptive timers inside the Macaulay stage, not new
exclusive online phases. In F6-IC only, time every invocation of the exact
node oracle and count those invocations, including refutations and witness
checks. This sum is also a nested diagnostic, not an exclusive online
phase. Keep the three solver algorithms and their limits unchanged.
Instrumentation overhead is charged to the observed online wall interval.
For each run compare `target_pdp_ns`, Macaulay build/reduce/readback sum,
oracle nanoseconds, oracle calls, attempts, reductions, additions, and
scalar replay. State any residual as unattributed rather than treating
it as a measured phase. Confirm that all variants return the same public
scalar and that the five primary phase costs sum to online wall.

The decision rule is diagnostic: if the oracle sum exceeds half of F6's
target-PDP interval on both T7 repetitions, prioritize oracle internals;
otherwise, if Macaulay build/reduce/readback exceeds half, prioritize the
inherited F4 kernel. If neither exceeds half, inspect specialization,
linear elimination, and setup on a later registered profile. No runtime
speedup gate is attached to this instrumentation. The Mac host is
unisolated, and this result cannot establish a controlled 2× F6 gain,
n83 transfer, an ordinary n83 relation, or IC-versus-rho speedup.
