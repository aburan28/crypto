# Subgroup-sized charged rho controls, frozen before V3 runs

The V2 `(p,r,round,index)=(61,19,0,0)` record in
[`RUN_61_R0_V2.jsonl`](RUN_61_R0_V2.jsonl) completed all 18 target solves,
but the first negation-rho attempt was unverified after 806,452 steps and
386,302,417 current-thread CPU ns. Retry 1 took 11 steps and 32,292 ns.
That one failed seed dominates the V2 reference; the V2 timing is retained as
a diagnostic and is not used to claim a batch crossover.

For V3 keep both public target batches, input and answer separation, field,
H selection, pairing, finite-field BSGS table, rounds and target order fixed.
For each target, cap each rho attempt at
`max(64, 16*ceil(sqrt(r)))` steps, using the V2 seed schedule
`seed_base + retry*0x01000000`, `retry=0,...,31`. Record and charge every
attempt, including failed walks. Add a six-class folded rho control using
negation and the public order-three `j=0` automorphism. Charge its subgroup
eigenvalue search once per round, then its separate attempt stream under the
same step cap and seeds. All three controls must return the same scalar for
each target. If 32 attempts exhaust either reference, preserve that run as
censored and do not assign a complete prefix comparison beyond it.

Report cold `N=1,4,16` in both groups and `N=60` in the order-61 group,
using complete current-thread CPU and wall intervals separately from the
counted primitives. Three rounds, no discarded samples. The shared-host
measurements are L0 descriptive; no isolated wall-time or scaling conclusion
follows. The original distinct `N=64` obligation remains pending because
these groups have only 18 and 60 nonidentity targets.
