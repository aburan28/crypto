# F6 inherited-basis specialisation breakdown

Registered before implementation or timing. The n17 prepared F6-IC profile in
[#1466](https://github.com/aburan28/crypto/pull/1466) charged 104.6–104.8 ms
of the eleven-attempt T7 target's 114.8–115.2 ms F4 build to
`ReducedBasis::specialise_shared`. The hypothesis is that row rewrite and
pivot reduction, rather than layout/bookkeeping or degree-drop completion,
dominate that interval. This is a diagnostic, not a faster algorithm.

Use the same exact n17 curve, 62 usable base points, 29 folded columns,
archived public T1/T7 targets, certified prepared logs, three summands,
degree three, 8,192 node budget, 32 target trials, one Rayon thread, and
default algorithm flags as #1466. Add opt-in nested timers and counts for
layout construction, child bookkeeping and pivot remap, displaced-row
rewrite, displaced-row reduction, completion and closure. Attribute every
return path, including early refutations, and make the reported components
sum to the `specialise_shared` timer within timing-call overhead. Keep the
default flag false. Freeze the exact source, binary, candidate manifest,
inputs, workload IDs and hashes before measurement. Run two fresh processes
per target and retain all stdout, stderr, exit status, timestamps, phase
costs, operation counts, and scalar replay receipts. The five exclusive
target-online IC costs must still sum exactly to charged online wall.

Before timing, run the focused F6 prepared-target correctness tests and the
inherited-basis exact row-space tests with the timer both disabled and
enabled. A changed scalar, replay status, attempt count, solver decision,
or word-operation count fails the correctness gate and remains a failure row.
If any phase is at least 50% of the nested specialisation interval in both
T7 repetitions, optimize that phase in a follow-on matched candidate. If no
phase dominates, retain the breakdown and investigate the largest two.
The 2× complete F6 stage and one-target IC-versus-rho objectives are not
evaluated by this diagnostic; their speedups remain unset.

The local host is unisolated. Wall-time attribution is exploratory and cannot
promote a CPU speedup. No n83 ordinary relation is inferred from n17 timing.
