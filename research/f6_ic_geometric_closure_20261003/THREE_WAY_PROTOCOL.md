# F4/F5/F6-IC matched one-target diagnostic: preregistration

The v3 F6-IC pilot compared only with inherited F4. This follow-on measures
the complete prepared one-target online call for inherited F4, matrix F5,
and F6-IC on **the same new public target**. Hypothesis: F6-IC's exact
geometric closure reduces complete online time versus both F4 and F5. The
engineering target is at least 2x less online time than each comparator in
every paired repetition. The prior v3 2x F4 gate failed; this experiment is
diagnostic and does not reset that result.

Freeze fixture seed `20261004001` and algorithm seed `20261004002` before
fixture generation. Generate one public hash-to-curve target without supplying
a scalar. Use the certified n17 prepared logarithm state already referenced
by the v3 pilot: 62 usable points, 29 folded columns, `m=3`, degree 3,
8192 nodes, one worker, one query per trial, at most 32 trials, a 600-second
cap per process, and the same cached symbolic preparation. F4 and F6 use
their production inherited-basis/highest-free route; F5 uses its production
matrix/lowest-free route. Preserve this solver-family distinction in the
candidate manifests. Use one freshly built binary and frozen source hashes
for every arm. Commit target, manifests, workload, inputs and hashes before
the first timed process.

Run nine separate processes, sequentially, in this fixed order: R1 F4, F5,
F6; R2 F6, F4, F5; R3 F5, F6, F4. Disable artifact cache and fix Rayon to
one thread. Keep raw stdout/stderr, start/end UTC, exit code, timeouts and
failures. Check each recovered scalar by independent general-group replay;
require the same target and scalar in all completed arms. Preserve every
target-dependent failed attempt. The five exclusive phases must sum to
`online_wall_ns`: target query, PDP, relation check, descent and recovery
check. Fixture construction, imported base logs and fixed symbolic setup are
outside this interval. Report attempt status mix, reductions, splits, F6
oracle operations, median/range of online costs, phase share and paired
ratios. Since the input target count is one, repetition quantifies local
runtime variation only; it does not estimate target-to-target variance.

The reference for this IC-variant diagnostic is inherited F4 on this exact
workload. F5 is the second comparator. No fresh rho arm is specified, so
the IC/rho ratio and any end-to-end attack speedup remain unknown. The host
is a non-isolated Mac; CPU wall ratios are exploratory even when all checks
pass. A controlled speedup requires the isolated benchmark service and a
paired same-target rho run before any IC-vs-rho claim. A timed-out or
unverified arm cannot count as a win. If F6 misses the 2x complete-call
target, do not promote it as an F4/F5 replacement; retain the result and
redirect optimization to the measured bottleneck.
