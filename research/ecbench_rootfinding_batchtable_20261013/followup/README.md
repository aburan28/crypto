# Target-pinned relation collection follow-up

The [sweep protocol](PROTOCOL.md) and [holdout protocol](HOLDOUT.md) froze the
registered 20-bit prime subgroup, public target seeds, arm definitions and
three measured rounds before execution. The [sweep report](reports/sweep/REPORT.md)
and [holdout report](reports/holdout/REPORT.md) link their sealed ecbench
sessions, full replay audits, exact paired checks, graph and PDF. The
[native-walk protocol](NATIVE-WALK.md) and [native report](reports/native/REPORT.md)
cover a separate fresh-target comparison, its exact stream check, replay,
graph and PDF.
The [validation receipt](VALIDATION.md) binds the release binary, full
replays, paired checks, focused tests and preserved full-suite failures.

The holdout's mean cold counted group-addition equivalents were 7,613.33
(B40), 6,980.42 (B48) and 7,970.50 (B64). B48/B64 was 0.8758 with a
target-cluster 95% interval of 0.8191–0.9370; B48/B40 was 0.9169 with
interval 0.8225–1.0150. The B48 relation phase averaged 4,603.17 of
6,980.42 counted units. Across the 24 measured holdout executions per arm,
every accepted relation was an independent row: 605/605 for B40, 662/662
for B48 and 896/896 for B64. The observed relation hit fractions were
0.575%, 0.854% and 1.466% respectively. The useful first optimization
target is therefore the cost and hit rate of relation probes, along with
the setup tradeoff, rather than elimination of dependent rows at this size.

The [pair-table occupancy diagnostic](pair-table-occupancy.txt) found 1,600,
2,304 and 4,084 distinct keys at B40, B48 and B64, with 40, 48 and 64
batched field inversions in construction. The heuristic uniform-point hit
fractions `2 * keys / 562427` are 0.569%, 0.819% and 1.452%, close to the
observed fractions. This calculation is a diagnostic, not a distribution
proof. It suggests that changing jump scheduling alone is unlikely to
increase hit yield substantially without changing table coverage or the
relation representation. The B48 holdout averaged 3,219 walk steps and
467 pooled-offset group operations inside relation collection per measured
execution; both matter when pricing a native implementation.

These are lower-bound counted costs on one curve size. Native operations
remain incompletely priced, every local wall measurement is L0, and no
cross-host validation is claimed. See the [complete-cost problem](../../../problems/rootfinding-complete-cost-calibration.md).
