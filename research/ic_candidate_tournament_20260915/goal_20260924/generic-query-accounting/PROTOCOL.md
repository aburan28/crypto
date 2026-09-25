# Generic query and terminal-descent accounting

This admission correction follows PRs 789 and 792. It does not change the
registered performance panel or consume an improvement round. Round one retained
the incumbent; two bounded improvement rounds remain. All controls use synthetic
degree-9 toy curves, with no external targets.

The generic collector currently discards decomposition statistics and retains
only witnesses. The individual-logarithm API additionally discards its entire
report when the trial budget expires. An observed execution must instead retain
one record per attempted query, including identity queries, unsupported encodings,
incomplete searches and misses. A witness is a solver result, not an independent
verification receipt. A partial pair-table miss is unresolved, never UNSAT.
Native algebra counters must remain attached to the actual attempted query.

Keep existing collection and successful-descent APIs compatible. Add opt-in query
records and an API returning the terminal descent report on both success and
failure. The tournament worker must request these records and retain failed
matrix/descent reports. Do not rerun queries to manufacture observations.

Freeze controls before execution: K_0 and K_1 over F_(2^9), factor-base index zero (and a 36-point-requested subgroup-orbit base
with base seed 43 for window controls),
query seed 2026092554, release mode, checked-in autolab dependency lock, one test
thread. Compare observed and ordinary paths on identical query ranges, including
a mid-run partition of a windowed pair collector. Replay every witness in group
arithmetic; cross-check completed negative answers with exhaustive enumeration.
Exercise zero-node budgets, unsupported layouts, failed terminal descents, and
successful descents. Require exact attempted-query closure, retained counters,
stable query order and matching certificates; stop on a failed control.

The worker's existing schema remains a legacy integration report. Query records
have their own version. Public-point input, exclusive phase accounting and an
independent generic stage auditor remain separate admission requirements. Online
and cold performance, Ir, S, rho ratios and speedup remain unknown in this change.
Do not use correctness-test elapsed time as a benchmark sample.
