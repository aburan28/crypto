# Cycle 2 — wire median solve_wall_ms into `ecbench table`

**Change:** `ecbench table` adds a `med wall ms` column (median
`solve_wall_ns` over verified counted records). Descriptive L0 host
timing; not an admitted cost unit. Makes oracle-time wins vs F4 visible
beside S without a second tool.

Verified on cycle-1 session table output: negfold ~0.24 ms vs F4 ~225 ms
on the same panel.
