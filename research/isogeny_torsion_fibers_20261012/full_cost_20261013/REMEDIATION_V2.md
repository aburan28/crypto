# Frozen V2 timing remediation

Date: 2026-10-13 local. Committed before editing the V2 runner or executing
the V2 run. V1 remains immutable in `raw-v1.tsv`, `run-v1.log`,
`summary-v1.tsv`, and `independent-v1.log`.

V1 completed 640/640 measured scalar checks and passed its independent
arithmetic checker. Its duplicate rho A/A arms performed identical group
work but showed aggregate wall ratios of 420.986730, 0.770853, 2.071269,
and 5.690051 across the four sizes. The first and last ratios demonstrate
that shared-host scheduling interrupted individual arms; V1 wall ratios are
not interpretable as algorithm performance.

V2 keeps the four exact pairs, all source targets, map, factor base, row
algorithm, seeds, caps, arm order, five measured rounds, warmups, and the
V1 raw schema. It adds a monotonic **current-thread CPU nanosecond** clock
from immediately before each arm's first operation through its final scalar
check. This directly prices the whole cold arm's native work in one unit,
excluding time descheduled by other processes. V2 retains wall and phase
clocks, operation counters, process RSS, and the duplicate rho A/A arm.
The independent checker must verify all V2 scalars and equal A/A counts.
The primary V2 table uses CPU time per complete solve; group operations
remain a separate diagnostic. Bootstrap paired target intervals by resampling
the eight targets and carrying all five rounds of a sampled target together.

Stop and retain V2 if the CPU A/A ratio lies outside `[0.8,1.25]` at any
size, if any row fails, or if a CPU clock is zero or exceeds its wall clock.
Meeting the A/A threshold makes a local CPU-cost observation, not an L2
wall-clock admission or independent-host reproduction. The original
isolation and second-host requirements stay open.
