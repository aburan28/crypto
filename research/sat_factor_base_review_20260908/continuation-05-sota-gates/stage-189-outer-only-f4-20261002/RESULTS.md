# Stage 189: outer-only parallel fixed-X1 F4

## Decision

`REJECTED_SCREEN`. Confirmation was not run. The repository default remains
adaptive inner Rayon sections within each F4 call plus the parallel fixed-X1
outer batch; the candidate scheduling knob was reverted and its patch is
preserved.

On the frozen public `n=59, ell=9, m=3` true-negative, the current arm took
`19.204775` wall-seconds, `192.819726` total core-seconds, and
`4,124,475,392` bytes peak RSS. The outer-only arm took `19.278470` wall,
`183.832989` core-seconds, and `3,904,471,040` bytes RSS.

| metric | outer-only / current | outcome |
|:---|---:|:---|
| wall | **1.003837** | 0.38% regression |
| total core-seconds | **0.953393** | 4.66% reduction |
| peak RSS | **0.946659** | 5.33% reduction |

The frozen screen required both wall and total core ratios below `1.00`.
Outer-only fails the wall gate, so the protocol prohibits confirmation. The CPU
and RSS reductions are retained as a scheduling tradeoff, not promoted as a
speed improvement and not used to move the threshold.

## Correctness and mechanism

Both arms used twelve outer Rayon workers, one complete fixed-X1 batch, dense
exact pair selection, the deterministic linear reducer scan, five-column
`BlockTables`, and full M4RI disabled. The candidate reported 242 internally
serial F4 calls while independent fixed-X1 calls remained in the outer pool;
the control reported zero disabled calls.

Both records authenticate the same source and algebraic factor base, visit all
512 masks, skip 270 non-rational masks, construct and complete all 242 rational
systems, find zero roots, and return exhaustive `UNSAT`. Equation, term, pair,
dense-pair, divisor, matrix, basis, extraction, degree, logical-XOR, performed-
XOR, table-memory, and full-M4RI counters agree exactly. Target-subgroup
enumeration and known log labels remain false; conflicts remain `null`.

The native phase verifier passes `8/8`. The final native replay passes `19/19`
with result SHA-256
`9b5e957686d9af59a8fad3bfaf0f0a58715876ab15019a9cd7476fbc31d05967`.

## Accounting and boundary

Stage 189 contributes a measured lower bound of 11 components,
`283.337375` wall-seconds, `1,258.115755` total core-seconds, and
`5,223,251,968` bytes peak RSS. The cumulative measured campaign lower bound is
540 components, `22,990.944527` wall-seconds, `56,916.084439` core-seconds,
and `6,310,576,128` bytes maximum RSS.

Complete campaign cost remains `null` because development compilation before
the exact native build and final writes after the charged composition control
are not outer-metered. No modelled values replace them.

The same-binary direct-MITM decomposition and full automorphism-aware rho
references remain unchanged. This one-target solver-scheduling result adds no
relation-yield, unknown-scalar, full-DLP, rho-crossover, external-reproduction,
novelty, or SOTA evidence. All-seven and Koblitz index-calculus SOTA remain
false.
