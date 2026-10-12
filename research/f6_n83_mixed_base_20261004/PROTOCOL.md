# K0 n83 mixed-base six-summand capacity screen

Registered before any new inventory run on 2026-10-04. This is a structural
screen for a possible F6-IC decomposition law, not an IC candidate, solver
timing, or relation-yield measurement. The parent n83 curve, subgroup, and
actual dimension-12 terminal base are fixed by PR #1341. This branch starts
from its head `76090c3fd12c8e50ce15796542d359c10aacd9cb`.

## Hypothesis and frozen inputs

Use `KoblitzCurve::known_n83_k0()` and standard polynomial-coordinate
subspaces of dimensions 12, 14, 15, and 16. Build each geometric base and
apply the curve's public cofactor exactly as in `K0_BASE_PROTOCOL.md`.
Run one native Rust inventory process with `RAYON_NUM_THREADS=1`. The
dimension-12 row must reproduce 4,054 usable points and 2,027 signed
columns. Keep all rows including errors. The dimension-16 base must contain
every projected dimension-12 point because the underlying subspaces nest;
check this explicitly.

For each actual projected base size `B`, compute the exact multiset count
`C(B+5,6)` for six summands. For a mixed six-summand law with four points
from that base and two from the fixed dimension-12 base, compute the exact
upper bound `C(B+3,4) C(4054+1,2)`. A uniformly random subgroup target has
coverage at most the minimum of that count divided by the exact prime
subgroup order and one. These counts do not predict yield; overlap across
roles and repeated group sums can only lower actual coverage. Preserve the
geometric point count, usable point count, signed columns, and a digest of
the projected point set for each row.

## Decision rule

Admit the mixed-base *capacity only* if the dimension-12 replay and subset
checks pass and the dimension-16 mixed count reaches at least 1% of the
subgroup order. Otherwise reject this base-size/arity combination. Even an
admitted count does not admit the solver: a separate exact, complete or
budget-labelled six-summand method must find an ordinary verified relation
within a frozen resource cap. No F4/F5/F6 or rho speedup is inferred here.

The ordinary n83 solver and one-target IC measurements remain unset in
either outcome. The source hash, command, all raw rows, failures, and exact
decision calculation will be added to `RESULT.md` after the run.
