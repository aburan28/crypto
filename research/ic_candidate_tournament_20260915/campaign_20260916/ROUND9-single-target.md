# Round 0009 pre-registration: a cheaper shared job for both arms

Written before any measured stage of this round was run.

## What changes, and for whom

Every job in this operation pays, before either algorithm starts, for the
curve construction (irreducible polynomial, point count, group order and its
factorisation, generator search, Frobenius eigenvalue) and the public target
(hash-to-curve lift and cofactor multiple), and at the end for the general
arithmetic final check. Through round 0008 the construction ran in the
library's general big-integer arithmetic: 3–11 million instructions and
0.5–1.8 ms per job, more than the whole index-calculus part of the winner and
more than rho's own solve on the small cells. It is identical for both arms,
which run the same executable, and it dilutes every ratio.

This round's baseline (`--source-root`) is the round-0008 incumbent source
plus a single-word path for exactly that shared work
([patch](round9-fastcurve.patch)): the generator search walks the same
abscissa sequence, lifts with the same tested exact half-trace lift
(`factor_base_points_with_x`, equal to `points_with_x` including order) and
applies the same tests in `FastCurve` arithmetic; the eigenvalue check tries
the same two candidates in the same order in that arithmetic; the group order
is trial-divided in one word when it fits; the target lift and cofactor
multiple use the same tested functions. The general path remains for fields
that do not fit and is the reference the fast path is tested against on every
cell the worker accepts and on a subfield curve. The final general-arithmetic
verification of every scalar is unchanged. Rho, the incumbent and the control
all run this executable.

Development evidence before freezing (not a claim): the new executable
reproduces, byte for byte, all 184 frozen fixtures of rounds 0007 and 0008
(curve, generator, eigenvalue, orders and targets); on one case per cell
under Callgrind, the shared construction phase fell from 0.9–2.3 to
0.5–1.4 million instructions after the earlier 3–11 million, and the
instruction ratio of the incumbent to rho was 0.51–0.61.

## Objective and gates

As rounds 0007 and 0008: `--objective rho`. The only challenger is a
configuration control (`batch_trials: 4`), so the expected outcome is
*retained*, with the retained incumbent's `winner_over_rho`,
`beats_rho_strict` and `rho_parity` recorded under the merged evaluator and
checker (every solution carries its descent relation). A promotion of the
control would be reported as such.

## Parent, seed, budget

Incumbent: `tiny2` plus certificate (round-0008 winner-config) on the
baseline above, `batch_trials: 1`. Fresh seed 2026091609; target count 1;
pilot profile; 1,356 paired jobs; one pinned CPU; 8 GiB cap; 60-second
watchdog; Valgrind 3.22.0 `Ir`; native progress recorded.

## Boundary, floor, class, honesty

Unit and boundary unchanged. Base support, `m = 3`, no direct relations,
every verification obligation, the worker's phase dumps and report fields
are unchanged. What moved is shared setup cost, equally for both arms; the
rho solve itself is the shipped implementation, untouched. Class:
engineering. The native ratio remains dominated by process creation on this
virtual machine; no arithmetic-complexity, family-wide or cryptographic-size
claim follows. Fresh fixtures from the new seed; every failure retained;
panels stay separate.
