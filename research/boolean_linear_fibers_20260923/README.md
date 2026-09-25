# Exact conditional linear fibers on generated Boolean systems

This standalone mathematical study implements the prospective contract from
PR #676. It accepts generated public Boolean fixtures only. Curve inputs,
imported targets, scalar recovery, production integration and IC/rho campaigns
are absent; their unmeasured costs remain null.

The primary is complete: the dramatic gate passes **0/18** groups, the incremental
gate **1/18**, and the matched scalar-column gate **17/18**. All universal gates
are rejected. All 168 systems complete with verified results. See
[CONCLUSION.md](CONCLUSION.md) for complete costs, retained exceptions and the
unmeasured next representation hypothesis.

## Algebra and selection

The original quadratic coefficient support defines a union interaction graph.
An independent set I has no term involving two selected variables. Fixing an
outside assignment z therefore leaves an exact linear system A(z)y=b(z) in I.
The graph is formed after parity cancellation; equation coefficients are never
combined by approximate equality. Row-basis changes alone cannot remove a
nonzero coordinate from the union support of a row space.

A deterministic meet-in-the-middle search selects a maximum independent set,
breaking cardinality ties by the smallest binary mask. A subset dynamic program
on the right half and exhaustive independent subsets on the left half certify
optimality for the bounded n<=24 domain. An edgeless graph selects every variable
directly. Input encoding, graph construction and selection are charged in every
cold solve. `selection_reference.py` supplies a separate Python implementation
for evidence replay; it is not part of the timed worker.

The remaining coordinates are split into up to four low outside bits and high
bits in reflected Gray order. The constant syndrome b(z) is quadratic in z;
each column of A(z) is affine. The worker updates these coefficients exactly,
screens a complete block of up to sixteen outside assignments, then solves each
surviving linear system in order. A returned inside witness is the smallest
binary solution; it is scattered through the selected original labels and
checked on every original equation.

## Exact screens and linear solvers

Equation values are packed into at most 32 u32 bits. For a given outside point,
if `b & ~(A_0 | ... | A_(k-1))` is nonzero, an equation has zero coefficients and
right-hand side one. The fiber is impossible.

The same test is applied to fixed redundant equation maps `T_d(v)=v XOR (v>>d)`
for d=1,2,3,4,5,7,8. A solution of A*y=b must solve T_d(A)*y=T_d(b), so every
rejection remains an exact contradiction. No rank assumption is used. Four-lane
groups stop screening only when every lane is rejected. All surviving fibers
receive a complete linear solve; passing a screen is not treated as consistency.

The row reference performs ordinary Gaussian elimination, including dependent
rows and inconsistent constants. The column implementation builds an independent
column basis with witness tags. In increasing variable order, every dependent
column uses only earlier columns, so setting the omitted variables to zero yields
the smallest binary witness. Full-word tests cover equation bit 31, changing
rank, zero columns and nonunique solutions. Small-column storage is specialized
through dimension five, with a bounded general fallback.

The scalar and NEON/SSE2 screening implementations must agree on masks, original
syndrome values, screen-round counts, models, semantic work and trace. Native
screening is specialized through dimension five; wider or partial blocks use
the scalar fallback. The `fiber_zero_simd` arm keeps only the original zero-row
screen and shares the optimized column solver, isolating the extra screen cost.
It is a separate policy and does not claim matching rejection counts.

## Discovery and frozen comparison

Three immutable discovery probes use n16/n24, seeds 17 and 937, three families,
and three repetitions. The first retained the initial zero-row design. The second
added redundant-row screens and small-column storage. On the two n24 unplanted
discovery fixtures, full linear queries fell from 132,868 to 1,576 and from
188,120 to 2,318. Complete cost still regressed on several inputs. Forced inlining
in the third probe gave no consistent improvement, so the second mathematical
implementation is selected. No holdout tuning is permitted.

The primary protocol fixes 168 systems: all 144 predecessor fixtures plus 24
unused holdouts at n12/16/20/24. All 28 predecessor methods remain. The row,
scalar-column, SIMD-column and zero-only SIMD arms bring the roster to 32.
Seven paired repetitions give 37,632 observations. Rotating the first position
by 7 times the fixture index gives nine or ten appearances at each position in
each 42-fixture size cohort. This is near balance across a cohort, not complete
balance within each individual fixture.

The reference contains all 28 predecessor methods and the row, scalar-column
and zero-only SIMD controls. The candidate must exceed the pointwise fastest
reference with a paired 95% lower bound above 2.0 in every n16/20/24, family and
regression/holdout group. All eighteen groups and every verified completion are
required. Separate >1.0 incremental and >1.05 matched-backend thresholds do not
replace the dramatic gate. Any eligible candidate requires unchanged timed-source
confirmation on unused holdouts. Earlier failures remain regressions.

## Accounting and checks

Cold totals include encoding, graph selection, coefficient-plan setup, all updates
and screens, linear solves, recovery, context destruction and independent result
checks. `selection_ns`, `fiber_setup_ns` and `fiber_ns` are exclusive subsets.
Uninstrumented phases remain null. Fixture generation, the completed retained
search reference and formatting are outside arm timers but within process receipts.
Whole-worker RSS includes all methods; candidate-specific memory and calibrated
operation costs are unmeasured.

Fiber prefixes are outside assignments screened, not fully evaluated original
assignments. Every lane in a screened block is charged even when an earlier lane
produces a model. `fiber_filter_rounds` counts lane-rounds, including lanes already
rejected when the rest of their four-lane group continues. `fiber_zero_rejected`
includes original and redundant zero-row contradictions. Full membership queries,
their rejected fibers and the represented extensions ruled out are counted
separately. Rank sums are ranks of repeated conditional matrices, not accumulated
independent relations. Caps return UNKNOWN, never UNSAT.

The 49 Rust tests retain the predecessor checks and add exhaustive graph selection
through five vertices, all small linear systems through four rows and columns,
direct checks of every transformed coefficient and screen, full equation words,
every SIMD lane, complete models and censored domains. Evidence replay uses the
Python selection oracle on every fixture. Producer checks are not external review.

Run only into a fresh directory:

```sh
python3 research/boolean_linear_fibers_20260923/run.py --out /tmp/boolean-fiber-replay
```

Sources, protocol, compiler/host data, raw observations, receipts and analysis are
frozen per run. The general generated-system driver does not invoke the repository
WDSat/full-curve suite and cannot establish a full index-calculus or rho crossover.
