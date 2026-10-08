# P-256 Dickson residual-depth scaling protocol, round 6

Date frozen: 2026-10-04

Round 5 lowered the maximum certified component solving degree from 4 to 3
when branching a depth-5 Dickson fibre until only one or two chain levels
remained.  It did not test whether that degree law survives deeper trees, and
therefore left the `4^16` P-256 branch-pair extrapolation resting on one depth.

This round measures the degree and fully charged component work across depths
5 through 8.  It is a scaling and falsification round, not a P-256 attack run.

## Hypotheses

Let `f(x) = x^2 - 2`.  For an original chain depth `D` and residual depth `r`,
fixing the last `D-r` Dickson levels creates `4^(D-r)` ordered two-summand
component pairs per target.

1. **Residual-degree law:** every certified component has maximum solving
   degree at most 3 for `r in {1, 2}`, while every frozen cell at `r = 3`
   reaches maximum degree 4.
2. **Depth invariance:** those maxima depend on `r`, not on `D`, for every
   `D in {5, 6, 7, 8}` and both deterministically selected terminal fibres.
3. **Correctness:** every positive target has at least one consistent
   component, every negative target has none, and every component agrees with
   exhaustive signed-point addition.

The round passes only if all three hold and every component completes.  A
single classification mismatch falsifies correctness; any timeout, degree-cap
stop, or unresolved pair makes the affected cell inconclusive.

## Frozen field, curve, fibres, and targets

- dependency: round 5 ending at local commit `d86f5ff93`;
- registered curve:
  `icv1-fp256-t89188191154553853111372247798585809583-f188c491`;
- toy prime: `p = 7681 = 15 * 2^9 + 1`;
- toy curve: `y^2 = x^3 - 3x + (b_P256 mod 7681)`;
- original depths: `D = 5, 6, 7, 8`;
- residual depths: `r = 1, 2, 3`;
- monomial order: graded reverse lexicographic;
- F4 degree bound: 12;
- per-component stop: 5 seconds.

Terminal selection is deterministic and solver-blind.  Scan `t = 0, 1, ...`
and retain the first two values for which `f^8(x) = t` has exactly 256 roots
in `F_p` and at least 96 of those x-coordinates lift to the toy curve.  Abort
without solver measurements if fewer than two terminals qualify.  The chosen
terminal values and fibre cardinalities are emitted in the raw result.

For each `(terminal, D)`, build the complete signed factor-base point set.  In
lexicographic affine-point order select:

- the first point in the exact two-point sum set as the positive target; and
- the first curve point outside that sum set as the negative target.

Abort the cell before solving if either target does not exist.  Target choice
may not depend on a solver result.

## Systems and exhaustive reference

Each component uses round 5's grouped direct-quadratic `S3` system: two
residual Dickson chains, curve equations with square auxiliaries, one product
auxiliary, and the quadraticized `S3` equation.  The component boundary is the
fixed value after `r` chain levels.  For every ordered boundary pair, the
reference enumerates all signed curve points in both component fibres and
tests exact affine addition before F4 runs.

The native harness must summarize every component in a canonical record and
SHA-256 the ordered record stream for each cell.  The committed JSON contains
the digest, complete degree and outcome histograms, total charged field
operations, extrema, every positive component, and every mismatch or
incomplete component.  This keeps all 114,240 component decisions committed
without a repository-sized pretty-printed transcript; deterministic replay is
the independent recovery path.

## Boundary, table, and classification

The algebraic boundary is the unbranched solving degree 4.  The enumeration
boundary is the exact ordered pair count `4^(D-r)` per target; it is not an
estimate.  The result table has one row per `(terminal, D, r, target class)`:

| terminal | D | r | target | branch pairs | correct / complete | max degree | degree / 4 | max columns | total field ops | ops / pair | stream digest | classification |
|---:|---:|---:|:--|---:|:--:|---:|---:|---:|---:|---:|:--|:--|

Degree 3 at a cost of enumerating every branch remains a solver-stage degree
advance and a branch-accounting liability, not an ECDLP speedup.  Operation
counts include all positive and negative components.  Wall time is recorded
only as a practicality note and is not a headline metric.

## Stop conditions and scope

Stop after the 48 frozen cells (two terminals, four depths, three residual
depths, and two target classes), or immediately on an input-eligibility abort.
Preserve failures and partial aggregates.  Do not tune the prime, terminals,
targets, order, degree cap, or residual depths after seeing solver output.

No full P-256 Gröbner basis is run.  No relation probability, ECDLP exponent,
rho comparison, or practical P-256 claim follows from this round.  If the
residual-degree law survives, the next separately frozen experiment may test a
mathematically sound pre-solver compatibility filter against the exact branch
pair boundary.  If it fails, do not pursue that extrapolation.
