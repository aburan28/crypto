# P-256 Dickson quadratic-S3 degree protocol, round 3

Date frozen: 2026-10-04

Rounds 1 and 2 established a stable small-prime solving degree of 4 for the
quartic two-summand `S3` system across the accepted terminal-zero Dickson base,
an affine-bitbox replacement, and 47 nonzero full Dickson fibres.  This round
keeps the best Dickson factor-base structure and quadraticizes `S3` directly,
rather than replacing it with the larger generic affine-addition incidence
system.

## Hypothesis

Writing `q = x1*x2` and `u_i = x_i^2` gives an equivalent quadratic system:

```text
u_i - x_i^2 = 0
y_i^2 - u_i*x_i - a*x_i - b = 0
q - x1*x2 = 0
(x1-x2)^2*xR^2
  - 2*((x1+x2)*(q+a) + 2b)*xR
  + (q-a)^2 - 4b*(x1+x2) = 0.
```

Together with the quadratic Dickson chains, every input equation has degree
at most 2.  The hypothesis is that at least one frozen variable layout has
completed F4 `solving_degree_max <= 3` at both depths 4 and 5, on both the
terminal-zero reference and the round-2 winning full fibre, with all answers
correct.

## Frozen inputs

- dependency: round 2 ending at local commit `c076118be`;
- curve identity:
  `icv1-fp256-t89188191154553853111372247798585809583-f188c491`;
- toy prime and curve: `p = 1151`,
  `y^2 = x^3 - 3x + (b_P256 mod 1151)`;
- depths and terminals:
  - depth 4: terminal zero and round-2 winner `c = 782`;
  - depth 5: terminal zero and round-2 winner `c = 369`;
- targets: the first three common-positive affine points from the corresponding
  round-2 cell, recomputed by exhaustive signed-point addition;
- monomial order: graded reverse lexicographic;
- F4 degree bound: 12;
- per-system stop: 30 seconds.

The ordinary quartic `S3` arm is rerun on every target as the paired reference.
Repeated summands and doubling are allowed.  Curve lifting and the full
Dickson fibre are exact under the same finite-field semantics as rounds 1–2.

## Frozen variable layouts

Each quadratic arm contains the same equations and differs only by the order
assigned to variables before grevlex:

1. `grouped`: `chain1, y1, u1, chain2, y2, u2, q`;
2. `aux-first`: `q, u1, u2, y1, y2, chain1, chain2`;
3. `level-interleaved`: alternating chain levels for points 1 and 2, followed
   by `y1, y2, u1, u2, q`;
4. `reverse-blocks`: `q, u2, y2, reverse(chain2), u1, y1, reverse(chain1)`.

Variable permutation changes the leading-term path but not the solution set.
The winner minimizes median completed solving degree, then median columns to
the last productive step, median field multiplications, and layout name.

## Correctness, boundary, and decision table

Exhaustive signed-point addition establishes every positive target.  A
completed F4 run must remain consistent; a certified inconsistency falsifies
the round.  Timeouts and degree truncations are preserved and excluded from
winning medians.

The boundary is the paired ordinary Dickson `S3` solving degree, expected to
remain 4.  The result table uses one row per `(depth, terminal, layout)`:

| depth | terminal | arm/layout | targets | correct | complete | input degree | solving degree min/median/max | columns median | field ops median | degree / quartic S3 | classification |
|---:|---:|:--|---:|:--:|:--:|---:|:--|---:|---:|---:|:--|

A strict ratio below 1 at both depths is a degree advance.  Equal degree with
lower secondary cost is engineering.  Merely replacing quartics with more
quadratic variables while retaining solving degree 4 is relabelling.

## Scope and stop conditions

Stop after the four layouts, two depths, two terminals, and three targets per
cell.  Preserve all regressions and incomplete runs.  Commit compact raw JSON,
the decision table, commands, and hashes.

This is a native Rust small-prime solver-stage experiment.  It neither
measures the P-256 Gröbner degree nor claims an ECDLP speedup or exponent
improvement.  The P-256 factor-base inventory remains round 2's verified
`FB1h2f8621cda105` unless a later, separately preregistered search replaces it.
