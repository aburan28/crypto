# P-256 Dickson indexed S3-image filter protocol, round 7

Date frozen: 2026-10-04

Round 6 established that residual depth 1 is the lowest fully charged
degree-3 toy strategy, but its naive ordered branch-pair frontier is
`4^(D-1)`.  This round tests whether the final two-summand `S3` gate can be
generated from one indexed pass over the left branches instead of scanning
every left/right pair.

## Hypothesis

For fixed left x-coordinate `u` and target x-coordinate `t`, the direct
P-256 `S3(u, v, t)` equation is quadratic in the right coordinate `v`:

```text
A v^2 + B v + C = 0
A = (t - u)^2
B = -2(t^2 u + t u^2 + t a + a u + 2b)
C = t^2 u^2 - 2t(a u + 2b) + a^2 - 4bu.
```

Solving this equation for every liftable `u` in each left residual-depth-1
component and looking each returned `v` up in an `x -> right boundary` index
should recover exactly the branch pairs that exhaustive signed-point addition
marks positive.

The round succeeds only if all of the following hold in every frozen cell:

1. the indexed retained-pair set equals the exhaustive positive-pair set, so
   false negatives and false positives are both zero;
2. every retained system completes under F4, is consistent, and has solving
   degree at most 3;
3. `filter modular multiplications + retained F4 row-reduction field ops` is
   below 0.05 times the round-6 all-pairs F4 field-ops baseline;
4. no identically-zero quadratic occurs and no cell times out or stops at the
   degree cap.

A false negative falsifies the method even if the cost gate passes.  A false
positive falsifies the exact-image hypothesis but may leave a sound necessary
filter; preserve and classify that outcome rather than retuning the corpus.

## Frozen corpus

- dependency: round 6 canonical artifact SHA-256
  `71d63031111ba48ff831e79430bc87e6c7bac34626d7d0de4acfb4eaa65d40f4`;
- registered curve:
  `icv1-fp256-t89188191154553853111372247798585809583-f188c491`;
- toy field and curve: `F_7681`,
  `y^2 = x^3 - 3x + 3506`;
- terminals: 1 and 272;
- original depths: 5, 6, 7, 8;
- residual depth: 1;
- monomial order: graded reverse lexicographic;
- F4 degree bound: 12;
- per-retained-component stop: 5 seconds.

Reuse the exact round-6 targets:

| terminal | D | positive target | negative target |
|---:|---:|:--|:--|
| 1 | 5 | `(9, 500)` | `(1, 1057)` |
| 1 | 6 | `(15, 1783)` | `(1, 1057)` |
| 1 | 7 | `(2, 1636)` | `(1, 1057)` |
| 1 | 8 | `(1, 1057)` | `(10, 323)` |
| 272 | 5 | `(1, 1057)` | `(8, 1810)` |
| 272 | 6 | `(19, 2530)` | `(1, 1057)` |
| 272 | 7 | `(2, 1636)` | `(1, 1057)` |
| 272 | 8 | `(1, 1057)` | `(53, 1838)` |

The harness must reject the input if any round-6 baseline cell, target,
component count, completion count, correctness count, maximum degree, or field
operation count differs from the committed artifact.

## Indexed filter and accounting

Build the right index from every liftable factor-base x-coordinate `v` to its
unique residual-depth-1 boundary `f(v)`.  For each left boundary:

1. enumerate its at most two `f(u) = boundary` roots;
2. discard non-liftable `u`;
3. form `(A, B, C)` and solve the quadratic over `F_7681`;
4. map returned liftable factor-base roots through the right index; and
5. deduplicate retained ordered boundary pairs.

Use a native counted Tonelli-Shanks implementation.  Do not use the toy
field's precomputed square-root table for the filter.  Count every modular
multiplication performed by coefficient construction, discriminants,
inversion exponentiation, Tonelli-Shanks, and boundary reconstruction.
Add the existing F4 row-reduction field-operation count for retained pairs.
Report index construction, lookups, returned roots, and exact-reference point
additions separately; they are not silently folded into the multiplication
counter.

The baseline frontier is `4^(D-1)` ordered pairs.  The candidate frontier is
the number of left boundaries plus liftable left roots and quadratic solves;
retained output pairs are charged separately.  The result table is:

| terminal | D | target | baseline pairs | left boundaries | quadratic solves | retained / exact pairs | FN / FP | baseline F4 ops | filter muls | retained F4 ops | candidate / baseline | max degree | correct | class |
|---:|---:|:--|---:|---:|---:|:--|:--|---:|---:|---:|---:|---:|:--:|:--|

## Boundary and scope

The structural boundary is the round-6 all-pairs frontier.  At original depth
18, residual depth 1 means `2^17` boundaries but `4^17 = 2^34` ordered pairs.
This round may claim a branch-frontier reduction only for the final
two-summand `S3` gate with known target x-coordinate and only if the exact-set
gate passes.  It may not claim the same reduction for an `S18` relation
system, unknown intermediate x-coordinates, relation collection, linear
algebra, or a complete P-256 ECDLP.

Stop after all 16 cells.  Preserve all filter misses, surplus pairs, degenerate
quadratics, retained F4 failures, and raw deterministic records.  No P-256
field computation or full P-256 Gröbner basis is run here.
