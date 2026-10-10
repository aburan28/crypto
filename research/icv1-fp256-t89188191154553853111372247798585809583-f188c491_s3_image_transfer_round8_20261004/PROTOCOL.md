# P-256 indexed S3-image transfer protocol, round 8

Date frozen: 2026-10-04

Round 7 replaced the toy residual-depth-1 all-pairs F4 frontier with one
quadratic `S3` solve per liftable left x-coordinate.  This round tests that
exact mechanism at the actual P-256 modulus and on the selected committed
factor base rather than extrapolating from `F_7681`.

## Hypothesis

On the exact P-256 factor base `FB1h2f8621cda105`, solving
`S3(u, v, target_x) = 0` as a quadratic in `v` for every left factor-base
abscissa recovers exactly the ordered column pairs found by an independent
signed-point subtraction oracle.

The round succeeds only if, for both frozen targets:

1. the algebraic and point-oracle ordered column-pair sets are equal;
2. false negatives and false positives are zero;
3. every emitted pair independently adds to either `Q` or `-Q` with the
   signed point rows in the factor base;
4. the indexed arm performs exactly one quadratic solve per factor-base
   column and never scans the `columns^2` pair frontier; and
5. the full factor-base rebuild and every frozen FB1 identity field verify.

Any pair-set mismatch falsifies transfer.  Do not change targets or factor
base after observing the result.

## Frozen curve and factor base

- curve:
  `icv1-fp256-t89188191154553853111372247798585809583-f188c491`;
- field: the exact NIST P-256 prime;
- factor-base spec:
  `dickson-torus:depth=18,root_exponent=0x2b6fdc73dc04e7667129`;
- FB1: `FB1h2f8621cda105`;
- FB1 SHA-256:
  `2f8621cda105a51f8703dabe38dec87f3715ab855770acd7b2af0a14906f9e42`;
- columns: 131,458;
- signed points: 262,916;
- points SHA-256:
  `70ab21cbda76205ff43ab6bf85c68607be9f6b89df6399e6fa98e3489b888cf1`;
- terminal:
  `0x5b17195299a3158b93389ad04c776fff2a8bb23ca8659b5b0c2b75f9009b65e5`.

The native builder must rebuild and verify the complete dump before target
work.  Abort before measuring if any identity or row differs.

## Frozen targets

1. **planted-positive:** add the low-y (`coef = 1`) rows of columns 0 and 1.
   Abort if the result is infinity or does not verify on the curve.
2. **hash-public:** let
   `k = SHA256("icv1-fp256-t89188191154553853111372247798585809583-f188c491/s3-image-transfer-round8/target/0") mod n`,
   replacing zero with one, and set `Q = [k]G`.

The hash target is not preregistered as positive or negative; its exact pair
count is an outcome.  Target coordinates and the scalar are emitted.

## Independent algorithms

### Indexed algebraic arm

Build an `x -> column` index.  For each factor-base x-coordinate `u`, form the
round-7 quadratic coefficients over the actual P-256 field, solve the
quadratic in `v`, and retain roots present in the index.  Because the P-256
prime is 3 modulo 4, compute square roots as
`d^((p + 1) / 4)` and verify the square.  Handle linear and universal
degeneracies explicitly.

Count modular multiplications in coefficient construction, discriminant
formation, exponentiation, inversion, root construction, and square checking.
Report index lookups and emitted pairs separately.  No Gröbner solver runs.

### Signed-point subtraction reference

Index all 262,916 exact signed affine rows.  For every signed left point `P`,
compute `Q - P` with the repository's public-point variable-time P-256 group
law and look it up in the signed index.  Fold hits to ordered factor-base
column pairs.  Charge all 262,916 group subtractions and lookups per target.
This oracle shares neither the quadratic coefficients nor modular-square-root
path with the algebraic arm.

## Boundary and table

The rejected frontier is `131458^2 = 17,281,205,764` ordered x-column pairs per
target.  The indexed boundary is 131,458 quadratic solves plus emitted pairs.
The reference boundary is 262,916 signed group subtractions.

| target | columns^2 | quadratic solves | returned roots / lookups | algebra pairs | reference pairs | FN / FP | algebra muls | reference group ops | exact | class |
|:--|---:|---:|:--|---:|---:|:--|---:|---:|:--:|:--|

Pair-set SHA-256 digests and every emitted pair are committed.  Wall time is a
practicality note only and may not enter the boundary ratio.

## Scope and stop

Stop after the two targets.  This is an actual-P-256 validation of a final
known-target, two-summand `S3` gate.  It is not an `S18` solver, relation
collector, linear-algebra run, logarithm recovery, rho comparison, or ECDLP
speedup.  A successful result removes generic Gröbner and squared branch scans
from this final gate only; unknown chained intermediate coordinates remain the
next research boundary.
