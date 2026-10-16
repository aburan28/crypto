# P-256 target-indexed sixteen-summand compression protocol, round 15

Date frozen: 2026-10-05

Round 14 proved that fixed sixteen-leaf images on P-256 reach the generic
32,768-x ceiling.  Materialising that final image costs 1 MiB of raw
coordinates and 16,536 univariate quadratic solves per fixed atom.  This round
tests the exact alternative already suggested by the toy image joins: retain
only the two eight-leaf images and query their known-target join directly.

## Hypothesis

For each of Round 14's four frozen P-256 samples, a target-indexed join of two
exact eight-leaf images:

1. agrees with exhaustive signed group addition on one planted-positive and
   one hash-public target;
2. stores exactly 256 affine x-coordinates, versus the 32,768-coordinate
   materialised sixteen-leaf image;
3. uses exactly 280 quadratic solves for a cold build plus one target query,
   versus 16,536 solves for the full materialised tree; and
4. retains local maximum degree 2, with zero linear or universal degeneracies.

The compression gates are therefore `stored_entries / full_entries = 1/128`
and `cold_quadratic_solves / full_quadratic_solves = 280/16536`.  Any
classification mismatch, intermediate image mismatch, identity error, or
larger count falsifies the round.  The public-hash target's positive/negative
outcome is not selected in advance and must be retained as observed.

This is an exact query representation, not a lossy fingerprint filter.  It is
also not the degree of regularity of an unsplit P-256 `S17` or `S18` ideal.

## Frozen dependency, curve, and factor base

- round-14 result:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_s17_image_width_round14_20261004/width-result.json`;
- required SHA-256:
  `6eefb6b27768023af5850cd785d75ef3729484ac85e5defef9d1d596834ebccc`;
- curve:
  `icv1-fp256-t89188191154553853111372247798585809583-f188c491`;
- field: the exact NIST P-256 prime;
- factor-base spec:
  `dickson-torus:depth=18,root_exponent=0x2b6fdc73dc04e7667129`;
- FB1: `FB1h2f8621cda105`;
- FB1 SHA-256:
  `2f8621cda105a51f8703dabe38dec87f3715ab855770acd7b2af0a14906f9e42`;
- columns / signed points: 131,458 / 262,916;
- points SHA-256:
  `70ab21cbda76205ff43ab6bf85c68607be9f6b89df6399e6fa98e3489b888cf1`.

Hash-check Round 14, regenerate its prefix and three hash-derived 16-column
samples exactly, and rebuild and verify the complete factor base.  Abort on
any identity, selected-index, coordinate, or cardinality mismatch.

## Frozen targets

For each ordered 16-column sample:

1. `planted-positive`: add all 16 stored low-y factor-base rows.  Abort if the
   result is infinity or off curve.
2. `hash-public`: let
   `k = SHA256(CURVE_SLUG + "/s17-target-join-round15/sample/" + sample_name +
   "/public") mod n`, replacing zero by one, and set `Q = [k]G`.

Emit both target coordinates and the public scalar.  Do not replace a public
target if it happens to be positive.

## Candidate: retained half-images plus target query

Build the left and right eight-leaf images independently with the same
balanced degree-2 composition as Round 14.  Verify each two-, four-, and
eight-leaf intermediate against exhaustive signed P-256 group addition.

Do not compose the two eight-leaf images.  For a finite target x-coordinate
`t`, solve `S3(u,v,t)=0` in `v` for each left-image x-coordinate `u` and test
the returned roots in the indexed right-image set.  Include both identity
routes.  Record every hit and independently verify that the full signed point
set contains the target or its negation.

The two half-images are target-independent and may be reused.  Report both:

- cold count: both image builds plus one query;
- warm count: one query after the half-images exist.

Count coefficient, discriminant, square-root, inversion and root-construction
field multiplications; quadratic solves; roots; index lookups; hits; and
identity/linear/universal cases.  Field multiplications and reference group
additions remain separate units.

## Exact reference

For every half-image, enumerate all 256 signed lifts and require exact affine
x plus identity equality with the algebraic image.  For each full sample,
enumerate the 65,536 signed sixteen-leaf tuples once, deduplicate exact points,
and classify both frozen targets by membership of `Q` or `-Q`.

Every algebraic hit must correspond to a reference-positive classification.
Every reference-positive target must have at least one algebraic hit.  Record
false negatives and false positives without replacement.

## Boundary table

| sample | target | reference / candidate | retained / full entries | cold / full solves | warm solves | roots / lookups | field muls cold | FN / FP | max local degree | class |
|:--|:--|:--|---:|---:|---:|:--|---:|:--|---:|:--|

At 32 raw bytes per x-coordinate, report retained raw bytes, materialised raw
bytes, and their ratio.  Exclude container and allocator overhead from both.

## Scope and stop

Stop after eight target checks, all intermediate and full-reference checks,
compact raw JSON, deterministic hit/image digests, and an immediate byte
replay.  Preserve the hash-public outcomes whether positive or negative.

This round compresses one fixed sixteen-leaf atom.  It does not compress the
residual-depth-1 branch-cell frontier, materialise or solve the complete
factor-base pair layer, collect a P-256 length-17 relation, recover a
logarithm, or establish end-to-end `S` or a rho speedup.  A pass moves the
remaining obstruction from local image state to shared work across branch
cells.
