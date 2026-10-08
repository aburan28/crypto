# P-256 sixteen-summand image-width transfer protocol, round 14

Date frozen: 2026-10-04

Round 13 kept every explicitly split toy-field join at degree 2 through
sixteen summands, while eight-leaf images grew to 117 x-coordinates and the
exact sixteen-leaf reference reached 298.  This round transfers the image tree
itself to the actual P-256 modulus and the committed factor base.  It measures
whether fixed factor-base atoms attain the generic `2^(k-1)` x-image ceiling,
which is the state-width cost hidden by the local degree reduction.

## Hypothesis

For each of four frozen 16-column samples from `FB1h2f8621cda105`, balanced
two-, four-, eight-, and sixteen-leaf algebraic images built only with
univariate quadratic `S3` solves match exhaustive signed P-256 group addition.
Every sample has no identity collision and attains the generic affine x-image
widths `2, 8, 128, 32768` at those four layers.

The round fails on any image or identity mismatch.  A width below the generic
ceiling is recorded as a structural collision and falsifies the no-collision
hypothesis; it is not discarded or replaced.  Every algebraic solve remains
degree 2 by construction.  That statement is not the degree of regularity of
an unsplit P-256 `S17` or `S18` ideal.

## Frozen dependency, curve, and factor base

- round-13 result:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_s17_image_growth_round13_20261004/image-result.json`;
- required SHA-256:
  `2008dcb659d3120157480b6096a4873d1f9c23ee30123f4dd3a32d450f53ba1d`;
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
  `70ab21cbda76205ff43ab6bf85c68607be9f6b89df6399e6fa98e3489b888cf1`.

Hash-check round 13 and rebuild and verify the complete factor base before
sample work.  Abort on any ICV1, FB1, point-set, terminal, row, or cardinality
mismatch.

## Frozen samples

Each sample is an ordered list of 16 distinct factor-base columns.

1. `prefix`: columns 0 through 15 in ascending order.
2. `hash-0`, `hash-1`, and `hash-2`: for sample `s` and leaf `j`, hash
   `CURVE_SLUG + "/s17-image-width-round14/sample/" + s + "/leaf/" + j +
   "/counter/" + c` with SHA-256, reduce the 256-bit big-endian digest modulo
   131,458, and take the first index not already used in that sample.  Start
   `c = 0` and increment only on a collision.

Emit every selected index and its x-coordinate.  Do not reorder, reject, or
replace a sample after seeing its image widths.

## Algebraic image tree

An image is a sorted affine x-coordinate set plus an independent identity
flag.  Starting from the 16 fixed factor-base x-coordinates:

- build eight two-leaf images by solving `S3(u,v,w)=0` in `w`;
- compose adjacent images into four four-leaf images, then two eight-leaf
  images, then one sixteen-leaf image;
- for every affine `u,v` pair, solve the same univariate quadratic in `w`;
- if either child contains identity, copy the other affine set;
- set identity when both children contain identity or their affine x sets
  intersect, since signed images are closed under negation.

Use the native P-256 field implementation.  Count coefficient,
discriminant, square-root, inversion, and root-construction multiplications;
quadratic solves; roots returned; linear and universal degeneracies; and the
width and identity bit of every node.

## Independent exact reference

For every two-, four-, eight-, and sixteen-leaf node, enumerate all signed
lifts of its fixed low-y factor-base points with the repository's public-point
P-256 group law.  Deduplicate exact affine points, project them independently
to affine x-coordinates plus identity, and require equality with the algebraic
image.  Count group additions separately.  This path shares neither the
summation-polynomial coefficients nor the modular-square-root solver.

## Boundaries and model

For `k` generic signed leaves, the exact x-image ceiling is `2^(k-1)`.  Report
one row per sample:

| sample | columns | widths at 2 / 4 / 8 / 16 | ceiling ratio at 16 | quadratic solves | roots | field muls | reference additions | exact | max local degree | class |
|:--|:--|:--|---:|---:|---:|---:|---:|:--:|---:|:--|

Also report two derived storage diagnostics from frozen cardinalities:

1. raw bytes for one collision-free 16-leaf image when each x-coordinate is
   stored in 32 bytes;
2. the unordered factor-base pair frontier `B(B+1)/2`, its generic two-root
   entry count, and raw 32-byte storage.

These are structural counts and storage lower bounds, not measured wall time,
an end-to-end relation cost, or a speedup.  Keep field multiplications and
group additions in separate units.

## Scope and stop

Stop after the four samples, their complete intermediate-node cross-checks,
the frozen storage model, and an immediate byte replay.  Commit selected
columns, compact widths and counts, final-image digests, failures, and the
factor-base receipt.

This round measures actual-P-256 fixed-atom image width.  It does not enumerate
the residual-depth-1 branch-cell frontier, solve a full P-256 `S18` system,
collect a length-17 relation, recover a logarithm, price linear algebra, or
establish `S` or a rho speedup.  A pass proves local degree-2 transfer while
quantifying the state introduced by that split; it does not make the split
cheap.
