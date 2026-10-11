# P-256 Dickson eight-summand image-growth protocol, round 12

Date frozen: 2026-10-04

Round 11 obtained an exact degree-2 decomposition for every frozen
four-summand component by materialising pair images and using indexed
quadratic joins.  This round tests the next boundary: whether composing two
four-summand images keeps the explicit image/atom count controlled or merely
moves the exponential frontier.

## Hypothesis

On the frozen depth-5 toy fibres, a balanced eight-summand image tree with
explicit infinity sentinels matches exhaustive signed point addition on every
selected parent cell, and its total quadratic-solve plus final-lookup count is
at most one quarter of the exhaustive signed-tuple count.

The local maximum polynomial degree is 2 by construction and must be reported
as such, not as the degree of regularity of an unsplit `S9` ideal.  The
hypothesis fails on any image mismatch, parent false negative/positive, or
per-terminal aggregate visit ratio above 0.25.

## Frozen dependency and field

- round-11 result:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_image_atomized_s5_round11_20261004/degree-result.json`;
- required SHA-256:
  `7266c05093c815734db0a0319b946f76d003e836cce0911541be5ab4d50f2e60`;
- registered curve:
  `icv1-fp256-t89188191154553853111372247798585809583-f188c491`;
- toy field and curve: `p = 1151`,
  `y^2 = x^3 - 3x + (b_P256 mod 1151)`;
- original Dickson depth: 5;
- terminals: 0 and 369;
- leaf residual depth: 1.

Hash-check round 11 before constructing cells.  Rebuild the 16 sorted depth-4
boundaries and exact liftable x/signed rows per terminal; abort if they differ
from the round-9/10 inventory.

## Frozen target and parent cells

For each terminal, take the first eight distinct nonempty boundaries.  Sum the
first signed row from each in boundary order; if the result is infinity,
enumerate signed-row tuples over those same boundaries lexicographically and
take the first finite sum.  This is the planted finite target.

The first two parent cells are the eight boundary IDs in forward and reverse
order; both must be positive.  For negative controls, enumerate octuples over
the first two nonempty boundary IDs in lexical order and retain the first two
that exhaustive signed addition classifies negative.  Abort rather than
change the boundary pool if two negatives are unavailable.

Thus each terminal has four frozen cells: two positive and two negative.

## Exact x atomization and balanced images

For each parent boundary octuple, enumerate the Cartesian product of its eight
exact liftable x lists.  Every x atom represents 256 signed lift tuples.

An image is `(sorted affine x set, identity flag)`.  Build four two-leaf images
with a univariate quadratic `S3` solve.  Compose adjacent images into two
four-leaf images:

- solve `S3(u,v,w)=0` for every affine `u,v` image pair;
- if either input contains identity, include the other's affine image;
- output identity when both inputs contain identity or their affine x sets
  intersect, because every signed image is closed under negation.

Cross-check both four-leaf images against exhaustive signed group addition for
that x atom.  Any affine-root or identity mismatch falsifies the round.

For the final known target, solve one quadratic in the right image coordinate
for every left four-leaf affine x and perform indexed right-image lookups.
Include the two identity routes exactly as in rounds 10–11.  Compare the result
with exhaustive addition of all 256 signed lifts.

## Counts and boundary table

Count coefficient, square-root, inversion and root-construction field
multiplications; quadratic solves; roots returned; final index lookups;
two-leaf and four-leaf affine image sizes; identity flags; x atoms; and all
signed tuples.  Do not scan the Cartesian product of left and right four-leaf
images at the final join.

| terminal | cells | x atoms | signed tuples | quadratic solves | final lookups | max 2-leaf image | max 4-leaf image | visits / signed tuples | FN / FP | max local degree | class |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|:--|---:|:--|

The structural visit count is `quadratic solves + final lookups`.  Field
multiplications and group additions remain separate units.  Report positive
and negative cells individually in the raw result and terminal aggregates in
the prose/dashboard.

## Scope and stop

Stop after all x atoms in the eight frozen cells.  Commit raw per-cell image
statistics, exact-reference checks, deterministic digests, failures, and an
immediate byte replay.

This is a toy eight-summand image-growth measurement.  Even a pass does not
establish acceptable growth at 16 or 17 summands, an actual P-256 `S18`
solver, a relation collector, logarithm recovery, S value, or speedup over
Pollard rho.
