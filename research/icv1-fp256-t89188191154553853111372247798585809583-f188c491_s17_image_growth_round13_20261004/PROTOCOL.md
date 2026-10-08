# P-256 Dickson sixteen-summand image-growth protocol, round 13

Date frozen: 2026-10-04

Round 12 kept every explicitly split algebraic step at degree 2 through eight
summands, but the largest affine x-image grew from 2 at two leaves to 8 at
four leaves.  This round measures the next balanced layer.  It builds complete
eight-leaf images and joins two of them for a known sixteen-summand target.

## Hypothesis

On the frozen depth-5 toy fibres, a balanced sixteen-summand image tree with
explicit infinity sentinels matches exact signed point addition on every
selected parent cell.  Its total quadratic-solve plus final-lookup count is at
most one sixteenth of the exhaustive signed-tuple count for each terminal.

The local maximum polynomial degree is 2 by construction.  This is not the
degree of regularity of an unsplit `S17` ideal.  The hypothesis fails on any
intermediate image mismatch, parent false negative or false positive, or
per-terminal aggregate visit ratio above 0.0625.

## Frozen dependency and field

- round-12 result:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_s9_image_growth_round12_20261004/image-result.json`;
- required SHA-256:
  `11512f946ee9256dcb978f11280f577f2fefff3f3d761feafded1b93349bbde5`;
- registered curve:
  `icv1-fp256-t89188191154553853111372247798585809583-f188c491`;
- toy field and curve: `p = 1151`,
  `y^2 = x^3 - 3x + (b_P256 mod 1151)`;
- original Dickson depth: 5;
- terminals: 0 and 369;
- leaf residual depth: 1.

Hash-check round 12 before constructing cells.  Rebuild the complete sorted
depth-4 boundary and liftable-x inventories and abort if either differs from
the hard-coded round-9--12 inventory.

## Frozen targets and parent cells

For each terminal, take the first eight nonempty boundaries in sorted order
and repeat that sequence once to obtain sixteen leaves.  Sum the first signed
row from each leaf; if the result is infinity, enumerate signed-row tuples in
lexicographic order and take the first finite sum.  This is the planted target.

The two positive cells are that doubled sequence and its reversal.  For the
negative controls, take the first two nonempty boundaries whose liftable-x
lists both have length one.  Enumerate their 16-bit boundary patterns in
lexicographic numeric order, classify each by exact signed point addition,
and retain the first two negatives.  Abort rather than change the pool if two
negative cells are unavailable.  This fixes four cells per terminal while
keeping each negative control to one x atom.

## Exact x atomization and balanced images

For each parent cell, enumerate the Cartesian product of its sixteen exact
liftable-x lists.  Each x atom represents `2^16` signed lift tuples on this
odd-prime toy curve.

Build eight two-leaf images with univariate quadratic `S3` solves.  Compose
adjacent images into four four-leaf images, then compose those into two
eight-leaf images:

- solve `S3(u,v,w)=0` for every affine `u,v` image pair;
- if either child contains identity, include the other child's affine image;
- output identity when both children contain identity or their affine x sets
  intersect, because each signed image is closed under negation.

Cross-check every four-leaf and eight-leaf image, including its identity bit,
against exact signed group addition for that x atom.  Any mismatch falsifies
the round.

For the final known target, solve one quadratic in the right eight-leaf image
coordinate for every left eight-leaf affine x and perform indexed right-image
lookups.  Include both identity routes.  Compare the classification with the
exact full sixteen-leaf image and with the independently selected parent-cell
label.  Never scan the Cartesian product of the two eight-leaf images.

## Counts and boundary table

Count coefficient, square-root, inversion and root-construction field
multiplications; quadratic solves; roots returned; final index lookups;
two-, four-, and eight-leaf affine image widths; identity flags; x atoms; and
all represented signed tuples.

| terminal | cells | x atoms | signed tuples | quadratic solves | final lookups | max 2-leaf image | max 4-leaf image | max 8-leaf image | visits / signed tuples | FN / FP | max local degree | class |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:--|---:|:--|

The structural visit count is `quadratic solves + final lookups`.  Field
multiplications and group additions remain separate units.  Report every cell
in the raw result and terminal aggregates in the prose and dashboard.

## Scope and stop

Stop after all x atoms in the eight frozen cells.  Commit raw per-cell image
statistics, exact-reference checks, deterministic digests, failures, and an
immediate byte replay.

This is a toy sixteen-summand image-growth measurement.  Even a pass does not
establish an actual P-256 `S18` solver, the degree of regularity of an unsplit
summation-polynomial system, a P-256 relation collector, logarithm recovery,
end-to-end `S`, or speedup over Pollard rho.  If eight-leaf images approach the
field's full x-domain, record that saturation even if the visit-ratio gate
passes; it is the structural stop signal for further tree-depth extrapolation.
