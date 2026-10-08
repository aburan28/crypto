# P-256 Dickson atomized S5 degree protocol, round 10

Date frozen: 2026-10-04

Round 9 lowered the terminal-0 balanced four-summand systems from solving
degree 3 to 2, but terminal 369 remained at 3 whenever a selected leaf
polynomial retained two stored roots.  It also proved that three exhaustive
group-law witnesses pass through the point at infinity, outside the affine
intermediate chart.  This round tests the two fixes together.

## Hypothesis

Splitting every residual-depth-1 leaf into its exact stored x-coordinate
lowers every affine atomized F4 subcomponent to solving degree at most 2 on
both frozen terminals.  Unioning those affine subcomponents with explicit
left- and right-identity charts recovers the exhaustive signed four-point
classification on all 16 round-9 cells with zero false negatives and zero
false positives.

The hypothesis fails if any F4 subcomponent is incomplete or incorrectly
classifies its affine chart, if maximum solving degree exceeds 2, or if the
chart union differs from the full group-law reference.

## Frozen dependency and inputs

- round-9 result:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_balanced_s5_round9_20261004/degree-result.json`;
- required SHA-256:
  `44f8412ff1c4dec61d0fd05bc2d66f8788e14da4063239e7431e6410da67565a`;
- registered curve:
  `icv1-fp256-t89188191154553853111372247798585809583-f188c491`;
- toy field and curve: `p = 1151`,
  `y^2 = x^3 - 3x + (b_P256 mod 1151)`;
- terminals, targets, boundary tuples, positivity labels, and exact liftable
  x lists: read unchanged from round 9;
- monomial order: graded reverse lexicographic;
- F4 degree bound: 8;
- per-atom stop: 20 seconds.

Abort before solving if the hash or any field/header/count differs.  Do not
reselect cells or targets.

## Atomized affine systems

For a round-9 cell with leaf root lists `R1,...,R4`, enumerate the complete
Cartesian product `R1 × R2 × R3 × R4` in lexical order.  For atom
`(r1,r2,r3,r4)`, keep the 15-variable round-9 specialised balanced system but
replace every leaf polynomial by the linear equation `xi-ri=0`.  Keep the two
variable-target quadratic `S3` blocks and final constant-target block
unchanged.  This is 15 variables and 16 equations of input degree 2.

For each atom, exhaust all 16 signed lifts of its four x-coordinates.  Record
three reference booleans:

1. `affine`: a sum equals `Q` or `-Q` and both pair intermediates are affine;
2. `identity`: a sum equals `Q` or `-Q` and at least one pair intermediate is
   the point at infinity;
3. `full = affine or identity`.

F4 is compared only with `affine`, the chart its variables represent.  The
cell classification is the OR of every atom's consistent affine F4 result and
its independently evaluated identity chart.

## Explicit identity charts

For an atom, the left pair can cancel exactly when `r1 = r2`; in that chart the
right pair must sum to `Q` or `-Q`, equivalently
`S3(r3,r4,target_x)=0`.  The right-identity chart is symmetric.  Both-pairs
identity cannot reach the frozen finite target.

Evaluate these constant `S3` conditions directly and compare them with the
signed group-law `identity` reference for every atom.  Count each direct `S3`
evaluation and all signed tuples; do not hide chart work inside the F4 count.

## Indexed staged control

Repeat the round-9 staged pair-image join on each parent cell, now carrying an
identity sentinel in each pair image.  The join succeeds by either:

- an affine left image root whose target quadratic hits the affine right
  image;
- left identity plus `target_x` in the affine right image; or
- right identity plus `target_x` in the affine left image.

This arm may use only univariate quadratics, direct identity checks, and index
lookups.  It must match the same full exhaustive reference but remains an
algorithmic local-degree statement, not a Gröbner regularity measurement.

## Boundary and table

Round 9's matched curve-lift F4 rows are the degree/operation reference.  The
new combinatorial boundary is the number of atomized x-products and their 16
signed tuples.  Report:

| terminal | arm | parent cells | atoms | signed tuples | complete / correct | max solving degree | max columns | field ops | identity checks | FN / FP | class |
|---:|:--|---:|---:|---:|:--|---:|---:|---:|---:|:--|:--|

Preserve round 9's degree 3 and operation totals as the “before” rows.  Report
the atom expansion relative to the 16 parent systems and to the rejected
unsplit leaf products.  Units remain separate: F4 field operations, staged
field multiplications, direct identity checks, and group additions are not
silently converted.

## Scope and stop

Stop after all atoms of the 16 frozen parent cells.  Commit every atom record,
chart classification, operation count, deterministic digest, and immediate
byte replay, including any failure.

Success establishes a degree-2 exact-chart decomposition for these toy
four-summand components.  It does not price the depth-18 P-256 branch frontier,
solve `S18`, collect a length-17 relation, recover a logarithm, or establish an
end-to-end improvement over Pollard rho.  The next boundary after success is
the measured growth of atom/image counts at 8 and 16 summands, not a claim that
the four-summand result scales for free.
