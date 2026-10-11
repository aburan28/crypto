# P-256 Dickson image-atomized S5 degree protocol, round 11

Date frozen: 2026-10-04

Round 10 made every leaf x-coordinate atomic and repaired the point-at-
infinity charts, but 44 negative affine atoms still required F4 solving degree
3.  Every such remainder is a join between the two affine pair-image sets.
This round atomizes the left pair image before the final join.

## Hypothesis

For every round-10 leaf atom, solving its two fixed-leaf `S3` quadratics first
and splitting the left affine pair image into individual x-coordinates reduces
the remaining F4 join to solving degree at most 2.  Unioned with round 10's
exact identity charts, all 16 parent-cell classifications must still match the
1,664 signed-tuple reference with zero false negatives and false positives.

The gate fails on any incomplete or incorrect image atom, any solving degree
above 2, or any changed parent classification.

## Frozen dependency

- round-10 result:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_atomized_s5_round10_20261004/degree-result.json`;
- required SHA-256:
  `1e14e88f624bb1dadee3ea398f33d4fe0c51b9ea0956d8a4377bbccd6810bc93`;
- curve, toy field, terminals, targets, parent cells, leaf atoms, signed
  references, and identity-chart outcomes: read unchanged from that file;
- monomial order: graded reverse lexicographic;
- F4 degree bound: 4;
- per-image-atom stop: 5 seconds.

Abort on any hash, schema, field, target, cell-count, atom-count, or signed-
tuple-count mismatch.  Do not reselect inputs.

## Pair-image construction

For leaf atom `(x1,x2,x3,x4)`, solve

```text
S3(x1,x2,zL)=0
S3(x3,x4,zR)=0
```

as univariate quadratics over `F_1151`, using the counted round-9/10 solver.
Deduplicate and sort each affine image.  Preserve the identity flags
`x1=x2` and `x3=x4` separately; infinity is not inserted as a field value.

Cross-check every image root against signed point addition.  Any surplus or
missing affine x-coordinate falsifies the run.

## One-variable F4 join

For every `zL` in the left affine image, make one image atom with a single
variable `zR` and two equations:

```text
S3(zL,zR,target_x) = 0
product(zR-r for r in right_affine_image) = 0.
```

Both equations have degree at most 2 because a fixed leaf pair has at most two
affine image roots.  Compare F4 consistency with direct set intersection for
that exact `zL`.  A parent affine chart is positive if any of its image atoms
is consistent.

Reuse, without alteration, the round-10 identity condition for the same leaf
atom.  The parent result is the OR of all image-atom F4 results and identity
charts.

## Boundary and table

Round 10's 104 leaf atoms, 44 degree-3 atoms, 135,927 terminal-369 F4
operations, and exact chart result are the fixed boundary.  Report:

| terminal | parent cells | leaf atoms | image atoms | signed tuples | complete / correct | max solving degree | max columns | image-build muls | F4 ops | identity checks | FN / FP | class |
|---:|---:|---:|---:|---:|:--|---:|---:|---:|---:|---:|:--|:--|

Report the image-atom expansion per leaf atom and per parent cell.  Field
multiplications, F4 row operations, identity checks, and group additions stay
separate units.

## Scope and stop

Stop after every image atom derived from the frozen 104 leaf atoms.  Commit
all per-atom records, pair-image roots and digests, operation counts, failures,
and an immediate byte-identical replay.

Success would establish an exact degree-2 decomposition of these toy
four-summand components by explicit image splitting.  It would also expose the
combinatorial price that replaces degree 3.  It does not show that image counts
stay controlled at 8, 16, or 17 summands, and it is not a P-256 `S18`
measurement, relation collector, logarithm recovery, S value, or rho speedup.
