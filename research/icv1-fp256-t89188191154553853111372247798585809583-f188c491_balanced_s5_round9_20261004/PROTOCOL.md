# P-256 Dickson balanced S5 degree protocol, round 9

Date frozen: 2026-10-04

Round 8 removed Gröbner/F4 entirely from the final known-target
two-summand `S3` gate.  The next place regularity can reappear is one level
earlier, where two pair sums have unknown x-coordinates and must join to the
known target.  This round isolates that four-summand `S5` boundary on frozen
small-prime Dickson components.

## Hypothesis

Specialising each residual-depth-1 factor-base leaf to its exact liftable
x-polynomial lowers the maximum solving degree of the balanced quadratic
four-summand F4 system by at least one relative to the curve-lift baseline,
for each frozen terminal, without changing consistency on any selected cell.

Independently, a staged pair-image solver using only univariate quadratic
`S3` solves must agree with exhaustive signed four-point addition on every
cell.  This staged arm is an algorithmic decomposition, not the Gröbner degree
of the original ideal; report its maximum local polynomial degree as 2 and do
not relabel it as a degree-of-regularity measurement.

## Frozen field, fibres, targets, and cells

- dependency: round-8 result commit `9c6a390ea`;
- registered curve:
  `icv1-fp256-t89188191154553853111372247798585809583-f188c491`;
- toy field and curve: `p = 1151`,
  `y^2 = x^3 - 3x + (b_P256 mod 1151)`;
- original Dickson depth: 5;
- terminals: 0 and 369, inherited from rounds 4–5;
- leaf residual depth: 1, hence 16 sorted depth-4 boundaries per terminal;
- monomial order: graded reverse lexicographic;
- F4 degree bound: 8;
- per-arm, per-cell stop: 20 seconds.

For each terminal, construct every boundary's exact signed affine rows.  Let
the planted target be the sum of the first signed row from each of the first
four distinct nonempty boundaries.  If that one tuple sums to infinity,
choose the lexicographically first signed-row tuple across those same four
boundaries whose sum is finite.  Abort if the result is infinity or off curve.

Enumerate ordered quadruples of nonempty boundary IDs lexicographically.
Classify each by exhaustive signed addition, accepting a cell when some
four-point sum equals `Q` or `-Q`.  Freeze the first four positive and first
four negative cells per terminal.  Abort rather than reduce a class if either
class has fewer than four cells.

## Balanced quadratic systems

Both F4 arms use the balanced tree

```text
S3(x1, x2, zL) = 0
S3(x3, x4, zR) = 0
S3(zL, zR, target_x) = 0.
```

Every variable-target `S3(x,y,z)` is quadraticised with
`s=x+y`, `q=xy`, `h=z^2`, and `k=zs`:

```text
k^2 - 4hq - 2k(q+a) - 4bz + (q-a)^2 - 4bs = 0.
```

The final constant-target `S3` uses the round-6/7 single-product
quadraticisation.  Both arms use the same variable order: four leaf blocks,
left pair block, right pair block, then the final product.

### Curve-lift baseline

For leaf boundary `c`, retain `x,y,u` and impose

```text
x^2 - 2 - c = 0
u - x^2 = 0
y^2 - ux - ax - b = 0.
```

Together with the balanced tree this is 23 variables and 24 quadratic
equations.

### Liftability-specialised candidate

Before F4, enumerate the exact liftable children `R_c` already represented by
the factor base and replace the three curve-lift equations by

```text
product(x-r for r in R_c) = 0.
```

Every selected boundary is nonempty and residual depth 1, so this polynomial
has degree 1 or 2.  Drop `y` and `u`; keep the identical balanced tree.  The
candidate has 15 variables and 16 equations.  This is an exact component
specialisation, not a relaxation: its roots are precisely the stored
factor-base abscissae for that boundary.

## Staged pair-image arm

For each `(x1,x2)` choice solve `S3(x1,x2,zL)` as a univariate quadratic and
deduplicate its roots; do the same on the right.  For every left image root,
solve `S3(zL,zR,target_x)` in `zR` and index the returned roots in the right
image.  Count coefficient, square-root, inversion, root-construction, and
index operations.  Emit a deterministic witness when a join exists.

Compare the staged classification with exhaustive signed four-point addition
on every selected cell.  The staged path may not call F4 or scan the Cartesian
product of left and right image sets.

## Boundary and result table

The regularity boundary is the curve-lift baseline on the same selected cell.
The staged structural boundary is exhaustive signed four-point products per
cell.  Report one row per terminal and arm:

| terminal | arm | cells | variables / equations | input degree | complete / correct | max solving degree | max columns | field ops | staged image solves / lookups | exhaustive signed tuples | class |
|---:|:--|---:|:--|---:|:--|---:|---:|---:|:--|---:|:--|

The degree hypothesis passes only if every F4 cell is complete and correct and
the candidate maximum is at least one below the baseline maximum for each
terminal.  A tie or any incomplete cell falsifies it.  The staged arm passes
only with zero false positives and false negatives; report it separately even
if the F4 hypothesis fails.

## Scope and stop

Stop after the 16 frozen cells.  Preserve timeouts, degree-cap stops,
mismatches, raw cell records, operation counts, pair-image digests, and an
immediate deterministic replay.

This is a small-prime, four-summand component screen.  It is not a P-256 `S18`
run, a length-17 relation, a full factor-base search, a logarithm recovery, or
an end-to-end comparison with Pollard rho.  Even if liftability specialisation
lowers degree, applying it at P-256 depth 18 still requires pricing the branch
and multi-pair image frontiers separately.
