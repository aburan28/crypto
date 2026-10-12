# P-256 affine-bitbox factor-base degree protocol

Date frozen: 2026-10-04

This is a bounded factor-base and algebraic-system experiment for
`icv1-fp256-t89188191154553853111372247798585809583-f188c491`
(standard name P-256).  It compares the accepted depth-18 Dickson-torus
factor base with an exact affine-bitbox factor base and measures solving
degree on small-prime analogues.  It does not run a P-256 Gröbner basis and
cannot establish a P-256 ECDLP speedup.

## Hypotheses

The candidate factor base is

```text
x = offset + sum(i=0..17, 2^i b_i),  b_i in {0,1},
y^2 = x^3 - 3x + b_P256.
```

Both signs of every liftable abscissa are materialised, then folded onto one
column exactly as in the accepted `FB1h2ea06bef7f7a` dump.  The hypotheses
are, in priority order:

1. In the paired small-prime `S3` cells, the affine-bitbox arm has a strictly
   lower median F4 `solving_degree_max` than the Dickson arm at at least two
   tested depths, with no correctness mismatch.
2. The quadratic incidence lift of the same two-point summation relation has
   a lower median `solving_degree_max` than the ordinary Dickson `S3` arm.
3. The selected P-256 interval has at least 131,239 negation-folded columns,
   the accepted Dickson count, while retaining the exact 18-bit membership
   encoding.

Hypotheses 1 and 2 are degree diagnostics on analogues.  Hypothesis 3 is a
P-256 factor-base cardinality result.  None transfers the measured toy-prime
degree to P-256.

## Frozen inputs

- repository base commit:
  `5ffc37b8d8f650a0833731805d09bfeb87e2947c`;
- registered curve and FB1 conventions: the P-256 artifacts merged by PR
  `#1327`;
- P-256 candidate family: `affine-bitbox`;
- bit width and interval width: 18 and `2^18`;
- candidate offsets: `k * 2^18` for `k = 0,1,...,7`;
- selection rule: largest number of liftable abscissae, ties resolved by the
  smaller offset;
- point-key encoding:
  `prime-affine-x-plus-one-shift-sign/be33/v1`;
- small prime: `p = 127`;
- small-prime curve: `y^2 = x^3 - 3x + (b_P256 mod 127)`;
- depths: `t = 2, 3, 4, 5`;
- ordinary summation equation: Semaev `S3(x1,x2,xR) = 0`;
- monomial order: graded reverse lexicographic;
- degree bound: 12;
- per-system wall-clock stop: 30 seconds;
- maximum targets per verdict class and depth: three.

The Dickson toy abscissae are all roots in `F_127` of
`D_(2^t)(x) = 0`, represented by the triangular recurrence
`z_(j+1) = z_j^2 - 2`.  The bitbox toy interval is selected by scanning every
non-wrapping offset `0 <= o <= p - 2^t`, minimizing first the absolute
difference from the Dickson lift count and then the offset.  This is a frozen
cardinality-matching rule, not a post-measurement choice.

## Systems and targets

For each factor point, the ordinary `S3` arm has exactly `t + 1` variables.
The Dickson arm uses `x, z1, ..., z_(t-1), y`; the bitbox arm substitutes the
linear form for `x` and uses `b0, ..., b_(t-1), y`.  Both enforce exact curve
lifting.  Their maximum input degree is four.

The incidence arms add one square variable per factor point and one slope
plus one denominator inverse.  They enforce curve membership and a generic
affine addition to a fixed target entirely with quadratic equations.  The
same non-exceptional target is used in both factor-base arms.  Targets with a
zero addition denominator are excluded by the inverse equation, and the
reference enumerator applies the identical condition.

Targets are enumerated in increasing affine point-key order.  A common
positive target is decomposable by both factor bases; a common negative target
is decomposable by neither.  Take the first three of each class for which the
generic-addition reference applies.  If a class has fewer than three targets,
run all that exist and record the shortfall.  No family-specific targets enter
the paired degree statistic.

Every system is checked two ways before its degree is accepted:

1. exhaustive enumeration of the small factor base establishes the expected
   verdict;
2. the F4 result must agree when it is certified (`inconsistent`, or a
   complete basis with no pairs above the degree bound).

A timeout or degree truncation is retained but excluded from completed-cell
medians.  A wrong certified verdict falsifies the round.

## Boundary, unit, and decision table

This is a solver-stage diagnostic, so no end-to-end `S` or rho speed ratio is
reported.  The frozen degree boundary is the matched Dickson arm in the same
cell.  The primary unit is F4 `solving_degree_max`; matrix columns to the last
productive step and field multiplications are secondary costs.

The result table has one row per `(depth, formulation, family, verdict)` and
these columns:

| depth | formulation | family | verdict | targets | correct | complete | solving degree min/median/max | max columns median | field ops median | degree / matched Dickson | classification |
|---:|:--|:--|:--|---:|:--:|:--:|:--|---:|---:|---:|:--|

Classification is `advance` only for a strict completed paired degree-ratio
reduction.  Equal degree with lower secondary cost is `engineering`; a lower
input degree without lower solving degree is `relabelling`; incomplete cells
are `inconclusive`.

## P-256 artifact and FB1 identity

The selected interval is regenerated natively with `CurveParams::p256()`.
Every stored point is checked on-curve, both signs share one column with
coefficients `1` and `n-1`, and sorted wide point keys are hashed.  The FB1
preimage is the canonical `ecbench.factor_base/v1` object with:

```text
curve = icv1-fp256-t89188191154553853111372247798585809583-f188c491
family = affine-bitbox
params.bit_width = 18
params.offset = selected hexadecimal offset
params.selection = max-lifts/k-times-width/k=0..7/v1
params.point_key_encoding = prime-affine-x-plus-one-shift-sign/be33/v1
columns = measured lift count
```

The full point dump is deterministic derived data and need not be committed.
The compact result records the FB1 preimage, its SHA-256, FB1 short name,
point-key SHA-256, signed-point and column counts, selected offset, all eight
candidate counts, regeneration command, and verification receipt.

## Success and stop conditions

Success for a lower solving-degree claim requires hypothesis 1, no certified
wrong verdict, and at least two depths with completed paired positive cells.
Hypothesis 2 is reported separately and cannot rescue a failed hypothesis 1.
The P-256 candidate is retained even if hypothesis 3 fails, but it is then not
called a cardinality improvement.  Stop after the four frozen depths and eight
P-256 intervals, or earlier on an identity/on-curve failure.  Preserve every
timeout, truncation, regression, empty target class, and failed hypothesis.

All implementation and experiment execution is native Rust.  Compact raw JSON,
the human-readable decision table, commands, and hashes are committed before
the result PR is merged.
