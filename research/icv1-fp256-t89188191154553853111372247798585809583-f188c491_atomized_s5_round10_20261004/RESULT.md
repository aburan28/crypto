# P-256 Dickson atomized S5 round 10: result

Date run: 2026-10-04

**Correctness passes; the degree hypothesis still fails.**  Atomization and
explicit identity charts recover all 16 frozen parent-cell classifications
with zero false negatives and zero false positives across 1,664 signed tuples.
All 104 affine F4 atoms complete and agree with their affine-chart reference.
However, 44 terminal-369 negative atoms still reach solving degree 3, so the
preregistered degree-at-most-2 gate fails.

The successful part is important: round 9's three missing positive cells are
exactly recovered by the left-identity chart, and the independently enhanced
indexed image solver also classifies all 16 cells correctly.  The point at
infinity is no longer an unpriced exception.

## Boundary table

Round-9 specialised F4 is the matched “before” row.  F4 field operations,
staged field multiplications, direct identity checks, and signed group-law
tuples remain distinct units.

| terminal | arm | parent cells | atoms | signed tuples | complete / correct | max solving degree | max columns | field ops | identity checks | FN / FP | class |
|---:|:--|---:|---:|---:|:--|---:|---:|---:|---:|:--|:--|
| 0 | round-9 specialised | 8 | 8 parent systems | 128 | 8 / 6 | 2 | 32 | 7,768 | 0 | 2 / 0 | incomplete chart |
| 0 | atomized affine + identity | 8 | 8 | 128 | 8 / 8 | **2** | 32 | 7,768 | 9 | 0 / 0 | correctness advance |
| 0 | enhanced staged image | 8 | local quadratics | 128 | 8 / 8 | local 2 | — | 945 | 9 index checks | 0 / 0 | stage advance |
| 369 | round-9 specialised | 8 | 8 parent systems | 1,536 | 8 / 7 | 3 | 447 | 3,011,675 | 0 | 1 / 0 | incomplete chart |
| 369 | atomized affine + identity | 8 | 96 | 1,536 | 96 / 96 atoms; 8 / 8 cells | **3** | 67 | 135,927 | 56 | 0 / 0 | rejected: degree 3 |
| 369 | enhanced staged image | 8 | local quadratics | 1,536 | 8 / 8 | local 2 | — | 3,690 | 9 index checks | 0 / 0 | stage advance |

Atomization expands the 16 parent cells to 104 affine atoms, 6.5 atoms per
parent overall and 12 per terminal-369 parent.  It cuts terminal-369 F4 work to
0.045133 of round 9's specialised count (and 0.029961 of the curve-lift
baseline) while reducing maximum columns from 447 to 67.  Those are meaningful
stage improvements, but maximum degree remains 3.

The degree histogram is:

| terminal | degree-2 atoms | degree-3 atoms | total |
|---:|---:|---:|---:|
| 0 | 8 | 0 | 8 |
| 369 | 52 | 44 | 96 |

Every degree-3 atom is negative.  They occur when a pair selects different
stored roots from a two-root boundary.  Fixing leaves with linear equations
therefore does not make the balanced affine ideal triangular enough for F4 to
refute every non-match at degree 2.

## Exact chart result

For each atom, the signed reference distinguishes witnesses with two affine
pair intermediates from witnesses where either pair cancels to infinity.  F4
is checked only against the affine subset.  Direct constant `S3` conditions
check the identity subsets:

- left identity: `x1=x2` and `S3(x3,x4,target_x)=0`;
- right identity: `x3=x4` and `S3(x1,x2,target_x)=0`.

All 104 F4 affine classifications and all 104 identity-chart classifications
match their respective signed references.  Their union matches every parent
cell.  The enhanced staged solver independently carries identity sentinels and
finds the three formerly missed cells by its left-identity route.

## Reproduction and evidence

```bash
cargo test --bin p256_dickson_balanced_s5
cargo clippy --bin p256_dickson_balanced_s5 -- -D warnings \
  -A clippy::mismatched_bit_width_type
cargo run --release --bin p256_dickson_balanced_s5 -- \
  --atomize-round9 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_balanced_s5_round9_20261004/degree-result.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_atomized_s5_round10_20261004/degree-result.json
```

The binary hash-checks round 9 before work.  It also replays round 9
byte-identically after the new mode is compiled.

`degree-result.json` is 113,724 bytes with SHA-256
`1e14e88f624bb1dadee3ea398f33d4fe0c51b9ea0956d8a4377bbccd6810bc93`.
An immediate round-10 replay was byte-identical.

## Decision and next iteration

Accept explicit identity charts as mandatory for balanced Semaev trees and
accept atomization as a large operation/column reduction.  Reject leaf
atomization as a uniform degree-2 Gröbner result.

The degree-3 remainder is now isolated to negative affine joins between fixed
pair images.  The next bounded lever is image atomization: solve each fixed
leaf pair's univariate `S3` quadratic first, split its affine image roots, and
leave only the final target-membership quadratic.  That must charge every image
atom and retain the identity charts.  It is expected to trade another explicit
component factor for maximum local degree 2; the experiment must measure that
factor rather than claiming the degree drop alone.

No P-256 `S18` degree, relation, logarithm, end-to-end S, or rho ratio is
established here.
