# P-256 Dickson indexed S3-image filter round 7: result

Date run: 2026-10-04

**PASS on every frozen gate.**  The indexed quadratic image recovered exactly
all 33 exhaustive-positive branch pairs from 87,040 ordered candidates.  It
produced zero false negatives, zero false positives, and no universal
degeneracies.  Every retained component completed consistently at solving
degree 3.

The candidate performed 914 quadratic solves and charged 137,887 modular/F4
field operations against 68,666,630 row-reduction field operations in the
all-pairs F4 baseline.  The aggregate ratio is 0.002008, a 497.99-fold
reduction in this frozen solver-stage unit.  The worst individual cell ratio
is 0.02157, below the preregistered 0.05 gate.

## Structural result

For fixed left coordinate `u` and known target coordinate `t`, the direct
`S3(u, v, t)` equation is quadratic in the right coordinate `v`.  Solving that
quadratic and looking its at most two roots up in an indexed factor base visits
one left factor-base x-coordinate at a time.  It does not enumerate the square
of the branch set.

| D | frozen cells | all-pairs frontier | left boundaries | liftable x / solves | retained = exact | baseline F4 ops | candidate ops | candidate / baseline | max retained degree | class |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:--|
| 5 | 4 | 1,024 | 64 | 86 | 4 | 807,800 | 13,536 | 0.016757 | 3 | stage advance |
| 6 | 4 | 4,096 | 128 | 94 | 4 | 3,230,804 | 14,714 | 0.004554 | 3 | stage advance |
| 7 | 4 | 16,384 | 256 | 250 | 6 | 12,925,397 | 34,981 | 0.002706 | 3 | stage advance |
| 8 | 4 | 65,536 | 512 | 484 | 19 | 51,702,629 | 74,656 | 0.001444 | 3 | stage advance |

The all-pairs frontier grows as `4^(D-1)`.  The indexed frontier grows as the
number of factor-base x-coordinates, at most `2^D`, plus the emitted matches.
At the depth-18 shape used in the earlier P-256 extrapolation, this replaces
`4^17 = 17,179,869,184` ordered branch pairs with `2^17 = 131,072` left
boundaries and at most two quadratic roots per liftable x-coordinate.

That is an algebraic frontier reduction for a final, known-target,
two-summand `S3` gate.  It is not a measurement at the P-256 modulus and does
not transfer automatically to an `S18` system whose intermediate
x-coordinates are unknown.

## Complete boundary table

Candidate operations are index-build multiplications, counted coefficient /
Tonelli-Shanks / inversion multiplications, and retained F4 row-reduction
field operations.  Index lookups and exhaustive-reference point additions are
reported in the raw result but are not silently converted into that unit.

| terminal | D | target | baseline pairs | boundaries | solves | retained / exact | FN / FP | baseline F4 ops | index muls | filter muls | retained F4 ops | candidate ops | ratio | max degree | correct |
|---:|---:|:--|---:|---:|---:|:--|:--|---:|---:|---:|---:|---:|---:|---:|:--:|
| 1 | 5 | positive | 256 | 16 | 22 | 2 / 2 | 0 / 0 | 202,028 | 22 | 2,600 | 1,734 | 4,356 | 0.021561 | 3 | yes |
| 1 | 5 | negative | 256 | 16 | 22 | 0 / 0 | 0 / 0 | 201,872 | 22 | 2,584 | 0 | 2,606 | 0.012909 | 0 | yes |
| 1 | 6 | positive | 1,024 | 32 | 22 | 2 / 2 | 0 / 0 | 807,862 | 22 | 2,470 | 1,734 | 4,226 | 0.005231 | 3 | yes |
| 1 | 6 | negative | 1,024 | 32 | 22 | 0 / 0 | 0 / 0 | 807,698 | 22 | 2,728 | 0 | 2,750 | 0.003405 | 0 | yes |
| 1 | 7 | positive | 4,096 | 64 | 58 | 4 / 4 | 0 / 0 | 3,231,527 | 58 | 6,950 | 3,468 | 10,476 | 0.003242 | 3 | yes |
| 1 | 7 | negative | 4,096 | 64 | 58 | 0 / 0 | 0 / 0 | 3,231,252 | 58 | 6,692 | 0 | 6,750 | 0.002089 | 0 | yes |
| 1 | 8 | positive | 16,384 | 128 | 114 | 8 / 8 | 0 / 0 | 12,925,155 | 114 | 13,700 | 6,936 | 20,750 | 0.001605 | 3 | yes |
| 1 | 8 | negative | 16,384 | 128 | 114 | 0 / 0 | 0 / 0 | 12,925,653 | 114 | 13,930 | 0 | 14,044 | 0.001087 | 0 | yes |
| 272 | 5 | positive | 256 | 16 | 21 | 2 / 2 | 0 / 0 | 202,028 | 21 | 2,422 | 1,734 | 4,177 | 0.020675 | 3 | yes |
| 272 | 5 | negative | 256 | 16 | 21 | 0 / 0 | 0 / 0 | 201,872 | 21 | 2,376 | 0 | 2,397 | 0.011874 | 0 | yes |
| 272 | 6 | positive | 1,024 | 32 | 25 | 2 / 2 | 0 / 0 | 807,540 | 25 | 2,966 | 1,734 | 4,725 | 0.005851 | 3 | yes |
| 272 | 6 | negative | 1,024 | 32 | 25 | 0 / 0 | 0 / 0 | 807,704 | 25 | 2,988 | 0 | 3,013 | 0.003730 | 0 | yes |
| 272 | 7 | positive | 4,096 | 64 | 67 | 2 / 2 | 0 / 0 | 3,231,403 | 67 | 7,769 | 1,734 | 9,570 | 0.002962 | 3 | yes |
| 272 | 7 | negative | 4,096 | 64 | 67 | 0 / 0 | 0 / 0 | 3,231,215 | 67 | 8,118 | 0 | 8,185 | 0.002533 | 0 | yes |
| 272 | 8 | positive | 16,384 | 128 | 128 | 11 / 11 | 0 / 0 | 12,925,855 | 128 | 15,145 | 9,530 | 24,803 | 0.001919 | 3 | yes |
| 272 | 8 | negative | 16,384 | 128 | 128 | 0 / 0 | 0 / 0 | 12,925,966 | 128 | 14,931 | 0 | 15,059 | 0.001165 | 0 | yes |

All positive and negative pair-set SHA-256 digests, retained component records,
index lookups, returned roots, linear-degeneracy counts, and exact-reference
addition counts are in `filter-result.json`.  There were no linear or universal
degeneracies in the frozen corpus.

## Reproduction and evidence

```bash
cargo run --release --bin p256_dickson_s3_image_filter -- \
  --round6 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_dickson_residual_scaling_round6_20261004/degree-result.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_dickson_s3_image_filter_round7_20261004/filter-result.json
```

`filter-result.json` is 29,348 bytes with SHA-256
`bd592033a4d9cbe9f1d28f825483c589b71d5e44949b3849cd5fde7863e090be`.
An immediate replay was byte-identical.

The native unit tests independently check the derived quadratic coefficients
against direct `S3` evaluation, Tonelli-Shanks on every square in `F_7681`, and
the quadratic solver against exhaustive roots, including the linear case
`u = t`.

## Decision and limits

Accept the indexed image as the replacement for all-pairs F4 at a final known-
target `S3` gate.  Calling generic Gröbner on every residual-depth-1 branch pair
is no longer the credible baseline for that gate.

This does **not** solve the user's 17-summand relation system.  For chained S3,
most intermediate target x-coordinates are variables, so their quadratic
images cannot simply be indexed against a known constant.  No full P-256
relation, logarithm, S value, rho comparison, or attack exponent is measured.
The next gate is an actual-P-256 replay on `FB1h2f8621cda105`, followed by a
separate protocol for propagating indexed images through unknown intermediate
coordinates without restoring a squared frontier.
