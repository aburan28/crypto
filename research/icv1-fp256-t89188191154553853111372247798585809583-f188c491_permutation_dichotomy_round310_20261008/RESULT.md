# P-256 same-base permutation image dichotomy, round 310: result

## Verdict

No unrestricted same-factor-base column permutation can provide both a
globally one-dimensional joint relation image and two independent useful
rows.  On the complete fixed-sum coefficient difference space, rank one is
**exactly** affine log transport.  That case either cannot hit an unaligned
product target or exposes the target logarithm directly; every non-affine
permutation has a two-dimensional image.

This closes the global linear-permutation escape left open by Round 309.  It
does not close nonlinear selectors restricted to a proper coefficient subset.
Such a selector remains admissible only if it accounts exhaustively for the
discarded domain and retains two-row information rank.

The complete `F_5` census checked 3,744 profile/permutation pairs and 468,000
affine-slice rows.  Among the 3,720 nonconstant cases, all 600 affine
transports had rank one and all 3,120 non-affine transports had rank two.
The 24 constant controls had rank zero.  Bucket, span, row-order, affine-rank,
false-positive, and false-negative discrepancies were all zero.

The P-256 labelled control independently obtains rank 16 for S17 sign-row
differences, rank two for the half-turn image, rank one for the identity
image, and rank two for the planted augmented rows over both `2^61-1` and the
P-256 subgroup order.  Both group equations replay exactly.  This is a
control, not a relation on `FB1h6255ce9746fe`.

No factor base is promoted.  Round 309's exact S17 product-target mean remains
`2^-271.270334`; its arity-40 balanced baseline remains `2^151.122587` times
rho, and structured degree remains unset.  No unplanted full-depth P-256
relation was attempted.

## Requirement-to-evidence status

| requested direction / gate | status | Round-310 evidence or gap |
|:--|:--:|:--|
| global same-base selector with all columns variable | closed for linear permutations | rank-one image iff transported log vector is affine on the full fixed-sum span |
| non-affine image near one group dimension | failed | every non-affine toy case has exact rank two; theorem identifies the general obstruction |
| two useful rows on one original log vector | verified as planted control | augmented rank two, no auxiliary log block, exact P-256 replay |
| actual non-labelled relation mechanism | failed | zero actual relations; P-256 scalar labels are explicitly controls only |
| zero false positives / negatives | verified | complete toy census reports zero of each and zero discrepancies |
| structured degree `<=5` | not established | Round-309 derived-base degree remains unset; no nonlinear restricted family exists to measure |
| complete cost below `2^120` and rho | failed | imported exact arity-40 balanced baseline is `2^151.122587` times rho |
| cost per usable relation below `2^103` | failed | no qualifying relation oracle |
| peak materialized storage below `2^50` | failed for imported exact baseline | Round 309 projects `2^285.358650` bytes |
| complete collection, sparse algebra, and recovery | unset | no qualifying relation source |
| no discarded branch counted exhaustive | verified | complete toy slice is enumerated; nonlinear restricted subsets receive no projection credit |
| unplanted P-256 relation | correctly not attempted | promotion gates fail |

## Linear-image certificate

Let `ell` be the factor-base log vector, `T` a column permutation, and

```text
A_s = {c : 1 dot c = s},
H   = A_s - A_s = {h : 1 dot h = 0}.
```

The difference image of the paired relation coordinates is

```text
h |-> (ell dot h, T(ell) dot h),   h in H.
```

Its rank is below two exactly when the two restricted linear functionals are
dependent.  Equivalently, for some `a`, the vector `T(ell)-a*ell` annihilates
`H`.  Since `H^perp=span(1)`, this holds exactly when

```text
T(ell) = a*ell + b*1.
```

For nonconstant `ell` the rank is one; otherwise it is two.  On `A_s`, an
affine transport forces the paired output law `y=a*x+b*s`.  It therefore
cannot produce an independent second row: an unaligned `(d,1)` target is
impossible, while an aligned target determines `d=(1-b*s)/a` when `a` is
nonzero.  This is the affine-log/direct-DLP boundary established in Round 39,
not a new index-calculus relation source.

S17's fixed eight-plus/nine-minus rows do not evade the theorem.  Differences
of sign patterns include `2*(e_i-e_j)`, and the odd P-256 subgroup order makes
two invertible.  These differences span the entire 16-dimensional sum-zero
hyperplane.  Both P-256 control fields reproduce rank 16 in both row orders.

The conclusion is intentionally scoped: it is a statement about the complete
linear span.  A proper nonlinear coefficient family could have a smaller
joint image without contradicting it, but the family must expose its retained
yield rather than treating rejected branches as searched.

## Complete toy census

The runner normalizes every nonzero vector in `F_5^4` projectively, preserving
all 156 profiles, and applies all 24 permutations.  For each pair it solves
the affine-transport condition, ranks the difference image in both row
orders, and enumerates all 125 rows in the fixed-sum affine slice.

| class | profile/permutation pairs | exact rank | expected image size | bucket occupancy |
|:--|--:|--:|--:|--:|
| constant controls | 24 | 0 | 1 | 125 |
| nonconstant affine | 600 | 1 | 5 | 25 |
| non-affine | 3,120 | 2 | 25 | 5 |

The separate six-pattern `2+2` sign control has difference-span rank three,
the complete sum-zero dimension at width four.  In total the runner performed
3,744 affine solves, 7,488 image-rank computations, 468,000 affine-row
enumerations, and 468,000 bucket insertions.  All expected image sizes and
occupancies match exactly.

## P-256 controls

Deterministic distinct scalar labels were generated independently modulo
`2^61-1` and the P-256 subgroup order.  Every matrix was ranked in forward
and reversed row order.

| certificate | rank in both fields | interpretation |
|:--|--:|:--|
| S17 fixed-sign difference span | 16 | the sign domain spans the full theorem space |
| half-turn paired image | 2 | non-affine transport does not collapse the joint image |
| identity paired image | 1 | affine transport collapses to one coordinate |
| planted `(Q,G)` augmented rows | 2 | information rank is real when a product hit is supplied |

The two exact group replays used 18 constructed points, 52 scalar
multiplications, 13,252 scalar bits, and 34 point additions.  All points were
on curve and scalar and group replay failures were zero.  The labelled
half-turn control introduces no auxiliary log vector, but its planted
coefficients are not a construction on the actual factor base.

## Cost boundary and next obligation

Round 310 changes the admissible mechanism set, not the executable cost
boundary.  It imports Round 309's hash-pinned projection without reclassifying
it as measurement:

| quantity | imported value |
|:--|--:|
| S17 product-target mean | `2^-271.270334` |
| first arity covering the product group | 40 |
| balanced arity-40 list width | `2^279.001098` |
| balanced-baseline ratio to rho | `2^151.122587` |
| structured degree | unset |

The remaining route is narrower:

> Construct a nonlinear restricted coefficient family whose searchable joint
> image is near size `n`, whose two augmented rows remain independent, whose
> retained yield is exhaustively accounted, and whose structured residual
> degree is at most five.

Without such a construction, additional tuning of an unrestricted column
permutation cannot reach rho parity.

## Isolation and artifact integrity

The protocol was committed as `1079b75ef` and the implementation as
`19e30de16`.  The canonical and independent runs were uncontended and emitted
byte-identical result and assessment files.

| run | exit | wall | user | peak RSS | contention |
|:--|--:|--:|--:|--:|:--|
| canonical-v1 | 0 | 0.132935 s | 0.132802 s | 2,228 KiB | none |
| independent-v2 | 0 | 0.149722 s | 0.149488 s | 2,048 KiB | none |

| artifact | bytes | SHA-256 |
|:--|--:|:--|
| `permutation-dichotomy-result.json` | 8,301 | `f8c39500b97394e633dcbcd36a92bf2cb896176aff7e4ffcb357e201b4dbfdc8` |
| `transfer-assessment.json` | 3,623 | `2dfe0762cc8e460f86ad8e24a45b223a9ccc88d7c6529a24dd21e245a4c80afe` |
| `isolation.jsonl` | 4,804 | `d40f494543628aa000eb70de6f779315c89b67a3deb5de9581bbdd953b1ba8d4` |
| `PROTOCOL.md` | 7,066 | `d7ca0e4433304386efc3304eb3ddef3e003e3a009a9c56ecd1c2273288008bda` |

Result semantic evidence:
`020616e5db516d7643ba50847726b3ee1ca8eb71d5ea84673591159dcbf9f3fa`.
Transfer-assessment semantic evidence:
`332fb3eb73e3a8f0e7e890215f037fcb7872f000d35a999a768a93fc5c7f61b8`.

The transfer workflow's methodology and assessment-template companion
resources were unavailable, and that limitation is recorded in the
assessment.

## Reproduction

```bash
cargo test --bin p256_permutation_dichotomy
cargo clippy --bin p256_permutation_dichotomy -- -D warnings
cargo build --release --bin p256_permutation_dichotomy --bin isolated_bench
target/release/isolated_bench run --wait --cpus 4 \
  --label p256-permutation-dichotomy-round310-canonical-v1 \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_permutation_dichotomy_round310_20261008/isolation.jsonl -- \
  target/release/p256_permutation_dichotomy \
  --round39 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_affine_log_orbit_round39_20261006/affine-log-orbit-result.json \
  --round304 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_algebraic_quotient_rigidity_round304_20261008/algebraic-quotient-rigidity-result.json \
  --round309 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_dickson_permutation_round309_20261008/dickson-permutation-result.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_permutation_dichotomy_round310_20261008/permutation-dichotomy-result.json \
  --assessment research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_permutation_dichotomy_round310_20261008/transfer-assessment.json
```
