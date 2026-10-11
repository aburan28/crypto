# P-256 Dickson image-atomized S5 round 11: result

Date run: 2026-10-04

**PASS on every frozen degree and correctness gate.**  Solving each fixed leaf
pair's affine `S3` image first and atomizing the left image reduces all 152
final-join F4 systems to solving degree 2.  Every image atom completes and
agrees with direct set intersection, maximum Macaulay columns fall to 3, and
the affine-plus-identity union still matches all 16 parent cells with zero
false negatives and false positives across 1,664 signed tuples.

This is the first exact four-summand result in the sequence with a uniform
degree-2 F4 remainder.  The degree belongs to the image-atomized decomposition,
not to the unsplit balanced `S5` ideal.

## Boundary table

Round 10's exact leaf-atomized systems are the matched before rows.  Candidate
total field work below is image-build multiplications plus final-coefficient
multiplications plus F4 row-reduction operations.  Identity checks and signed
group additions remain separate units.

| terminal | parent cells | leaf atoms | image atoms | signed tuples | complete / correct | max solving degree | max columns | image-build muls | coefficient muls | F4 ops | identity checks | FN / FP | class |
|---:|---:|---:|---:|---:|:--|---:|---:|---:|---:|---:|---:|:--|:--|
| 0 | 8 | 8 | 8 | 128 | 8 / 8 | **2** | 3 | 585 | 88 | 130 | 9 | 0 / 0 | degree advance |
| 369 | 8 | 96 | 144 | 1,536 | 144 / 144 | **2** | 3 | 7,800 | 1,584 | 2,454 | 56 | 0 / 0 | degree advance |

All 152 image atoms have solving degree exactly 2.  The image split expands
104 leaf atoms to 152 final joins, a factor of 1.46154 overall and 1.5 on
terminal 369.  The original 16 parent cells therefore average 9.5 image atoms.

For terminal 369, candidate counted field work is
`7,800 + 1,584 + 2,454 = 11,838`, or 0.087091 of round 10's 135,927 F4
operations (11.48 times lower even before converting units).  F4 row reduction
alone falls 55.39 times, from 135,927 to 2,454 operations.  For terminal 0,
the corresponding total is 803 versus 7,768, a ratio of 0.103373.

## Exact image and chart checks

For every fixed leaf pair, the algebraic quadratic roots are compared with all
signed point sums before use.  There are no missing or surplus affine image
roots.  Each left image root then creates one single-variable F4 system:

```text
S3(zL,zR,target_x) = 0
product(zR-r for r in right_image) = 0.
```

Both polynomials have degree at most 2.  F4 consistency matches direct right-
image membership on all 152 systems.  Round 10's 65 direct identity checks are
carried unchanged and continue to recover the three decompositions whose pair
intermediate is infinity.

## Reproduction and evidence

```bash
cargo test --bin p256_dickson_balanced_s5
cargo clippy --bin p256_dickson_balanced_s5 -- -D warnings \
  -A clippy::mismatched_bit_width_type
cargo run --release --bin p256_dickson_balanced_s5 -- \
  --image-atomize-round10 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_atomized_s5_round10_20261004/degree-result.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_image_atomized_s5_round11_20261004/degree-result.json
```

The binary hash-checks round 10 before solving and independently verifies every
pair image with the group law.  Its round-9 and round-10 modes remain
byte-identical after the round-11 extension.

`degree-result.json` is 179,215 bytes with SHA-256
`7266c05093c815734db0a0319b946f76d003e836cce0911541be5ab4d50f2e60`.
An immediate replay was byte-identical.

## Decision and limit

Accept image atomization plus explicit identity charts as an exact degree-2
decomposition for these frozen toy four-summand components.  This achieves the
requested regularity reduction at the local four-summand boundary: the
curve-lift and unsplit specialised systems reached degree 3, while every final
image join is degree 2.

The exchange is explicit component growth.  The next credible experiment must
measure image-set and atom growth when two exact four-summand blocks are joined
into eight summands, retaining infinity sentinels and exhaustive correctness
on a bounded corpus.  A degree-2 local join is not useful if image counts
restore the exponential frontier.

No actual P-256 `S18` degree, length-17 relation, logarithm, end-to-end S, or
rho speedup is established.
