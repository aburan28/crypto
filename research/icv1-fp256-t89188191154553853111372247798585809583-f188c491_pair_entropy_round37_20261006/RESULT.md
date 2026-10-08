# P-256 cutoff-conditioned pair entropy, round 37: result

Date run: 2026-10-06

Capacity conditioning does not concentrate signed pair accesses enough for a
small hot tier.  At the largest registered depth, a perfect 256-MiB oracle
table contains only `0.128826` signed pairs expected to be present among all
136 pairs of a retained 17-column start.  Even allowing the selector to choose
its eight pairs adaptively after seeing the whole start, expected cache hits
are at most `0.128826`; at least `7.871174` cold pair constructions remain.

Granting every hot access zero cost, every cold construction exactly one
addition, and ignoring factor-base reads gives an optimistic lower bound of
`1.032250` rho.  Merely lowering this bound to parity requires at least
272,190,818 signed entries, or 17,420,212,352 bytes (`2^34.020`).  That
threshold is only necessary: the calculation ignores matching conflicts, so
it does not establish that a realizable selector could obtain those hits.

This exact negative result closes capacity-only hot tiering at the registered
depths.  No hot/cold benchmark and no full-depth unplanted relation were
attempted.

## Frozen boundary

The audit imports round 33 by SHA-256
`9931bccbd9e2821f4498f465ce65f83efb400e7fb3740a239bdb40d3275ffbd8`
and round 36 by SHA-256
`f145231c098da83878e890d340f8ae921598519b0f2626bcc7c4cd8f9446d015`.
It preserves:

```text
factor-base columns                       131,458
variable columns per start                17
cutoff capacity                            219
signed pair lookups per segment             8
signed pair entries             34,562,148,612
Montgomery-affine entry                     64 bytes
corrected base ratio          0.998569034150286 rho
mean segment capacity         224.361249307511
remaining parity budget         0.334410508423 additions
```

The capacity histogram remains exactly 6,935 columns in every class 0 through
17 and 6,628 columns in class 18.  Recomputing its without-replacement
generating function gives exactly round 33's
`4,299,882,272,811,428,257,776,224,473,756,934,611,512,306,612,347,523,812,324,992,780,973,981`
retained starts.

## Optimistic oracle frontier

For every one of the 34,562,148,612 signed entries, the audit computes its
exact conditional probability of appearing among the 136 signed pairs of a
retained start.  It then fills each table with the highest-probability entries.
The hit column below is an upper bound because an actual selector must choose
eight disjoint pairs; the oracle only counts cached edges present and ignores
all matching conflicts.

| table | signed entries | bytes | expected cached pairs present / hit upper bound | cold constructions lower bound | projected / rho lower bound |
|:--|--:|--:|--:|--:|--:|
| 2 MiB | 32,768 | 2,097,152 | 0.001006 | 7.998994 | 1.032797 |
| 64 MiB | 1,048,576 | 67,108,864 | 0.032207 | 7.967793 | 1.032664 |
| **256 MiB** | **4,194,304** | **268,435,456** | **0.128826** | **7.871174** | **1.032250** |
| necessary parity threshold | 272,190,818 | 17,420,212,352 | 7.665590 | 0.334410 | 1.000000 |
| complete table | 34,562,148,612 | 2,211,977,511,168 | 136.000000 / 8.000000 | 0.000000 | 0.998569 |

The 17.42-GB threshold is 0.787540% of the complete table.  It is far above
the registered cache-resident depths and assumes free hot access.  Its
`0.999999999911` numerical lower bound reflects the exact integer ceiling at
the frozen 12-decimal budget; it is not a measured or achievable sub-rho row.

The densest block is capacity `(18,18)`.  It alone contains 87,847,512 signed
entries (5.62 GiB), and each entry is present with conditional probability
only `3.071454e-8`.  The cutoff biases classes, but not individual pairs within
a class enough to make a small pair cache effective.

## Exact controls

```text
round-33 retained count match          true
full tuple population match            true
full expected signed-pair count         136 exactly
capacity-class blocks                   190
signed entry population      34,562,148,612
density monotonicity failures              0
population failures                        0
toy exhaustive cells                       2
toy fixed pairs checked                    43
toy count failures                         0
false positives                            0
false negatives                            0
```

Deterministic digests:

```text
ordered class rows  97c1524e0528fb2a196d418382fd25838331f44db4088d49797dce4a4e2e3ab5
frontier            0ca2c6224a29961da1a84a20fb3572706025301bf8524b91a2087f03c8e2d505
semantic            60a9284e34c5a39e01bb422badd1a92489f961bed5dda353c6e8d3c59f0f9fab
```

An independent release replay reproduced the complete artifact byte for byte.

## Requirement status

| requirement | status | evidence or gap |
|:--|:--|:--|
| exact conditional signed-pair distribution | verified | all 190 blocks; exact round-33 count and 136-pair normalization; zero toy failures |
| registered hot tier at or below rho | failed | 256-MiB optimistic lower bound is `1.032250` rho |
| complete measured time at or below rho | failed | analytical screen fails before implementation |
| storage below `2^50` | verified for the necessary threshold | 17.42 GB is `2^34.020` bytes, but is not cache resident or sufficient |
| proved P-256 usable-relation probability | not established | pair co-presence does not prove relation coverage |
| structured residual degree at most 5 | not established | no algebraic residual solve in this round |
| collection below `2^120` and per row below `2^103` | not established | no usable full-depth relation or collector |
| demonstrably non-generic end-to-end method | not established | imported restart segments remain representation-based |

## Decision

Reject capacity-only hot tiering at 2, 64 and 256 MiB.  Do not build the
hot/cold selector: even a stronger adaptive oracle cannot cover enough pairs.
A further iteration must correlate pair generation itself across segments or
change the underlying start distribution; another cache-size or prefetch
tuning round cannot bridge this bound.

## Reproduction

```bash
cargo test --bin p256_pair_entropy
cargo run --release --bin p256_pair_entropy -- \
  --round33 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_biased_restart_round33_20261006/biased-restart-result.json \
  --round36 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_mont_affine_round36_20261006/mont-affine-result.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_pair_entropy_round37_20261006/pair-entropy-result.json
```

Canonical result artifact: `pair-entropy-result.json`, 81,754 bytes, SHA-256
`8b86bbc1488e2c3f9d6bc8a6c4fc9cf0c8c21a1a440adcc8fecba5569c7ed687`.
