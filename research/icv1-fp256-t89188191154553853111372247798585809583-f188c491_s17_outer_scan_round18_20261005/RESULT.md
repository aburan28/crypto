# P-256 charged S17 outer-relation scan round 18: result

Date run: 2026-10-05

**CORRECTNESS PASS; ATTACK-PROMOTION FAIL.**  The charged outer scanner agrees
with an independent full-point reference on all 4,096 actual P-256 oracle
calls.  It recovers and directly verifies the frozen length-17 planted
relation, reports zero false negatives and zero false positives, and reports
no relation for the hash-public target in the bounded prefix.

The newly priced enumeration boundary rejects this fixed-atom route.  A full
signed seventeenth-column scan around one sixteen-column atom projects to
`2^33.280` field multiplications, while one such atom has only
`2^-221.996` modeled relation success.  Even the explicitly optimistic
geometric projection is therefore `2^255.276` field multiplications per
relation and `2^272.404` for a 105%-sized relation matrix: **144.404 bits above
the `2^128` rho reference**, before charging atom generation or failed-target
exhaustion.

Low degree survives but does not rescue the search.  The local batched image
oracle remains degree 2, and the independent frozen residual systems remain at
degrees 3, 3, and 4 for residual depths 1, 2, and 3.  The degree-at-most-5 gate
passes; the below-`2^120` enumeration gate fails.

## Measured outer-prefix boundary

The `hash-0` packet is fixed at 36 bytes.  Its two eight-leaf images are built
once for 52,010 field multiplications.  Each eligible outer column contributes
two signed target adjustments and two exact batched oracle calls per target.

| outer depth | stream positions | public oracle calls | public query muls | build + query muls | planted relations | public relations | exact |
|---:|---:|---:|---:|---:|---:|---:|:--:|
| 4 | 16 | 32 | 1,269,728 | 1,321,738 | 1 | 0 | yes |
| 5 | 32 | 64 | 2,539,456 | 2,591,466 | 1 | 0 | yes |
| 6 | 64 | 128 | 5,078,912 | 5,130,922 | 1 | 0 | yes |
| 7 | 128 | 256 | 10,157,824 | 10,209,834 | 1 | 0 | yes |
| 8 | 256 | 512 | 20,315,648 | 20,367,658 | 1 | 0 | yes |
| 9 | 512 | 1,024 | 40,631,296 | 40,683,306 | 1 | 0 | yes |
| 10 | 1,024 | 2,048 | 81,262,592 | 81,314,602 | 1 | 0 | yes |

No packet column occurs in the first 1,024 permutation positions.  Every query
returns 256 roots from 128 quadratic solves, uses one batch of 128 inversions,
and charges exactly 39,679 multiplications.  The oracle-call slope is exactly
1.000 bits of work per added outer depth.  Including the one-time image build,
the fitted multiplication slope is 0.991543.

The planted witness occurs at the preregistered stream position 7, column
123,427, positive outer sign, and atom sign mask zero.  Direct addition of all
17 signed points reproduces the target.  The public target is SHA-256-selected
and is not constructed from factor-base points; its zero hits are a bounded
observation, not a claim that it has no relation in the complete factor base.

## Independent correctness boundary

The reference enumerates all 65,536 signed affine sums of the sixteen packet
points and indexes complete `(x,y)` points, not image x-coordinates.  It charges
131,070 group additions to build.  Each adjusted public or planted target is
looked up in that index independently of the quadratic oracle.  Candidate and
reference classifications agree on all 4,096 trials.  The sole emitted
relation is then replayed with 17 fresh group additions.

At depth 10, per target:

| quantity | value |
|:--|---:|
| signed branch cells / oracle calls | 2,048 / 2,048 |
| quadratic solves | 262,144 |
| roots / index lookups | 524,288 / 524,288 |
| batch inversions / scalar fallbacks | 2,048 / 0 |
| field multiplications | 81,262,592 |
| target-adjustment group additions | 2,048 |
| false negatives / false positives | 0 / 0 |

The planted target has one candidate and reference hit.  The hash-public target
has zero in both paths.

## Enumeration projection

The selected factor base has 131,458 columns.  There are

```text
C(131458,16) =
379818644896935897531886053293722324029315526042580642955413667491880
```

distinct sixteen-column atoms (`2^227.816`) and `2^240.733` distinct
seventeen-column sets.  Including 17 independent signs gives a complete signed
domain of `2^257.733`, Poisson mean 3.3242414 and the previously reported
96.400018% whole-base relation probability.

That high global probability is not the probability for one fixed atom.  One
complete atom scan has 131,442 eligible outer columns and 262,884 signed oracle
calls.  Its signed candidate domain is only 17,228,365,824, giving:

| projected boundary | value | log2 |
|:--|--:|--:|
| success probability per complete atom scan | `1.4878707e-67` | -221.996 |
| expected atom scans per relation | `6.7210141e66` | 221.996 |
| complete atom-scan multiplications | `1.0431026e10` | 33.280 |
| multiplications per relation | `7.0107074e76` | **255.276** |
| targets for 138,031 rows at 96.400018% | 143,185.66 | 17.127 |
| collector multiplications | `1.0038328e82` | **272.404** |

This is deliberately optimistic: it treats atom scans as independent
geometric trials and charges neither atom generation nor the complete
exhaustion of the approximately 3.6% of targets with no relation.  Those
omissions can only increase the fixed-atom route's work.

## Memory and sparse linear algebra

| item | logical bytes |
|:--|---:|
| persistent atom packet | 36 |
| reconstructed half-images | 8,192 |
| candidate peak including batch scratch | 75,008 |
| complete factor-base points, 96 raw bytes per column | 12,619,968 |
| independent-reference final table | 4,521,984 |
| independent-reference build peak, excluded from candidate | 13,074,432 |

For 138,031 rows of weight 17, the sparse matrix contains 2,346,527
nonzeros.  The frozen CSR model uses 12,836,891 bytes; three 32-byte field
vectors use 12,619,968 bytes, for 25,456,859 modeled working bytes
(24.278 MiB).  A conservative `2*columns` Wiedemann schedule charges
616,939,492,732 sparse nonzero additions (`2^39.166`), and the separate
Berlekamp--Massey model charges 17,281,205,764 scalar operations (`2^34.008`).
These units are not added to field multiplications.  Both are negligible next
to the `2^272.404` collection projection; relation acquisition, not matrix
storage or sparse linear algebra, is the rejection boundary.

## Reproduction and evidence

```bash
cargo test --bin p256_s3_image_transfer
cargo clippy --bin p256_s3_image_transfer -- -D warnings \
  -A clippy::mismatched_bit_width_type
cargo build --release --bin p256_s3_image_transfer
target/release/p256_s3_image_transfer \
  --outer-scan-round17 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_s17_batch_inversion_round17_20261005/batch-result.json \
  --round6 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_dickson_residual_scaling_round6_20261004/degree-result.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_s17_outer_scan_round18_20261005/outer-result.json
```

`outer-result.json` is 25,338 bytes with SHA-256
`b9081ea3cd3957dcae7805852371e395a2a63c8d4a2062e5c8e2c0b23980fce3`.
Both dependency hashes, the complete factor-base receipt, the 36-byte packet,
and the residual-degree envelope are checked before scanning.

## Decision and next credible lever

Reject enumeration of fixed sixteen-column atoms as a P-256 S17 relation
collector.  Improving the 39,679-multiplication microkernel by another constant
factor cannot close a 144.404-bit gap over rho.

Any follow-on must change the atom-selection exponent: it must propagate or
index compatibility across variable column sets without visiting roughly
`2^222` fixed atoms per relation.  A proposal that retains fixed-atom
enumeration is stopped by this result.  The measured local degree 2 and
residual degrees 3--4 remain useful constraints for such a different global
solver, but they are not an end-to-end degree-of-regularity result for the
unsplit P-256 S17 ideal.
