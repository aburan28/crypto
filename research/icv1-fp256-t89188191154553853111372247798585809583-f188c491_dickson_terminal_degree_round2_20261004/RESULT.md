# P-256 Dickson terminal-fibre degree round 2: result

Date run: 2026-10-04

The terminal-fibre hypothesis is **falsified**.  Every one of the 48 eligible
full fibres—31 at depth 4 and 17 at depth 5—had paired common-positive `S3`
solving degree 4, exactly the terminal-zero reference.  All 288 F4 runs
(candidate plus paired reference) completed and agreed with exhaustive
enumeration.  Changing the Dickson terminal constant does not lower the
observed degree barrier.

The P-256 coset search did find a larger factor base.  Hash-derived coset 4
gives 131,458 negation-folded columns, 219 more than the accepted terminal-zero
base and 18 more than round 1's affine-bitbox candidate.  It is the best
factor-base inventory found so far, but it has no lower-degree claim.

## Degree sweep

Every fibre has exactly `2^t` roots before curve lifting.  `targets` is the
number of common positive affine targets available; the frozen first three
were measured.  Every row is 3/3 correct and complete.  `degree` is
min/median/max, and columns and field operations are medians to the last
productive F4 step.  The degree ratio to the paired `c=0` reference is 1.00
on every row.

| depth | terminal | liftable roots | targets | degree | columns | field ops | classification |
|---:|---:|---:|---:|:--|---:|---:|:--|
| 4 | 0 | 8 | 124 | 4/4/4 | 309 | 1,183,930 | reference |
| 4 | 28 | 11 | 22 | 4/4/4 | 309 | 1,105,856 | engineering only |
| 4 | 48 | 9 | 12 | 4/4/4 | 309 | 1,182,339 | equal degree |
| 4 | 64 | 11 | 16 | 4/4/4 | 309 | 1,183,559 | equal degree |
| 4 | 109 | 8 | 10 | 4/4/4 | 309 | 1,202,409 | regression |
| 4 | 165 | 8 | 12 | 4/4/4 | 309 | 1,183,309 | equal degree |
| 4 | 173 | 8 | 8 | 4/4/4 | 309 | 1,183,645 | equal degree |
| 4 | 294 | 11 | 20 | 4/4/4 | 309 | 1,183,709 | equal degree |
| 4 | 341 | 11 | 26 | 4/4/4 | 309 | 1,176,204 | engineering only |
| 4 | 358 | 9 | 10 | 4/4/4 | 309 | 1,687,256 | regression |
| 4 | 369 | 8 | 10 | 4/4/4 | 309 | 1,175,822 | engineering only |
| 4 | 383 | 8 | 18 | 4/4/4 | 309 | 1,183,132 | equal degree |
| 4 | 401 | 8 | 14 | 4/4/4 | 309 | 1,183,853 | equal degree |
| 4 | 428 | 9 | 16 | 4/4/4 | 309 | 1,188,021 | regression |
| 4 | 510 | 8 | 10 | 4/4/4 | 309 | 1,202,347 | regression |
| 4 | 548 | 8 | 4 | 4/4/4 | 309 | 1,202,059 | regression |
| 4 | 603 | 8 | 8 | 4/4/4 | 309 | 1,183,933 | equal degree |
| 4 | 641 | 8 | 10 | 4/4/4 | 309 | 1,202,454 | regression |
| 4 | 675 | 9 | 22 | 4/4/4 | 309 | 1,183,821 | equal degree |
| 4 | 750 | 8 | 20 | 4/4/4 | 309 | 1,202,956 | regression |
| 4 | 768 | 8 | 8 | 4/4/4 | 309 | 1,202,673 | regression |
| 4 | 782 | 11 | 24 | 4/4/4 | 309 | 1,105,052 | frozen winner; engineering only |
| 4 | 793 | 11 | 20 | 4/4/4 | 309 | 1,202,829 | regression |
| 4 | 810 | 8 | 22 | 4/4/4 | 309 | 1,183,227 | equal degree |
| 4 | 857 | 9 | 24 | 4/4/4 | 309 | 1,183,268 | equal degree |
| 4 | 978 | 8 | 16 | 4/4/4 | 309 | 1,105,380 | engineering only |
| 4 | 986 | 8 | 8 | 4/4/4 | 309 | 1,183,240 | equal degree |
| 4 | 1042 | 8 | 18 | 4/4/4 | 309 | 1,152,085 | engineering only |
| 4 | 1087 | 9 | 22 | 4/4/4 | 309 | 1,183,728 | equal degree |
| 4 | 1123 | 8 | 8 | 4/4/4 | 309 | 1,686,445 | regression |
| 4 | 1150 | 14 | 60 | 4/4/4 | 309 | 1,183,731 | cardinality only |
| 5 | 0 | 13 | 304 | 4/4/4 | 647 | 36,725,681 | reference |
| 5 | 1 | 16 | 84 | 4/4/4 | 647 | 36,202,665 | engineering only |
| 5 | 28 | 19 | 138 | 4/4/4 | 647 | 36,176,809 | engineering only |
| 5 | 109 | 20 | 178 | 4/4/4 | 647 | 36,177,187 | engineering only |
| 5 | 173 | 13 | 64 | 4/4/4 | 647 | 36,223,598 | engineering only |
| 5 | 341 | 19 | 150 | 4/4/4 | 647 | 32,409,965 | engineering only |
| 5 | 369 | 16 | 138 | 4/4/4 | 642 | 31,566,344 | frozen winner; engineering only |
| 5 | 401 | 20 | 152 | 4/4/4 | 647 | 36,154,616 | engineering only |
| 5 | 510 | 16 | 114 | 4/4/4 | 647 | 36,175,307 | engineering only |
| 5 | 641 | 20 | 164 | 4/4/4 | 647 | 36,205,344 | engineering only |
| 5 | 750 | 16 | 110 | 4/4/4 | 647 | 36,161,928 | engineering only |
| 5 | 782 | 19 | 120 | 4/4/4 | 647 | 34,797,195 | engineering only |
| 5 | 810 | 16 | 98 | 4/4/4 | 647 | 31,440,737 | engineering only |
| 5 | 978 | 13 | 62 | 4/4/4 | 647 | 37,523,579 | regression |
| 5 | 1042 | 16 | 130 | 4/4/4 | 647 | 36,174,932 | engineering only |
| 5 | 1123 | 16 | 92 | 4/4/4 | 647 | 36,172,270 | engineering only |
| 5 | 1150 | 20 | 144 | 4/4/4 | 647 | 36,728,728 | cardinality only |

The frozen ranking selects terminal 782 at depth 4 and terminal 369 at depth
5 because it compares degree, matrix columns, field operations, lift count,
then terminal in that order.  Neither selection is a degree advance.

## Direct P-256 coset search

The eight root exponents are the low 78 SHA-256 bits specified in the
protocol, made odd.  Candidate 4 won the frozen cardinality rule:

| k | root exponent | columns | FB1 |
|---:|:--|---:|:--|
| 0 | `0x23e6cba8834be09f5b1d` | 131,343 | `FB1h6b7699a7d171` |
| 1 | `0x25ed4362f26d2637b345` | 131,034 | `FB1h085b2065fa5b` |
| 2 | `0x4575da009fcbf9ee1fd` | 131,136 | `FB1h6837e846b36c` |
| 3 | `0xa8c2d4b6e375f181469` | 130,966 | `FB1h4d26bec094de` |
| 4 | `0x2b6fdc73dc04e7667129` | 131,458 | `FB1h2f8621cda105` |
| 5 | `0x212c0b7da2df8d57187f` | 131,201 | `FB1hf3f83ccb1b6b` |
| 6 | `0x3f4e18fd61903c9c42bb` | 131,123 | `FB1h142ab37c877c` |
| 7 | `0x2830d05994353952ad75` | 130,955 | `FB1h0a4989de61a1` |

Selected `FB1h2f8621cda105` has full digest
`2f8621cda105a51f8703dabe38dec87f3715ab855770acd7b2af0a14906f9e42`,
262,916 signed points, terminal trace
`0x5b17195299a3158b93389ad04c776fff2a8bb23ca8659b5b0c2b75f9009b65e5`,
and point-key digest
`70ab21cbda76205ff43ab6bf85c68607be9f6b89df6399e6fa98e3489b888cf1`.
The full native rebuild and row comparison passed.

| base | columns | m=17 Poisson success | delta from terminal zero |
|:--|---:|---:|---:|
| selected nonzero Dickson coset | 131,458 | 96.400018207451% | +219 columns; +0.35049 percentage points |
| accepted terminal-zero Dickson | 131,239 | 96.049526975874% | reference |
| round-1 affine-bitbox | 131,440 | 96.372082631861% | +201 columns |

This makes the selected nonzero coset the best P-256 relation-domain count in
the current search, while leaving the measured toy solving degree at 4.

## Reproduction and evidence

```bash
cargo run --release --bin p256_dickson_terminal_degree -- \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_dickson_terminal_degree_round2_20261004/degree-result.json

cargo run --release --bin p256_dickson_coset_search -- \
  --curve icv1-fp256-t89188191154553853111372247798585809583-f188c491 \
  --verify \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_dickson_terminal_degree_round2_20261004/factor-base-result.json
```

| artifact | bytes | SHA-256 |
|:--|---:|:--|
| `degree-result.json` | 206,970 | `5c682c6d11feb72f86b74f8b1d5a4b92d0728c8631131a32b1d2361e30e21b37` |
| `factor-base-result.json` | 4,710 | `dcaa1f66f5f89757a3670a4a8e158ca6d5c360b7ca119d384d2c4a0b2794c2c8` |

The next degree iteration should stop tuning the factor-base constant.  The
data show that the quartic `S3` input itself pins solving degree 4 across a
broad set of full fibres.  The credible next target is a Dickson-specific
quadratic incidence or bilinearized addition chain whose degree stays below 4
at depths 4 and 5, with matrix width and field operations reported beside the
degree.
