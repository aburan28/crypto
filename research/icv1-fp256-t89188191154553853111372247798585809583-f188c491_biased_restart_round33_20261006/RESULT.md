# P-256 biased common-edge restart selector, round 33: result

Date run: 2026-10-06

The registered biased selector crosses rho in group-addition accounting only.
At cutoff `D>=219`, the complete exact start distribution gives

```text
0.9942899645 rho  without the right-colour target correction,
0.9985690342 rho  with that correction.
```

This is the first sub-rho stage row in the P-256 factor-base workstream.  It is
not end-to-end parity.  The margin with the target correction is only 0.1431%;
the calculation grants every retained path visit as distinct P-256 support,
does not convert random reads from a 1.140-TB pair table into group-equivalent
work, has no proved P-256 relation probability or structured degree, and has
not established a non-generic algorithm.

No full-depth unplanted relation was attempted.

## Exact start distribution

Round 31's primitive `R=6935` mechanical cycle has the following residual
common-edge capacities:

```text
capacity 0..17: 6,935 columns each
capacity 18:    6,628 columns
```

The arbitrary-precision without-replacement generating function sums to

```text
binom(131458,17)
= 2936713077796649896669774506884320571474193610240522639490910781321628880
```

exactly.  Its distribution digest is
`bf0076371d24d569e7a73f94ee1e80d2215f77f87907190e1177d005a322066e`.

The selected cutoff retains 0.146418% of unsigned starts.  Their exact mean
capacity is 224.361249 common transitions.  Crediting all `2^17` sign choices
and all `D+1` states gives 256.133431 support bits, just 0.133431 bits above the
P-256 group order.  Increasing the cutoff to 220 drops the support upper bound
below the group and immediately imposes a 1.062584 target-retry lower bound.

| cutoff | retained starts | mean D | support upper bits | target retry | / rho, no target correction | / rho, with correction |
|--:|--:|--:|--:|--:|--:|--:|
| 0 | 100% | 152.643 | 256.000 capped | 1 | 1.008272 | 1.014548 |
| 180 | 11.7790% | 190.376 | 256.000 capped | 1 | 0.999609 | 1.004648 |
| 200 | 1.84124% | 207.330 | 256.000 capped | 1 | 0.996739 | 1.001368 |
| **219** | **0.146418%** | **224.361** | **256.000 capped** | **1** | **0.994290** | **0.998569** |
| 220 | 0.125115% | 225.274 | 255.912423 | 1.062584 | 1.056388 | 1.060917 |
| 230 | 0.022348% | 234.469 | 253.484837 | 5.716623 | 5.676630 | 5.700042 |

All 307 cutoffs were evaluated.  The ordered sweep digest is
`55bd716f5bf006282926613c75e5fa065860f5c9417b47b6866a8ef346d94586`.
Discarded starts are represented by the support upper bound and retry charge;
they are not reported as exhaustively searched.

## Pair-table accounting

The registered implementation materializes every signed pair of distinct
columns:

| item | value |
|:--|--:|
| pair entries | 34,562,148,612 |
| bytes per compressed point | 33 |
| materialized bytes | 1,140,550,904,196 |
| log2 materialized bytes | 40.052868 |
| table-build group additions | 34,562,148,612 |
| build / sqrt(n) | 1.016e-28 |
| lookups / segment | 8 |
| bytes read / segment | 264 |
| setup additions / segment | 8 |

The table passes the registered `2^50`-byte gate.  Eight lookups replace 16
input points by eight stored pair sums; combining those with the remaining
point costs eight additions.  The headline includes this setup and the full
table build.  The 0.998569 row additionally charges one right-colour target
addition per segment.  Neither row prices random table traffic in the same
unit as rho, so the complete measured-time gate remains false.

## Exact references and native replay

The two complete toy cells enumerate 4,965 unsigned tuples, 78,480 signed
segments, 965,664 visited path states and 887,184 exact transitions.  Every
cutoff agrees with the generating-function distribution and exact support is
never larger than the registered visit bound.  Edge failures, support-bound
violations, distribution failures, false positives and false negatives are
all zero.

The native replay reconstructs the complete 131,458-column P-256 coefficient
cycle and matches round 31's digest
`980917981827d813e60484abb0655e8bd527b0beb540ea146d2974ff17303a53`.
It then accepts 4,096 deterministic cutoff-219 starts from 2,803,238 charged
draws, rejecting 2,799,142.  Their measured mean capacity is 224.402100.  All
919,151 transitions replay in exact P-256-order coefficient arithmetic with
zero failures; 32,768 pair lookups read 1,081,344 projected bytes.  The sample
digest is
`52ffdecf9f4cf3343af24b6bd8c30f8082beb5ecd2479019c1362a4af05c9166`.

## Requirement status

| requirement | status | evidence or gap |
|:--|:--|:--|
| exact deterministic distribution and replay | verified | exact binomial identity, complete toy references and 919,151 native transitions |
| discarded starts charged | verified for the registered bound | retained visits cap support; cutoff 220 and above pay target retries |
| group-addition projection at or below rho | verified as a stage relaxation | 0.994290; 0.998569 with target correction |
| peak materialized storage below `2^50` | verified | `2^40.053` projected bytes |
| measured complete time including table traffic | not established | 1.140-TB table was not materialized or benchmarked |
| proved P-256 usable-relation probability | not established | 256.133 support bits are an occurrence upper bound, not achieved coverage |
| structured residual degree at most 5 | not established | imported degree remains unknown |
| collection and per-usable-relation gates | not established | no full relation or end-to-end collector |
| non-generic end-to-end algorithm | not established | segments are translations with representation restarts |

## Decision

Preserve cutoff 219 as a stage candidate; do not promote it as rho parity.
The next registered experiment must measure random pair-table traffic in the
same unit as P-256 group additions and replace the visit-count support upper
bound with exact or statistically bounded P-256 coverage.  The 0.1431% margin
with target correction is too small to absorb an unmeasured memory term.

## Reproduction

```bash
cargo test --bin p256_biased_restart_selector
cargo run --release --bin p256_biased_restart_selector -- \
  --round31 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_low_delta_round31_20261006/low-delta-result.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_biased_restart_round33_20261006/biased-restart-result.json
```

Canonical result artifact: `biased-restart-result.json`, 273,383 bytes,
SHA-256
`9931bccbd9e2821f4498f465ce65f83efb400e7fb3740a239bdb40d3275ffbd8`.
