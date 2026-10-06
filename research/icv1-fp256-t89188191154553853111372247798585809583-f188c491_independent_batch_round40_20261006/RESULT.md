# P-256 independent-start batch routing, round 40: result

Date run: 2026-10-06

Independent batching is the first screened reuse mechanism that preserves all
start entropy and makes the **pair-construction term alone** touch the frozen
rho boundary.  It does not establish parity end to end.

The construction threshold is a batch of 103,352,459,747 independent starts
(`2^36.589`).  Its 826,819,677,976 pair requests cover essentially the entire
34,562,148,612-entry signed-pair universe, reducing expected distinct pair
constructions to `0.334410508421` per start.  That exactly consumes the frozen
`0.334410508423`-addition allowance.  It leaves no routing budget while
requiring 28.11 TB of materialised state and 132.29 TB of logical radix traffic
for one batch.

At the largest batch below the registered `2^50`-byte storage gate, pair
construction falls to `0.008349680434` additions per start and the optimistic
arithmetic ratio is `0.998604763014` rho.  But that row materialises almost one
petabyte, moves 5.30 PB per batch, and requires routing throughput equivalent
to 3,925.65 bytes per native group-addition time.  No measurement was made in
that storage tier.

The actual complete RAM router at `2^16` starts costs 3.43258 native-addition
equivalents per start—10.2646 times the permitted headroom—and gives the
stage projection `1.013257282` rho.  It is therefore not promoted.

## Scope correction

The initial protocol named `FB1h2f8621cda105` as the measured base.  Before
implementation or execution, amendment 1 corrected that invalid composition:

- rounds 33–38's near-rho constants belong to the known-log low-delta/scalar
  selector family;
- round 19's structured residual maximum of 4 belongs only to
  `FB1h2f8621cda105`;
- this round measures batching for the former and retains the latter only as a
  comparison boundary.

Accordingly, no sub-rho number here is a projection for collecting 138,031
Dickson-base rows.  The Dickson relation-collection cost, per-row cost, and
sparse linear algebra remain unset.

## Exact depth ladder

Every start contributes all eight requests.  No request, start, or
probabilistic branch is discarded.

| starts | requests | distinct | repeats | observed / expected distinct | logical traffic | peak RAM | FP / FN |
|--:|--:|--:|--:|--:|--:|--:|--:|
| `2^12` | 32,768 | 32,768 | 0 | 1.000000474 | 5.24 MB | 1.11 MB | 0 / 0 |
| `2^14` | 131,072 | 131,072 | 0 | 1.000001896 | 20.97 MB | 4.46 MB | 0 / 0 |
| `2^16` | 524,288 | 524,283 | 5 | 0.999998048 | 83.89 MB | 17.83 MB | 0 / 0 |
| `2^18` | 2,097,152 | 2,097,078 | 74 | 0.999995052 | 335.54 MB | 71.30 MB | 0 / 0 |
| `2^20` | 8,388,608 | 8,387,565 | 1,043 | 0.999997010 | 1.342 GB | 285.21 MB | 0 / 0 |

The observed occupancy agrees with the exact uniform expectation at every
registered depth.  Memory and logical traffic are exactly linear in starts;
reuse is negligible at tractable depths because the pair universe is 34.56
billion entries.

The fixed traffic accounting is 160 bytes per request, or 1,280 bytes per
start:

```text
record generation write        16 bytes
three radix read/write passes  96 bytes
sorted scan read                16 bytes
accumulator read/write          32 bytes
```

Nine exhaustive toy cells over pair universes 17, 31, and 63 matched a direct
ordered-set reference exactly, with zero route failures, false positives, or
false negatives.

## Native group controls

For 4,096 deterministic routed pair identifiers, the binary decoded both
endpoints and signs, constructed both endpoint points in native P-256 group
arithmetic, and compared their signed sum with an independently combined
scalar construction.

```text
requests checked              4,096
decode failures                   0
group replay failures             0
negative-control false positives  0
false negatives                   0
scalar multiplications        12,288
group additions            1,589,154
group doublings            3,133,368
```

Control digest:
`98c8f56a06257256bb88f6c6d7cffce79391c461170b91287f8fd193fe5711e6`.
These controls validate routing and signed-pair decoding on synthetic
scalar-derived endpoints.  Rounds 34–36 remain the evidence for the real
factor-base pair representation.

## Isolated RAM timing

The release binary ran with one thread pinned to CPU 4 on an AMD EPYC 9V74.
The isolation wrapper moved 365 threads, observed 0.01 seconds of other CPU
work during the final 17.47-second run, and marked the run uncontended.  Peak
RSS was 416,896 KiB.

The binary performed a warm-up, an A/A control, and five interleaved
baseline/candidate repetitions.  The A/A ratio was 1.01957.

| quantity | median |
|:--|--:|
| native P-256 addition | 899.104 ns |
| complete RAM route per start | 3,086.245 ns |
| complete route / native addition | 3.432580 |
| allowed additions per start | 0.334411 |
| measured / allowed | **10.264569** |
| measured logical traffic rate | 414.744 MB/s |
| stage projection including measured router | **1.013257282 rho** |

The candidate timing includes SHA-256 request generation, allocation, three
radix passes, and reconstruction.  Its checksums were identical in every
repetition.  It is a matched RAM diagnostic at `2^16` starts, not a disk
measurement and not a completed selector or DLP.

Both isolated runs are retained in `isolation.jsonl`; the first successful
artifact before adding derived ratio fields is preserved under
`runs/prederived-result.json` rather than overwritten.

## Projected boundary

| batch | starts | distinct pairs / start | pair-only / rho | peak materialisation | traffic / batch | required routing bytes / addition-time |
|:--|--:|--:|--:|--:|--:|--:|
| `2^34` | 17.18 B | 1.974061 | 1.007016180 | 4.67 TB | 21.99 TB | no positive headroom |
| `2^36` | 68.72 B | 0.502945 | 1.000721173 | 18.69 TB | 87.96 TB | no positive headroom |
| parity threshold | 103.35 B | 0.334411 | 1.000000000 | 28.11 TB | 132.29 TB | effectively infinite |
| `2^38` | 274.88 B | 0.125736 | 0.999107069 | 74.77 TB | 351.84 TB | 6,133.97 B/add-time |
| `2^40` | 1.100 T | 0.031434 | 0.998703543 | 299.07 TB | 1.407 PB | 4,224.75 B/add-time |
| largest below `2^50` | 4.139 T | 0.008350 | **0.998604763** | 1.126 PB | 5.298 PB | 3,925.65 B/add-time |

The table's sub-rho rows are arithmetic-only projections.  Every projected
batch is far beyond the measured 17.83-MB timing cell, and no external-memory
throughput is imported or guessed.  The full known-log one-target stream is
about `2^128.326` starts and `2^138.648` logical routing bytes.  Those figures
do not apply to a 138,031-row Dickson collection.

## Gate status

| requirement | status | evidence or gap |
|:--|:--|:--|
| zero false positives and false negatives | passed | all exhaustive toys, five ladder cells, regenerated accumulators, and native controls |
| exact group replay | passed for routing controls | 4,096/4,096 native synthetic endpoint checks |
| structured residual degree at most 5 | failed for measured family | unknown for round-33 family; degree 4 belongs only to comparison base |
| relation collection below `2^120` | failed / unset | no valid composition with 138,031 Dickson rows |
| cost per usable relation below `2^103` | failed / unset | same composition gap |
| storage below `2^50` | passed narrowly in projection | selected row is one byte below the registered limit |
| measured complete time below rho in applicable tier | failed | RAM router is 10.26× over budget; projected tier unmeasured |
| no discarded branch counted as exhaustive | passed | every request and start is retained |

## Decision

Do not promote and do not attempt a full-depth relation.  Independent batching
survives as a pair-construction idea because it preserves start entropy and its
pair-only arithmetic can fall below rho.  The complete measured router does
not: it moves the RAM-stage projection to `1.013257282` rho, and the only
arithmetic-only sub-rho rows require tens of terabytes to one petabyte of
materialisation with unmeasured external routing.

The next bounded iteration is a two-pass 18-bit radix router on a pre-generated
request stream.  It should separate routing from the SHA-256 uniform-stream
control and reduce logical traffic from 160 to 128 bytes per request.  It must
retain the complete-router row as the baseline and may not transfer any result
to `FB1h2f8621cda105` without a new composition proof.

## Reproduction

```bash
cargo test --bin p256_independent_batch
cargo clippy --bin p256_independent_batch -- -D warnings -A clippy::mismatched-bit-width-type
cargo build --release --bin p256_independent_batch
RAYON_NUM_THREADS=1 python3 tools/isolated_bench.py run --wait --cpus 4 \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_independent_batch_round40_20261006/isolation.jsonl \
  --label round40-independent-batch-derived -- \
  target/release/p256_independent_batch \
  --round19 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_s17_multilevel_selector_round19_20261005/selector-result.json \
  --round36 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_mont_affine_round36_20261006/mont-affine-result.json \
  --round38 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_codebook_support_round38_20261006/codebook-support-result.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_independent_batch_round40_20261006/independent-batch-result.json
```

Final result artifact: 18,115 bytes, SHA-256
`52154eb455ea7edd187d83944ba4b6d953d92014d4f05919e4bd8b7bc097ecbf`.
Semantic evidence SHA-256:
`382b3caab85c1bd71eacf6f12bd37b0117fd9fb5c369d37ae2ed37fb719c6e49`.
Isolation JSONL SHA-256:
`c04812b438ea26ed6b0f94fee6a7425730fd49b395c9f92c70db0c61915b59a9`.

