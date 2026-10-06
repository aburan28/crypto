# P-256 fixed-universe direct-address replay, round 42: result

Date run: 2026-10-06

Direct addressing removes the catastrophic start-proportional storage term,
but exact replay is far too expensive to approach rho.

The full 34,562,148,612-entry signed-pair universe needs exactly
557,314,646,392 algorithmic bytes: 4,320,268,584 bytes of seen bitmap,
552,994,377,792 bytes of 16-byte pair values, and one 16-byte streaming
accumulator.  This is `2^39.020` bytes, versus round 41's 299-TB to 1.126-PB
start-proportional projections.

At the registered primary cell, however, complete direct-address replay took
1.00235 seconds, 2.496 times the complete materialised reference.  It cost
8.87777 native P-256 addition equivalents per start, or 26.55 times the
0.334411 routing allowance.  The direct candidate's fitted time exponent was
1.1066, worse than the reference's 1.0198.

Even an optimistic projection that transfers this scaled-table timing to the
full table remains 1.03656 to 1.03710 times rho.  Such transfer is not valid:
the largest measured table was 67.63 MB, while the projected table is 557.31
GB and has unmeasured random-access behavior.

## Preregistration corrections

Both corrections were committed before source implementation or execution.

- Amendment 01 removed an erroneous `16S` accumulator-vector term.  Replay
  visits each start's eight requests contiguously, so the candidate streams
  one accumulator at a time.
- Amendment 02 corrected the decimal evaluation of the already-correct state
  formula to 557,314,646,392 bytes.

The harness retains reference accumulators only for exhaustive comparison.
That memory is reported separately and is not counted as candidate state.

## Scope

This is a representation-stage experiment for the round-33 known-log
low-delta/scalar selector family.  It imports and verifies round 41's request
digest, native-control digest, exactness gates, isolation status, and
non-promotion decision.

`FB1h2f8621cda105` remains a separate comparison base with an observed
structured residual maximum of 4.  This experiment does not transfer that
degree, project a 138,031-row Dickson relation collection, price sparse linear
algebra, or solve a P-256 discrete logarithm.

## Exact algorithm

The matched materialised reference generates all eight requests per start
once, applies three stable 12-bit radix passes, constructs one synthetic
16-byte value per distinct pair, and reconstructs every start.

The candidate retains no request records:

1. generate every request and set its pair-universe bitmap bit;
2. on first occurrence, construct and store the pair's 16-byte value;
3. generate every request again;
4. direct-address all eight values for one start, emit its accumulator, and
   advance to the next start.

No request, start, or probabilistic branch is discarded.  The synthetic value
function checks routing semantics; projected native pair constructions are
charged separately.

## Complete checked cells

Nine exhaustive toy cells over universes 17, 31, and 63 matched direct set
references exactly.  The load-one ladder also matched every materialised
accumulator:

| universe | starts | requests | distinct observed / expected | candidate state | reference peak | semantic / cache-line traffic | FP / FN |
|--:|--:|--:|--:|--:|--:|--:|--:|
| `2^12` | 512 | 4,096 | 2,593 / 2,589.35 | 66.1 KB | 139.3 KB | 238.9 KB / 987.3 KB | 0 / 0 |
| `2^14` | 2,048 | 16,384 | 10,435 / 10,356.85 | 264.2 KB | 557.1 KB | 956.7 KB / 3.96 MB | 0 / 0 |
| `2^16` | 8,192 | 65,536 | 41,387 / 41,426.84 | 1.057 MB | 2.228 MB | 3.821 MB / 15.78 MB | 0 / 0 |
| `2^18` | 32,768 | 262,144 | 165,640 / 165,706.80 | 4.227 MB | 8.913 MB | 15.29 MB / 63.14 MB | 0 / 0 |
| `2^20` | 131,072 | 1,048,576 | 663,064 / 662,826.63 | 16.91 MB | 35.65 MB | 61.15 MB / 252.64 MB | 0 / 0 |
| `2^22` | 524,288 | 4,194,304 | 2,650,160 / 2,651,305.97 | 67.63 MB | 142.61 MB | 244.58 MB / 1.010 GB | 0 / 0 |

All 1,395,284 negative controls, all selected stored-value regenerations, and
all reference regenerations passed.  Candidate and reference BLAKE3
accumulator digests were identical in every cell.

At fixed universe `2^18`, the saturation ladder confirmed the state remains
4,227,088 bytes as request load rises from 1/4 to 16:

| load | starts | distinct / universe | candidate state | reference peak | semantic / cache-line traffic |
|--:|--:|--:|--:|--:|--:|
| 0.25 | 8,192 | 58,096 / 262,144 | 4.227 MB | 2.228 MB | 4.091 MB / 17.92 MB |
| 1 | 32,768 | 165,640 / 262,144 | 4.227 MB | 8.913 MB | 15.29 MB / 63.14 MB |
| 4 | 131,072 | 257,308 / 262,144 | 4.227 MB | 35.65 MB | 54.61 MB / 200.71 MB |
| 16 | 524,288 | 262,144 / 262,144 | 4.227 MB | 142.61 MB | 206.08 MB / 704.64 MB |

The semantic model counts bitmap bits and 16-byte values.  The cache-line
model charges 64-byte random bitmap and value accesses plus value write
allocation/writeback.  At load one, direct addressing moves fewer semantic
bytes but more cache-line bytes than the sequential radix reference.

## Isolated timing

The release binary ran with one thread pinned to CPU 4 on an AMD EPYC 9V74.
The isolation wrapper moved 368 threads, observed 0.06 seconds of other CPU
work over 71.56 seconds, and marked the run uncontended.  Peak RSS was
152,040 KiB.  The native-addition A/A ratio was 1.00436.

The primary cell used seven balanced rotating-order repetitions at universe
`2^20` and `2^17` starts:

| variant | request-generation passes | median | minimum | ns/start | additions/start | routing-budget multiple |
|:--|--:|--:|--:|--:|--:|--:|
| complete materialised 3x12 | 1 | 401.610 ms | 377.955 ms | 3,064.041 | 3.557042 | 10.6368x |
| complete direct-address replay | 2 | 1,002.350 ms | 957.313 ms | 7,647.323 | 8.877771 | **26.5475x** |

All repetition checksums were identical.  The direct/reference median ratio
was 2.49583.

The five-repetition load-one timing ladder was:

| universe / starts | materialised median | direct median | direct / materialised |
|:--|--:|--:|--:|
| `2^16 / 2^13` | 23.209 ms | 46.011 ms | 1.982 |
| `2^18 / 2^15` | 94.905 ms | 213.069 ms | 2.245 |
| `2^20 / 2^17` | 408.959 ms | 965.709 ms | 2.361 |
| `2^22 / 2^19` | 1.588 s | 4.623 s | 2.912 |

The fitted time exponents are 1.01979 for materialised radix and 1.10663 for
direct replay.  Direct algorithmic memory scales with exponent 0.999997 in
universe size and with exponent zero in start count.

## Optimistic full-universe projection

These rows add the measured scaled-table replay equivalents to round 40's
expected distinct-pair construction count.  They are deliberately marked
inapplicable because the 557-GB random-access tier was not measured.

| batch | pairs/start | combined additions/start | ratio to rho | fixed state | semantic traffic | cache-line traffic | tier measured? |
|:--|--:|--:|--:|--:|--:|--:|:--:|
| `2^38` | 0.125736 | 9.003507 | 1.037095669 | 557.31 GB | 106.39 TB | 356.27 TB | no |
| `2^40` | 0.031434 | 8.909205 | 1.036692143 | 557.31 GB | 423.87 TB | 1.412 PB | no |
| saturated limit | approximately 0 | 8.877771 | **1.036557634** | 557.31 GB | linear in starts | linear in starts | no |

The `2^50` storage cap does not bind this representation.  Traffic and
request regeneration remain linear in starts, so increasing the batch cannot
amortise the measured 8.878-addition replay term.

## Promotion gates

| requirement | status | evidence or gap |
|:--|:--|:--|
| zero false positives and false negatives | passed | all toy, load-one, and saturation cells |
| exact equality with materialised reference | passed | all accumulators, distinct counts, and digests |
| structured residual degree at most 5 | failed / unset | measured family unknown; comparison base's degree 4 is not transferable |
| relation collection below `2^120` | failed / unset | no valid 138,031-row Dickson composition |
| cost per usable relation below `2^103` | failed / unset | same composition gap |
| materialised storage below `2^50` | passed in formula | fixed state is `2^39.020` bytes |
| measured complete selector below rho in applicable tier | failed | 26.55x routing budget at scaled tier; full tier unmeasured |
| no discarded branch counted as exhaustive | passed | every request and start is replayed |

## Decision

Do not promote and do not attempt a full-depth unplanted relation.  Fixed
direct addressing solves round 41's storage obstruction but replaces it with
two complete request-generation passes and increasingly unfavorable random
access.  The best optimistic saturated projection remains 1.036558 times rho.

The remaining useful direction is not another table layout.  It must remove
the second generation pass or make pair reuse local in the generator's native
order, while preserving start entropy and exactness.  Otherwise the linear
replay term dominates even after pair construction is fully amortised.

## Reproduction

```bash
cargo test --bin p256_direct_address
cargo clippy --bin p256_direct_address -- -D warnings -A clippy::mismatched-bit-width-type
cargo build --release --bin p256_direct_address
RAYON_NUM_THREADS=1 python3 tools/isolated_bench.py run --wait --cpus 4 \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_direct_address_round42_20261006/isolation.jsonl \
  --label round42-direct-address -- \
  target/release/p256_direct_address \
  --round41 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_two_pass_router_round41_20261006/two-pass-router-result.json \
  --isolation research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_two_pass_router_round41_20261006/isolation.jsonl \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_direct_address_round42_20261006/direct-address-result.json
```

Final result artifact: 24,498 bytes, SHA-256
`82adf1c47c60a9004a83ac061d9a75ce908f3bb9c28f537260b43b7d5016aa2a`.
Semantic evidence SHA-256:
`79d7fdffb583b3a5f378821399a5fb076e79d25393848e61c6e6bfd1177c824a`.
Isolation JSONL SHA-256:
`7a6d11512161cee3020f599e8e705f6ea8fa68ba4e14aada73e0716230d4dad8`.
