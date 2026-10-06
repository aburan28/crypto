# P-256 two-pass independent router, round 41: result

Date run: 2026-10-06

The two-pass 18-bit radix router is exact and reduces registered logical
traffic by 20%, from 160 to 128 bytes per request.  It does **not** provide a
stable timing improvement over the three-pass 12-bit router, and it does not
establish parity with rho.

At the registered primary cell of `2^16` starts, the pre-generated 2x18
candidate cost 0.224747 native P-256 addition equivalents per start.  That
fits the isolated routing allowance of 0.334411, but its 13.208-ms median was
1.526 times the matched pre-generated 3x12 control's 8.653 ms.  In the
separate depth ladder the ordering changed: 2x18 was faster at `2^14` and
`2^18`, but slower at `2^16`.  Its fitted time exponent was 1.376, versus
1.337 for 3x12.  The result is therefore a traffic reduction, not a robust
RAM-time reduction.

An optimistic stage projection combines the pre-generated routing kernel with
round 40's pair-construction model.  It first falls below rho at the reported
`2^40` batch, where it would materialise 299.07 TB and move 1.126 PB of
logical traffic.  The largest registered batch below `2^50` bytes would
materialise 1.126 PB and move 4.239 PB.  Neither storage tier was measured,
and pre-generated routing excludes upstream request derivation.  These rows
are not end-to-end speedups.

## Scope

This round measures only the independent-pair routing kernel for the
round-33 known-log low-delta/scalar selector family.  It imports and verifies
round 40's request digest, native-control digest, exactness gates, pair-only
threshold, selected storage row, and isolated-run status.

`FB1h2f8621cda105` remains a separate comparison base.  Its maximum observed
structured residual degree of 4 is not transferred to the measured family.
No result here prices the collection of 138,031 Dickson-base rows, sparse
linear algebra, or a complete P-256 discrete logarithm.

## Exactness and traffic

Every start contributes all eight requests.  Both radix implementations sort
the entire 36-bit pair identifier, scan every distinct pair, and reconstruct
one accumulator for every start.  Independent regeneration and one-unit
negative controls were applied at every depth.

| starts | requests | distinct pairs | sorted / accumulator mismatches | route failures | FP / FN | 3x12 traffic | 2x18 traffic | 3x12 / 2x18 peak RAM |
|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| `2^12` | 32,768 | 32,768 | 0 / 0 | 0 / 0 | 0 / 0 | 5.24 MB | 4.19 MB | 1.11 / 3.21 MB |
| `2^14` | 131,072 | 131,072 | 0 / 0 | 0 / 0 | 0 / 0 | 20.97 MB | 16.78 MB | 4.46 / 6.55 MB |
| `2^16` | 524,288 | 524,283 | 0 / 0 | 0 / 0 | 0 / 0 | 83.89 MB | 67.11 MB | 17.83 / 19.92 MB |
| `2^18` | 2,097,152 | 2,097,078 | 0 / 0 | 0 / 0 | 0 / 0 | 335.54 MB | 268.44 MB | 71.30 / 73.40 MB |

The depth-`2^16` sorted digest is
`b4178238067a37febeeb46116dffc633311549e6325461670bd12db7383e0dd6`,
identical to round 40's frozen request digest.  Nine complete toy cells over
pair universes 17, 31, and 63 also matched direct references exactly.

The fixed logical-traffic models are:

```text
3x12: initial clone/write 16 + three read/write passes 96
      + sorted scan 16 + accumulator read/write 32 = 160 bytes/request
2x18: initial clone/write 16 + two read/write passes 64
      + sorted scan 16 + accumulator read/write 32 = 128 bytes/request
```

The 2x18 count table adds 2,097,152 bytes of peak RAM.  It therefore has
higher peak RAM at every checked depth despite moving fewer bytes.

Round 40's 4,096 native P-256 group replays remain the group-level control.
Their digest,
`98c8f56a06257256bb88f6c6d7cffce79391c461170b91287f8fd193fe5711e6`,
was verified before this routing-only experiment ran.

## Isolated timing

The release binary ran with one thread pinned to CPU 4 on an AMD EPYC 9V74.
The isolation wrapper moved 362 threads, observed zero seconds of other CPU
work, marked the run uncontended, and recorded 193,132 KiB peak RSS.  Total
wall time was 9.128 seconds.  The native-addition A/A ratio was 1.01633.

The primary cell used seven balanced rotating-order repetitions at `2^16`
starts:

| variant | request generation | median | minimum | ns/start | additions/start | routing budget multiple |
|:--|:--:|--:|--:|--:|--:|--:|
| complete 3x12 | yes | 226.665 ms | 207.046 ms | 3,458.631 | 3.857005 | 11.5337x |
| pre-generated 3x12 control | no | 8.653 ms | 6.302 ms | 132.035 | 0.147243 | 0.4403x |
| pre-generated 2x18 candidate | no | 13.208 ms | 8.969 ms | 201.534 | 0.224747 | 0.6721x |

All repetition checksums were identical.  The complete row includes the
SHA-256 rejection sampler used to create the deterministic experimental
stream.  The pre-generated rows clone the frozen input inside the timed
interval but do not derive it; they may price only the routing stage.

The registered five-repetition timing ladder was:

| starts | 3x12 median | 2x18 median | 2x18 / 3x12 |
|--:|--:|--:|--:|
| `2^14` | 2.265 ms | 1.835 ms | 0.8102 |
| `2^16` | 9.804 ms | 10.618 ms | 1.0830 |
| `2^18` | 92.132 ms | 83.310 ms | 0.9042 |

The log-log fitted exponents over these three depths are 1.33659 for 3x12 and
1.37619 for 2x18.  The changing ordering and superlinear fits show a
cache/allocation transition within the measured range.  They do not support
extrapolating either RAM timing to multi-terabyte batches.

## Optimistic stage projection

The following rows combine the primary 2x18 kernel's 0.224747
addition-equivalents per start with round 40's expected distinct-pair
construction counts.  They exclude request derivation and do not import a
measured external-memory cost.

| batch | pair constructions/start | combined additions/start | ratio to rho | peak materialisation | traffic/batch | measured tier? |
|:--|--:|--:|--:|--:|--:|:--:|
| `2^38` | 0.125736 | 0.350483 | 1.000068776 | 74.77 TB | 281.47 TB | no |
| `2^40` | 0.031434 | 0.256181 | **0.999665250** | 299.07 TB | 1.126 PB | no |
| largest below `2^50` | 0.008350 | 0.233096 | **0.999566470** | 1.126 PB | 4.239 PB | no |

The sub-rho values are arithmetic-only stage projections.  The measured
complete experimental 3x12 route remains above rho, and no complete 2x18
selector or applicable external-memory run exists.

## Promotion gates

| requirement | status | evidence or gap |
|:--|:--|:--|
| zero false positives and false negatives | passed | all toy and depth-ladder controls |
| exact equality with 3x12 reference | passed | sorted records and all reconstructed accumulators |
| structured residual degree at most 5 | failed / unset | measured family unknown; comparison base's degree 4 is not transferable |
| relation collection below `2^120` | failed / unset | no valid 138,031-row Dickson composition |
| cost per usable relation below `2^103` | failed / unset | same composition gap |
| materialised storage below `2^50` | passed only as a projection | largest row is just below the registered cap |
| measured complete selector below rho in applicable tier | failed | primary complete row is 11.53x routing budget; projected tier unmeasured |
| no discarded branch counted as exhaustive | passed | every request and start retained |

## Decision

Do not promote and do not attempt a full-depth unplanted relation.  The 2x18
router is an exact 20% traffic reduction, but it is not a stable RAM-time
improvement and its arithmetic-only crossover requires an unmeasured
hundreds-of-terabytes tier.  Local routing headroom cannot be presented as
end-to-end parity.

The next bounded screen should target representation, not a wider in-memory
radix: either eliminate materialisation with an exact generator-compatible
partition, or prove a lower bound showing that exact global deduplication
necessarily moves linear traffic.  Any candidate must be measured as a
complete selector and must establish its own structured degree before it can
be composed with relation collection.

## Reproduction

```bash
cargo test --bin p256_two_pass_router
cargo clippy --bin p256_two_pass_router -- -D warnings -A clippy::mismatched-bit-width-type
cargo build --release --bin p256_two_pass_router
RAYON_NUM_THREADS=1 python3 tools/isolated_bench.py run --wait --cpus 4 \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_two_pass_router_round41_20261006/isolation.jsonl \
  --label round41-two-pass-router -- \
  target/release/p256_two_pass_router \
  --round40 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_independent_batch_round40_20261006/independent-batch-result.json \
  --isolation research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_independent_batch_round40_20261006/isolation.jsonl \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_two_pass_router_round41_20261006/two-pass-router-result.json
```

Final result artifact: 16,786 bytes, SHA-256
`8bf3dd1d2e7942858c6cbf8486170cb6275aa0ec5e5380e686daeebcd9fd26a6`.
Semantic evidence SHA-256:
`d6544e4ca6470482003b8d0badd40ed2bf6dc151987bbda41d2329fb9104ab44`.
Isolation JSONL SHA-256:
`6e6dedb3f9367229808b1d3eb7c9978c9445fd01643a72078a8489215d9ebf32`.
