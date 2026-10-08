# P-256 exact radix portfolio, rounds 43–290: result

Date run: 2026-10-06

No round is promoted.  The finite sweep closes radix widths 6 through 36 for
the existing round-33 known-log selector representation: the fastest exact
routing kernel already consumes 1.34447 times the entire routing allowance,
before request generation, pair construction, relation collection, or sparse
linear algebra.

Round 108, width 9 at `2^16` starts, is the robust measured winner.  Its
seven-run median is 399.337 ns/start, or 0.449605 native P-256 addition
equivalents per start.  The registered routing allowance is only 0.334411,
and the resulting routing-stage lower-bound ratio is 1.000493 to rho.  This
cannot establish a complete attack cost below rho because request generation
was deliberately excluded from the matched width comparison.

Round 290, width 36 at `2^40` starts, minimizes formula traffic by using one
radix pass.  It is not an executed result.  Its counter table alone is
549,755,813,888 bytes, projected peak materialization is
299,616,918,568,960 bytes (`2^48.090`), and projected logical traffic is
846,074,197,573,632 bytes.  Those exact formulas pass the `2^50` storage cap,
but the memory tier was not allocated or timed.

## Scope and preregistration

Rounds 43–290 are one preregistered 248-cell Cartesian portfolio, not 248
adaptive discoveries.  The frozen map combines all radix widths `6..36`
with start depths `12, 14, 16, 18, 20, 24, 32, 40`:

    round = 43 + 31 * depth_index + (radix_width - 6)

Every round has a separate deterministic JSON artifact and a path/size/hash
entry in the manifest.  An independent post-run pass verified all 248 file
hashes and confirmed 248 unique round numbers covering exactly 43 through
290.

The measured family is the round-33 known-log low-delta/scalar selector.
`FB1h2f8621cda105`, with 131,458 columns and an observed structured-residual
maximum of 4, remains only a comparison factor base.  Its degree does not
transfer to this routing family.  The curve name and comparison-base identity
are preserved in every round artifact using the repository's ICV1 naming and
factor-base conventions.

## Algorithm

For each width `w`, the router applies `ceil(36 / w)` stable
least-significant-digit radix passes to the complete frozen stream of eight
requests per start.  The last pass allocates counters for only the remaining
key bits.  Executed cells compare:

- the entire sorted `(pair, owner_slot)` sequence with an independently
  canonicalized reference;
- every reconstructed accumulator with an independently regenerated value;
- distinct-pair counts and deterministic BLAKE3 digests; and
- one-unit negative controls for every start.

No request, start, or probabilistic branch is discarded.  The candidate
traffic formula charges generation, every radix read/write, the sorted scan,
accumulator movement, and histogram scans.  Routing performs no native field
or group operations; native additions are used only as the timing-equivalent
control.

## Complete and formula-only cells

All 76 registered complete cells executed: widths 6 through 24 at each of
`2^12`, `2^14`, `2^16`, and `2^18` starts.  Every complete cell had zero
sorted mismatches, zero accumulator mismatches, zero regeneration failures,
zero false positives, and zero false negatives.

The other 172 cells are projection-only.  This includes widths 25 through 36
at the four smaller depths and every width at depths `2^20` through `2^40`.
Their wall-time, exactness, false-positive, false-negative, and applicable-tier
fields are null.  None inherits an execution result from another cell.

At Round 108 the exact counts are:

| quantity | value |
|:--|--:|
| starts / requests | 65,536 / 524,288 |
| distinct pairs | 524,283 |
| radix passes / counters per pass | 4 / 512 |
| histogram increments / scatter writes | 2,097,152 / 2,097,152 |
| sorted scans / accumulator updates | 524,288 / 524,288 |
| logical traffic | 100,712,448 bytes |
| peak algorithmic state | 17,829,888 bytes |
| negative controls / FP / FN | 65,536 / 0 / 0 |

## Isolated timing

The release binary ran with one Rayon thread pinned to CPU 4 on an AMD EPYC
9V74.  The isolation wrapper moved 371 threads, observed zero seconds of other
CPU work during the 14.822-second run, and marked it uncontended.  Peak RSS
was 290,644 KiB.  The native-addition A/A ratio was 1.04785 and the
seven-repetition median control cost was 888.194 ns/addition.

The five widths selected by their one-shot `2^16` timing were measured for
seven balanced rotating-order repetitions on the identical pre-generated
request vector:

| round | width | median ns/start | additions/start | routing-budget multiple | routing lower bound / rho |
|--:|--:|--:|--:|--:|--:|
| 107 | 8 | 424.806 | 0.478280 | 1.43022x | 1.000616 |
| **108** | **9** | **399.337** | **0.449605** | **1.34447x** | **1.000493** |
| 112 | 13 | 424.942 | 0.478434 | 1.43068x | 1.000616 |
| 111 | 12 | 417.056 | 0.469556 | 1.40413x | 1.000578 |
| 114 | 15 | 437.599 | 0.492684 | 1.47329x | 1.000677 |

All repetition checksums were identical.  These timings include clone,
allocation, radix passes, scan, and reconstruction, but exclude shared
deterministic request generation.  They therefore compare routing kernels;
they are not complete-selector or relation timings.

## Scaling fits

Observed time exponents were fitted independently for each executed width
over the four complete depths.  Selected fits are:

| width | time exponent | memory exponent |
|--:|--:|--:|
| 6 | 1.12473 | 0.99999 |
| 9 | 1.11572 | 0.99990 |
| 12 | 1.15940 | 0.99920 |
| 18 | 1.05988 | 0.96779 |
| 24 | 0.32406 | 0.79226 |
| 36 | not measured | 0.27790 |

The declining high-width exponents are a small-depth fixed-histogram effect,
not evidence of sublinear full-depth work.  No time fit for widths 25–36 was
measured, and none of the fits is extrapolated across an unmeasured memory
tier.  Exact memory formulas, rather than fitted timings, define the
projection-only cells.

## Round 290 formula boundary

| quantity | value |
|:--|--:|
| starts / requests | 1,099,511,627,776 / 8,796,093,022,208 |
| radix passes | 1 |
| 36-bit histogram counters | 68,719,476,736 |
| histogram counter bytes | 549,755,813,888 |
| record traffic | 844,424,930,131,968 bytes |
| histogram traffic | 1,649,267,441,664 bytes |
| total logical traffic | 846,074,197,573,632 bytes |
| peak algorithmic state | 299,616,918,568,960 bytes |
| executed / applicable tier measured | no / no |

The width-36 cell saves passes but cannot remove the two full materialized
record arrays or the linear scan and accumulation traffic.  More radix-width
screening within this representation cannot overcome the already-exceeded
routing allowance.

## End-to-end accounting and promotion gates

There is no valid projection for collecting 138,031 independent P-256 rows
from this experiment.  The measured family has no established structured
residual degree, relation probability, duplicate/rank allowance, or sparse
linear-algebra composition.  Filling those fields from
`FB1h2f8621cda105` would transfer evidence between different families.
Moreover, the best measured routing component alone exceeds its rho allowance
before any omitted complete-selector work.  Reporting a numeric collection
cost under those conditions would manufacture an extrapolation rather than
measure one.

| requirement | status | evidence or gap |
|:--|:--|:--|
| all 248 registered rounds present once | passed | exact range 43–290; all artifact hashes rechecked |
| all 76 registered complete cells executed | passed | widths 6–24 at four depths |
| zero FP and FN on complete cells | passed | all 76 cells |
| structured residual degree at most 5 | failed / unset | measured family unknown; comparison-base degree 4 is not transferable |
| relation collection below `2^120` | failed / unset | no valid 138,031-row composition |
| cost per usable relation below `2^103` | failed / unset | same composition gap |
| materialized storage below `2^50` | passed in formula | Round 290 peak is `2^48.090` bytes |
| complete selector below rho in applicable tier | failed | best routing-only kernel is already 1.34447x its allowance |
| no discarded branch counted as exhaustive | passed | complete cells route every request; projections claim no exactness |

## Decision

Do not promote any round and do not attempt a full-depth unplanted relation.
Round 108 is the implementation choice only if this router is needed for a
scaled diagnostic; it is not an attack candidate.  Round 290 is a resource
formula, not a benchmark.

The dominant obstruction is the linear materialized routing term, not radix
digit choice.  A credible next design must coalesce reuse in the generator's
native order or avoid materializing both record arrays, while retaining all
start entropy and exact replay.  Another width in the same stable-radix family
is now a closed direction.

## Reproduction

```bash
cargo test --bin p256_radix_portfolio
cargo clippy --bin p256_radix_portfolio -- -D warnings -A clippy::mismatched-bit-width-type
cargo build --release --bin p256_radix_portfolio
RAYON_NUM_THREADS=1 python3 tools/isolated_bench.py run --wait --cpus 4 \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_radix_portfolio_round43_290_20261006/isolation.jsonl \
  --label round43-290-radix-portfolio -- \
  target/release/p256_radix_portfolio \
  --round42 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_direct_address_round42_20261006/direct-address-result.json \
  --isolation research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_direct_address_round42_20261006/isolation.jsonl \
  --round-dir research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_radix_portfolio_round43_290_20261006/rounds \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_radix_portfolio_round43_290_20261006/manifest.json
```

Manifest: 109,072 bytes, SHA-256
`41cdec6ad5265464cc512b2f026f6fa77b333aa209f163e7e77ab78ba665b876`.
Semantic evidence SHA-256:
`fe852d9a50e29cc3331e901858ac4b2c1d2e5c00a710e0bd392087b67137bd2b`.
Isolation JSONL: 2,126 bytes, SHA-256
`ccf48d7d77d848da677cc7eb211cd4f738b61955e8b75d1921abf428fbec2335`.
The 248 per-round artifacts total 652,943 bytes; each file's byte count and
SHA-256 are bound by the manifest.
