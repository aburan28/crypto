# P-256 serpentine sign-orbit streaming, round 291: result

Date run: 2026-10-06

Serpentine sign ordering removes the pair table and the radix router while
preserving every sign and path state of each selected 17-column tuple.  The
complete `2^17` native P-256 orbit emitted 29,884,416 states with zero
transition, endpoint, false-positive or false-negative failures.

The conservative counted cost is 1.004439 P-256 addition equivalents per
state, including one target correction per sign segment, initialization, and
charging every endpoint doubling as two additions.  Applied to the frozen
round-33 local-oracle boundary, this is `0.9686171464` rho.  Peak live state is
9,808 bytes instead of a 1.141-TB pair table or the Round 290 radix arrays.

This is a verified construction-stage crossing, not end-to-end parity.
Achieved P-256 coverage of the correlated stream, the same-family structured
degree, the 138,031-row collection, sparse linear algebra, and non-genericity
remain unestablished.  No full-depth unplanted relation was attempted.

## Requested scope and actual execution

| requirement | execution |
|:--|:--|
| retain all 17 variable columns | complete 17-column tuples; no column fixed by the algorithm |
| retain every sign branch | complete `2^17` Gray orbit executed for the primary P-256 tuple |
| remove materialized routing | no pair requests, pair table, radix records, or routing pass |
| exact controls | exhaustive toy multiplicity equality and native P-256 endpoint/transition replay |
| progressively increasing depths | primary depths 8, 10, 12, 14, 16 and 17; holdouts at 8, 10 and 12 |
| preserve promotion gates | incomplete coverage, degree and collection fields remain false rather than inferred |

The experiment imports Round 33, its transitively bound Round 31 coefficient
cycle, and the Round 43–290 manifest by exact SHA-256.  It verifies the curve,
family, cutoff 219, arity, mean capacity, coefficient digest, Round 108/290
identities, and both predecessors' non-promotion state.

## Construction

For each accepted unsigned tuple, the algorithm enumerates all sign masks in
binary-reflected Gray order.  Even sign ordinals traverse the registered
common-edge path forward; odd ordinals traverse the same path backward.  The
endpoint of one segment is consequently the start endpoint of the next.

A Gray transition flips one signed column.  Adding or subtracting the doubled
current endpoint point updates the entire 17-term sum with one group addition.
Tuple initialization and 17 endpoint doublings are amortized over 131,072 sign
segments.  Nothing is looked up in a pair table.

For capacity `D` and `L=2^17`, the charged operations are:

```text
path additions             L*D
sign-boundary additions    L-1
right-colour corrections   L
initial 17-term sum        16 additions
endpoint doubles           17, charged as 34 additions
states emitted             L*(D+1)
```

At the exact retained-tuple mean `D=224.36124930751055`, this gives
1.004438979 additions/state and `0.9686171464` rho.

## One boundary table, one unit

The unit is a native P-256 addition equivalent.  Stage rows do not include the
missing coverage, relation, or linear-algebra phases.

| variant | additions/state or start | ratio to rho | storage / traffic scope | correctness and classification |
|:--|--:|--:|:--|:--|
| Pollard rho | complete reference | 1.000000 | matched reference boundary | reference |
| Round 33 pair-table stage, with target correction | 1.035499/state | 0.998569 | 1.141-TB table; traffic unpriced | projected stage relaxation |
| Round 108 width-9 routing | 0.449605/start extra | 1.000493 lower bound | 17.83-MB measured cell; generation excluded | exact stage regression |
| direct table-free sign reconstruction | 1.071000/state | 1.032802 | small live state | exact baseline |
| **Round 291 serpentine counted** | **1.004439/state** | **0.968617** | **9,808 bytes live state** | **exact stage advance** |

The ratio moved because the algorithm removed work while preserving the
registered per-tuple state multiset.  It has not yet shown that those
correlated states achieve the occurrence support credited by Round 33.

## Exhaustive toy equality

Natural sign enumeration and serpentine enumeration produced identical state
multiplicities, not merely identical distinct sets:

| group / shape | tuples | sign segments | states/reference | distinct | multiplicity mismatches | FP / FN |
|:--|--:|--:|--:|--:|--:|--:|
| `p=257`, `B=10`, `R=3`, arity 3 | 120 | 960 | 4,416 | 73 | 0 | 0 / 0 |
| `p=65537`, `B=20`, `R=3`, arity 4 | 4,845 | 77,520 | 961,248 | 177 | 0 | 0 / 0 |

Direct and serpentine multiset digests match in each cell.  Transition and
endpoint failures are zero.  The small exact support—especially 177 states in
the larger toy group—is also a warning: ordering preserves the Round 33 model,
but the model's visit-count coverage upper bound is not an achieved-coverage
measurement.

## Native P-256 replay

The primary tuple is deterministic hash attempt 240 and has capacity 227.
Every cell compares incremental projective states and scalar coefficients with
direct 17-term reconstructions at both path endpoints.  Every common-edge
transition is executed in coefficient arithmetic and with the native P-256
complete group law.  One-unit mutations are rejected.

| class | sign bits | segments | emitted states | additions/state | stage / rho | direct endpoint checks | failures |
|:--|--:|--:|--:|--:|--:|--:|--:|
| primary | 8 | 256 | 58,368 | 1.005225 | 0.969376 | 512 | 0 |
| primary | 10 | 1,024 | 233,472 | 1.004596 | 0.968768 | 2,048 | 0 |
| primary | 12 | 4,096 | 933,888 | 1.004438 | 0.968617 | 8,192 | 0 |
| primary | 14 | 16,384 | 3,735,552 | 1.004399 | 0.968579 | 32,768 | 0 |
| primary | 16 | 65,536 | 14,942,208 | 1.004389 | 0.968569 | 131,072 | 0 |
| **primary** | **17** | **131,072** | **29,884,416** | **1.004388** | **0.968568** | **262,144** | **0** |
| holdouts | 8 / 10 / 12 | 5,376 total | 1,300,480 | 1.004081–1.005135 | 0.968271–0.969289 | 10,752 | 0 |

All nine native cells have zero transition failures, endpoint failures,
false positives and false negatives.  The full primary depth completed in
31.178 seconds including direct validation.  This validation time is not an
attack-cost measurement.

## Isolated kernel timing and fits

The successful release run used one thread pinned to CPU 4 on an AMD EPYC
9V74.  It was uncontended, with 0.04 seconds of other CPU work during 65.954
seconds and peak RSS 20,608 KiB.  The isolation log also preserves the first
0.005-second dependency-parser failure rather than overwriting it.

At sign depth 12, seven alternating-order repetitions compared table-free
direct reconstruction with serpentine generation.  Both variants include one
target-correction addition per segment and exclude validation:

| variant | median | minimum | relative |
|:--|--:|--:|--:|
| direct reconstruction | 859.939 ms | 847.865 ms | 1.000 |
| serpentine | 819.475 ms | 793.673 ms | 1.04938x faster by medians |

Per-variant checksums were identical across repetitions.  The two-run native
addition A/A ratio was 1.02291, so the wall result is descriptive and has no
registered 95% paired interval.  Counted operations remain the primary
evidence.

Across primary sign depths 8–17, the counted-operation exponent is 0.999890,
emitted-state exponent 1.000000, validation-included time exponent 1.000836,
and live-memory exponent zero.  These are within-tuple sign-orbit fits; none is
transferred to global tuple coverage or a complete DLP.

## Promotion gates

| gate | status | evidence or gap |
|:--|:--|:--|
| zero FP/FN | passed | both exhaustive toys and all nine native cells |
| exact direct and native replay | passed | complete toy multiplicities, 29.88 million full-depth primary native states, endpoint reconstructions |
| counted construction stage below rho | passed | `0.9686171464` rho at the exact mean capacity |
| measured complete selector below rho | failed / unset | timing covers generation only; collision detection and other phases absent |
| proved P-256 coverage or usable-relation probability | failed / unset | Round 33 visit count remains an upper bound |
| same-family structured degree at most 5 | failed / unset | comparison-base degree 4 is not transferable |
| collection below `2^120` | failed / unset | no valid 138,031-row composition |
| per usable row below `2^103` | failed / unset | same missing composition |
| peak materialized storage below `2^50` | passed for this stage | 9,808 live bytes |
| demonstrably non-generic end-to-end | failed / unset | no collision/coverage proof or complete recovery |

## Decision and next falsification target

Retain Round 291 as an exact stage candidate.  Do not claim rho parity and do
not attempt a full-depth relation.

The radix-routing obstruction is removed for this ordering.  The immediate
test is now exact within-orbit support and collision multiplicity at the full
P-256 `2^17` sign depth, followed by overlap across independently selected
tuples.  If the achieved distinct-state rate or cross-tuple overlap consumes
the 3.14% arithmetic margin, this route fails before degree or collection.
If it survives, same-family degree and end-to-end collection remain mandatory.

## Reproduction

```bash
cargo test --release --bin p256_serpentine_sign_orbit
cargo clippy --release --bin p256_serpentine_sign_orbit -- -D warnings
cargo build --release --bin p256_serpentine_sign_orbit
RAYON_NUM_THREADS=1 python3 tools/isolated_bench.py run --wait --cpus 4 \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_serpentine_sign_orbit_round291_20261006/isolation.jsonl \
  --label round291-serpentine-sign-orbit -- \
  target/release/p256_serpentine_sign_orbit \
  --round33 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_biased_restart_round33_20261006/biased-restart-result.json \
  --round31 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_low_delta_round31_20261006/low-delta-result.json \
  --round290 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_radix_portfolio_round43_290_20261006/manifest.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_serpentine_sign_orbit_round291_20261006/serpentine-result.json
```

Result artifact: 18,095 bytes, SHA-256
`068be64f508ea2fb2996c265426dd8213088d1cad905859470792fcbc80972da`.
Semantic evidence SHA-256:
`453d14bc5a9e68ab5c6ba9d73ffa9040bc3c8557be7d73b0eb42d5395d939710`.
Two-record isolation JSONL: 4,319 bytes, SHA-256
`922d71e39e960ed60d6514882037c71e7fc2c15a1c79599a717ef6f7c6b2d0a0`.
