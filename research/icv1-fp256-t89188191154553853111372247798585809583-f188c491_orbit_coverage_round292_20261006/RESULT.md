# P-256 serpentine-orbit coverage, round 292: result

## Verdict

The exact coverage test falsifies the Round 291 independence assumption.
The full primary orbit emits `29,884,416` states but covers only `5,845,754`
distinct P-256 generator multiples, a distinct fraction of `0.1956121210`.
Charging the same construction work per distinct state changes the Round 291
stage ratio from `0.968617` to `4.951724` rho (`4.951470` using the primary
tuple's actual capacity).

Combining four deterministic frozen tuples makes the fraction worse, not
better: `121,634,816` emitted states cover `10,380,468` distinct elements, or
`0.0853412562`.  The corrected four-tuple stage ratio is `11.349928` rho
(`11.348493` with actual capacities).  This fails before summation-polynomial
degree, relation probability, row collection, sparse linear algebra, or scalar
recovery are charged.  No full-depth unplanted relation was attempted.

## Requirements and evidence

| requirement | execution | result |
|:--|:--|:--|
| frozen ICV1 identity and `FB1h2f8621cda105` | Round 291 artifact and transitive Round 31 coefficient source are hash-bound | passed |
| all 17 columns variable | imports the four frozen Round 291 17-column tuples without fixing a column | passed |
| exact support, not a hash estimate | inclusive modular intervals use exact `BigUint` endpoints and exact sort/merge | passed |
| exhaustive references where feasible | two toy groups plus native depths 8, 10, and 12 | zero FP/FN |
| Round 291 stream identity | all six primary boundary digests, including depth 17, are replayed | matched |
| progressively increasing depths | primary depths 8, 10, 12, 14, 16, 17 | completed |
| cross-tuple coverage | cumulative full-depth portfolios of 1, 2, 3, and 4 frozen tuples | completed exactly |
| group replay of reported relations | no relations are reported; coefficients are exact modulo the generator order | not applicable; zero claims |
| resource and operation evidence | prefix steps, signed sums, sort/merge comparisons, logical bytes, RSS, isolated CPU/wall | recorded |
| structured degree and end-to-end projection | not produced by this coverage experiment | remain failed/unset |

## Why the interval result is exact

Every common-edge path transition changes the scalar coefficient by exactly
`+1` or `-1`.  A nearest-neighbour integer walk visits every integer between
its minimum and maximum prefix sum.  Thus one sign segment's support is exactly
one inclusive interval translated by its signed starting coefficient modulo
the P-256 subgroup order; the implementation splits the interval only if it
wraps through zero.  Exact sorting and union of all intervals gives the exact
set of covered generator multiples.  No fingerprint decides equality.

The primary path has 228 emitted positions per sign segment, but the maximum
interval width is 228 and the average full-depth width is only 83.05.  Path
backtracking creates the first loss.  Overlap among translated intervals
creates the second.

## One boundary, one unit

The unit is counted P-256 addition equivalents per distinct covered group
element, normalized to the frozen Round 291 selector-stage rho boundary.
These are construction-stage ratios, not complete DLP ratios.

| variant | emitted | exact distinct | distinct / emitted | corrected ratio to rho | correctness |
|:--|--:|--:|--:|--:|:--|
| Pollard rho | reference | reference | 1.000000 | **1.000000** | complete-method boundary |
| Round 291, uncorrected emitted-state row | 29,884,416 | unset | assumed 1.000000 | 0.968617 | superseded by exact coverage |
| primary, depth 8 | 58,368 | 36,856 | 0.631442 | 1.533977 | direct exact |
| primary, depth 10 | 233,472 | 118,606 | 0.508010 | 1.906691 | direct exact |
| primary, depth 12 | 933,888 | 411,359 | 0.440480 | 2.199004 | direct exact |
| primary, depth 14 | 3,735,552 | 1,290,842 | 0.345556 | 2.803069 | interval exact; boundary replay |
| primary, depth 16 | 14,942,208 | 3,712,464 | 0.248455 | 3.898564 | interval exact; boundary replay |
| **primary, depth 17** | **29,884,416** | **5,845,754** | **0.195612** | **4.951724** | **interval exact; boundary replay** |

The depth fit is `1.000000` for emitted states but only `0.819206` for exact
distinct states.  It describes this frozen orbit only and is not extrapolated
to global relation probability.

## Duplicate decomposition

| full-depth portfolio | emitted | per-segment support sum | exact union | within-segment duplicates | cross-segment duplicates | corrected ratio |
|--:|--:|--:|--:|--:|--:|--:|
| 1 tuple | 29,884,416 | 10,885,760 | 5,845,754 | 18,998,656 | 5,040,006 | 4.951724 |
| 2 tuples | 60,293,120 | 21,673,300 | 8,147,348 | 38,619,820 | 13,525,952 | 7.168093 |
| 3 tuples | 89,128,960 | 32,214,740 | 9,338,950 | 56,914,220 | 22,875,790 | 9.244277 |
| **4 tuples** | **121,634,816** | **43,607,440** | **10,380,468** | **78,027,376** | **33,226,972** | **11.349928** |

The second, third, and fourth tuples add only `2,301,594`, `1,191,602`, and
`1,041,518` new elements respectively, despite emitting roughly 29--32 million
states each.  The shared two-delta coefficient structure aligns the orbit
supports strongly; these tuples do not behave as independent random samples.

## Correctness controls

The `p = 257` toy emits 80 states and both direct enumeration and interval
union find 46 distinct states.  The `p = 65,537` toy emits 672 states and both
find 314.  Duplicate counts match exactly; both cells have zero false
positives and zero false negatives.

Native direct scalar sets at depths 8, 10, and 12 contain 36,856, 118,606, and
411,359 elements.  The interval unions are identical with zero false
positives and zero false negatives.  All six primary start-coefficient digests
match the independently frozen Round 291 digests.  Because the repository's
P-256 generator has order `n`, equality of coefficients modulo `n` is equality
of the represented group elements.

## Operations and resources

The isolated native run completed in 1.821 seconds (1.737 user, 0.084 system)
on reserved CPU 4 of an AMD EPYC 9V74.  It was uncontended and recorded zero
other-process CPU seconds.  Peak RSS was 171,884 KiB.

The four-tuple cell performs 121,110,528 path-prefix steps, 8,912,896 signed
sum terms, 10,266,207 sort comparisons, and 524,287 merge comparisons.  It
materializes 524,288 source intervals, 71,390 merged intervals, and 33,554,432
logical endpoint bytes; no interval wrapped modulo `n`.  The primary
full-depth cell uses 8,388,608 logical endpoint bytes.  Both actual and
projected materialized storage are below `2^50` bytes.

The interval audit itself is additional selector work.  The corrected ratios
above conservatively retain the Round 291 construction charge and therefore
already fail the boundary before charging sorting, merging, relation testing,
or any later phase.

## Gates

| gate | status | evidence or obstruction |
|:--|:--|:--|
| zero FP/FN on complete direct controls | passed | two toys and three native depths |
| exact toy support and duplicate counts | passed | both complete groups agree |
| native stream identity | passed | six of six primary boundary digests match Round 291 |
| primary distinct fraction at least 0.968617 | **failed** | 0.195612 |
| primary corrected stage ratio at most 1 | **failed** | 4.951724 |
| four-tuple corrected stage ratio at most 1 | **failed** | 11.349928 |
| materialized storage below `2^50` | passed | 32 MiB logical endpoints; 171,884 KiB peak RSS |
| same-family structured residual degree at most 5 | failed / unset | no same-family polynomial solve |
| relation collection below `2^120` | failed / unset | coverage gate already fails; no usable-row model |
| per usable relation below `2^103` | failed / unset | same obstruction |
| demonstrably non-generic complete DLP below rho | failed / unset | no relation or scalar recovery |
| promoted | **no** | coverage hypothesis falsified |

## End-to-end consequence and decision

Any projection that counted all Round 291 emissions as independent trials
understated construction work per distinct group element by `5.112157x` on the
primary orbit and `11.717662x` on the four-tuple portfolio.  This is not yet a
usable-relation probability, so those factors are not promoted into a complete
138,031-row collection estimate.  They are sufficient to reject this
serpentine two-delta orbit as a rho-parity route: its construction stage alone
is already 4.95--11.35 times the boundary under exact coverage accounting.

Retain Round 291 as a table-free enumeration identity, but withdraw its
sub-rho reading as soon as distinct coverage is required.  Do not attempt an
unplanted full-depth relation for this candidate.  Further factor-base
screening must break both sources of correlation: nearest-neighbour path
backtracking and the cross-tuple alignment induced by the global two-delta
coefficient lattice.  Local degree improvements cannot repair this measured
coverage loss.

## Reproduction and hashes

```bash
cargo test --release --bin p256_orbit_coverage
cargo clippy --release --bin p256_orbit_coverage -- -D warnings
cargo build --release --bin p256_orbit_coverage --bin isolated_bench
target/release/isolated_bench run --wait --cpus 4 \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_orbit_coverage_round292_20261006/isolation.jsonl \
  --label round292-orbit-coverage -- \
  target/release/p256_orbit_coverage \
  --round291 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_serpentine_sign_orbit_round291_20261006/serpentine-result.json \
  --round31 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_low_delta_round31_20261006/low-delta-result.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_orbit_coverage_round292_20261006/orbit-coverage-result.json
```

Result artifact: 20,655 bytes, SHA-256
`19d546fd2e53728f029d349080edbf7337ac2ada785446978b91a62eabbe3821`.
Semantic evidence SHA-256:
`301a4ebd03aa1e7ea7c02ceab328729889a120d8c27b0980e6172bc86863efb1`.
Isolation JSONL: 2,073 bytes, SHA-256
`335a55e9b2c1ae20e22721a8ff201971b930856fc402257ea823ee463336fcdf`.
