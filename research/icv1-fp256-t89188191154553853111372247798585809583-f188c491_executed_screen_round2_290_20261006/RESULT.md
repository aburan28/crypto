# P-256 factor-base screening, 289 executed rounds: result

## Verdict

All 289 requested screening rounds, numbered 2 through 290 after the existing
Round 1 baseline, actually executed.  There are 289 separate result artifacts,
289 complete full-`2^17` sign-orbit evaluations, and zero projection-only
entries.  Every factor base is primitive, nonidentity, unique up to sign, and
has all 131,458 edges replayed with zero failures.

No candidate passes the selector gate.  The best is Round 189,
`LD2R7021h58521a0ef3de`.  Its monotone sign-oriented paths emit 29,753,344
states and exactly cover 13,096,682 P-256 generator multiples, a distinct
fraction of `0.4401751279`.  Its ideal-local construction-stage ratio is
`2.200455` rho and its R-specific factor-base-corrected lower bound is
`2.320918` rho.  Both already exceed rho before sorting, polynomial solving,
relation verification, row collection, sparse linear algebra, or recovery.

No structured-degree claim transfers to these candidates, no relation is
reported, and no unplanted full-depth relation was attempted.

## What “289 rounds” means here

Round 1 remains the existing hash-bound baseline.  This experiment performs
one fresh, uniform, deterministic screen for every integer round 2–290:

```text
R(round) = 6647 + 2(round - 2)
```

Thus it constructs 289 distinct exact known-log factor bases with odd rare-edge
counts 6647, 6649, ..., 7223.  Each round writes its receipt immediately and is
marked `execution_status: complete` only after factor-base validation, tuple
selection, the complete 131,072-sign orbit, exact interval union, accounting,
and timing finish.

This does not relabel the historical rounds 43–290 radix matrix.  Its 172
formula-only large cells remain formula-only; in particular, the historical
Round 290 configuration still requires about 300 TB.  The new sweep is the
requested set of 289 actual tractable screening iterations.

## Registered construction

Every factor base uses 131,458 columns, common delta one, and the exact rare
delta satisfying

```text
(131458 - R) + R * delta_rare = 0 mod n.
```

All registered `R` values are coprime to 131,458.  The implementation still
constructs every coefficient, verifies the complete cycle, searches for the
first valid deterministic anchor, checks nonidentity and uniqueness up to
sign, replays every edge, and hashes the stored coefficient sequence.  IDs use
the established `LD2R{R}h{digest-prefix}` format.  All 289 anchors passed on
their first attempt.

For each factor base, the first deterministic 17-column tuple with total
common-edge capacity at least 219 is selected.  Observed capacities range from
219 to 244, with median 223.

The selector reverses the coordinate direction for negative signs.  Every
step within a segment therefore changes the signed group sum monotonically by
one generator.  This removes every within-segment revisit seen in Round 292.
Gray-ordered segments alternate forward and reverse, retaining one precomputed
boundary addition per sign transition.  Exact modular interval union then
measures the remaining cross-sign overlap.

## One boundary, one unit

The table uses charged P-256 addition equivalents per exact distinct covered
group element, normalized to rho.  These are selector-stage lower bounds, not
complete DLP ratios.

| candidate | rare edges | emitted states | exact distinct | distinct fraction | ideal-local / rho | factor-base-corrected / rho | result |
|:--|--:|--:|--:|--:|--:|--:|:--|
| Pollard rho | reference | reference | reference | 1.000000 | **1.000000** | **1.000000** | complete-method boundary |
| Round 292 ordinary serpentine, R=6935 | 6,935 | 29,884,416 | 5,845,754 | 0.195612 | 4.951470 | not computed there | rejected |
| executed Round 146, monotone, R=6935 | 6,935 | 29,884,416 | 12,473,518 | 0.417392 | 2.320521 | 2.445966 | rejected |
| **executed Round 189, best** | **7,021** | **29,753,344** | **13,096,682** | **0.440175** | **2.200455** | **2.320918** | **rejected** |
| median of 289 executed rounds | — | — | — | 0.414058 | — | 2.605617 | rejected |
| worst corrected round, Round 12 | 6,667 | 28,966,912 | 10,766,008 | 0.371666 | 2.606377 | 5.334667 | rejected; finite-support retries dominate |

Sign-oriented traversal more than doubles exact coverage for the frozen
R=6935 base and removes all within-segment duplication, but it cannot remove
cross-sign alignment.  Even the best candidate loses 16,656,662 of its
29,753,344 emissions to interval overlap.  Its distinct fraction misses the
registered `0.968617` gate by more than a factor of two.

Across all 289 rounds, distinct fractions range from `0.339440` to `0.440175`;
the first quartile, median, and third quartile are `0.406582`, `0.414058`, and
`0.421178`.  Factor-base-corrected stage ratios range from `2.320918` to
`5.334667`, with median `2.605617`.  None is below one.

## Exactness controls

Two toy groups are exhaustively enumerated.  For `p=257`, direct enumeration
and interval union both find 66 distinct states from 80 emissions.  For
`p=65,537`, both find 592 from 672.  Both have zero false positives and zero
false negatives.

Native R=6935 controls compare direct scalar sets with exact interval union at
sign depths 8, 10, and 12.  The respective distinct counts are 58,368,
233,472, and 909,392, with zero false positives and false negatives.  The
R=6935 coefficient digest and its accepted attempt-240, capacity-227 tuple
match the frozen Round 31 / Round 291 identity.

An independent post-run pass reopens every artifact and checks its byte count,
SHA-256, round number, and completion status.  The manifest proves that every
integer 2–290 occurs exactly once.

## Aggregate work and resources

The 289 rounds execute:

| quantity | total |
|:--|--:|
| factor-base coefficients constructed | 37,991,362 |
| factor-base edges exactly replayed | 37,991,362 |
| complete sign segments | 37,879,808 |
| emitted path states | 8,523,743,232 |
| sum of per-round exact distinct states | 3,515,111,820 |
| charged P-256 addition equivalents | 8,561,632,288 |
| exact sort comparisons | 665,574,788 |
| exact merge comparisons | 37,879,519 |
| reported relations | 0 |

The isolated run completed in 59.721 seconds wall time (58.820 user, 0.898
system) on reserved CPU 4 of an AMD EPYC 9V74.  It was uncontended; other
processes consumed 0.03 CPU seconds.  Peak RSS was 100,172 KiB.  The maximum
per-round logical materialization remained far below `2^50` bytes.  The 289
round artifacts total 1,102,484 bytes.

## Gates

| gate | status | evidence or obstruction |
|:--|:--|:--|
| exactly 289 rounds covering 2–290 once | passed | manifest range and independent replay |
| every round actually complete | passed | 289 complete, 0 projection-only |
| every factor base exact and primitive | passed | 37,991,362 coefficients and edges checked |
| zero FP/FN on complete controls | passed | two toys plus three native cells |
| manifest hashes and byte counts replay | passed | all 289 artifacts |
| materialized storage below `2^50` | passed | 100,172 KiB observed peak RSS |
| any distinct fraction at least 0.968617 | **failed** | maximum 0.440175 |
| any ideal-local selector ratio below 1 | **failed** | minimum 2.200455 |
| any factor-base-corrected ratio below 1 | **failed** | minimum 2.320918 |
| same-family structured residual degree at most 5 | failed / unset | no polynomial solve for these bases |
| cost per usable relation below `2^103` | failed / unset | selector gate already fails; no relation model |
| 138,031-row collection below `2^120` | failed / unset | no usable-row or rank evidence |
| demonstrably non-generic complete DLP below rho | failed / unset | no relation, linear solve, or recovery |
| promoted | **no** | zero selector-gate passes |

## End-to-end consequence and decision

The 289-run screen establishes a useful negative boundary.  Sign-aware
direction fixes path backtracking, but factor-base coefficient structure still
causes enough cross-sign overlap that the cheapest construction stage is
2.32 times rho before any summation-polynomial or linear-algebra work.  Degree
reduction cannot repair a stage that has already crossed the full-method
boundary.

Do not promote any candidate and do not attempt an unplanted full-depth P-256
relation from this family.  A further candidate needs genuinely decorrelated
sign-segment translations, not another small adjustment of rare-edge density.
The next credible screen should change the coefficient geometry itself while
preserving exact known logs and the global low-cost boundary transition.

## Reproduction and hashes

```bash
cargo test --release --bin p256_executed_screen
cargo clippy --release --bin p256_executed_screen -- -D warnings
cargo build --release --bin p256_executed_screen --bin isolated_bench
target/release/isolated_bench run --wait --cpus 4 \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_executed_screen_round2_290_20261006/isolation.jsonl \
  --label round2-290-executed-screen -- \
  target/release/p256_executed_screen \
  --round1-factor-base research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_affine_bitbox_degree_round1_20261004/factor-base-result.json \
  --round1-degree research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_affine_bitbox_degree_round1_20261004/degree-result.json \
  --round292 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_orbit_coverage_round292_20261006/orbit-coverage-result.json \
  --round-dir research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_executed_screen_round2_290_20261006/rounds \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_executed_screen_round2_290_20261006/manifest.json
```

Manifest: 193,193 bytes, SHA-256
`8723757f965eec5cfb1310de1091c8668815212c31823fd30b976de48ebadf45`.
Semantic evidence SHA-256:
`3576baaeb3fb89b7b1188726da1266a5afb01e03efe82509f212b186acf6d46e`.
Isolation JSONL: 2,423 bytes, SHA-256
`afe5e62ba3e721a29139e3bc63a018e2063273bf0a5919662a4cf4dbd31495d9`.
The sorted per-file SHA-256 listing hashes to
`ee35b00e5ecfeec9c00cda179ba6eb211a3541d92e18c22e6b3007ed8589e332`.
