# P-256 Dickson-aware Riemann--Roch selector, round 296: result

## Verdict

The factor base has an exact compact global description, and the unordered
109-column selector has an exact 273-variable quadratic presentation.  The
degree gate nevertheless fails almost immediately.

The complete `m=3`, `B=5` controls solve at degree 3.  At only `m=5`, `B=8`
(13 variables), productive F4 work reaches **degree 6**, the degree-7 scan
leaves 27,628 higher pairs, and each target costs 9,994,651,362 field
operations.  The `m=7` and `m=9` cells time out with matrices as wide as
173,538 and 373,917 columns.  Only two of eight target runs are complete, so
no time, field-operation, width, or degree exponent is fitted and no P-256
cost projection is licensed.

All four planted coefficient witnesses satisfy every equation.  Exhaustive
group replay covers 2,138,272 signed distinct-column candidates, finds and
replays all 3,776 relations with zero failures, and agrees with both complete
F4 classifications.  The six incomplete classifications remain censored;
they are not counted as false positives, false negatives, or exhaustive
searches.  No P-256 relation was attempted and promotion is rejected.

## Exact P-256 formulation

The selected Round-295 union is exactly

```text
u_1 = x^2 - 2
u_{i+1} = u_i^2 - 2
(u_8 - t_8)(u_6 - t_6) = 0.
```

All 164 stored abscissae pass this one-chain quadratic circuit.  Complete
rebuilds recover the same 136-column depth-8 and 28-column depth-6 fibres,
with zero overlap.  Replacing the second factor by the tempting lifted
depth-8 terminal is not exact: that fibre has 129 liftable columns and its
union with the first component has 265 columns, **101 unwanted columns**.
The mixed-depth terminal product is necessary.

The monic support polynomial

```text
F(X) = product_j (X - x_j)
```

has degree 164 and coefficient-stream SHA-256
`ce4be042ef3622fbe9bdcf5be60f59b9ef6db30acf740c748187d0540442d5ec`.
All 164 roots evaluate to zero, all 164 derivative evaluations are nonzero,
and an independent reconstruction is byte-identical.

For odd arity 109, normalize the Riemann--Roch function by its nonzero top
coefficient.  The exact selector is

```text
F = g q
A^2 - (X^3 + aX + b) B^2 = (X - x_R) g,
```

with monic `deg(g)=109`, monic `deg(q)=55`, monic `deg(A)=55`, and
`deg(B)<=53`.  This gives 273 variables, 274 quadratic coefficient equations,
and 12,592 input monomials.  It folds permutations, forces a square-free
109-column divisor of the exact support, and lets the norm choose signs.
The unnormalized 275-by-275 presentation has a scaling orbit and is not the
measured system.

At degree 5 the generic dense Macaulay dimensions are 13,345,524,830 columns
by 949,711,400 rows, or 405,580,706,241,089,984,000 bytes at 32 bytes per
entry.  This is a generic dense comparison, not a lower bound on a specialized
subresultant solver.  Consequently it cannot by itself fail the `2^50`-byte
gate; that gate remains unset and therefore does not pass.

## Complete toy boundary table

The unit is one `F_1151` multiplication in F4 row reduction.  `D_sr` is the
highest productive solving degree.  A `>=6` entry on a censored row is an
observed lower bound, not a completed degree of regularity.

| `m,B` | target | signed candidates | exact relations / replay failures | vars / eqs | complete | `D_sr` | max columns | field operations | decision |
|:--|:--|--:|--:|--:|:--:|--:|--:|--:|:--|
| 3,5 | planted | 80 | 2 / 0 | 8 / 9 | yes | 3 | 130 | 28,224 | exact positive |
| 3,5 | unplanted | 80 | 0 / 0 | 8 / 9 | yes | 3 | 79 | 22,740 | exact negative |
| 5,8 | planted | 1,792 | 10 / 0 | 13 / 14 | no: 27,628 pairs above bound | >=6 | 30,369 | 9,994,651,362 | degree gate fails |
| 5,8 | unplanted | 1,792 | 10 / 0 | 13 / 14 | no: 27,628 pairs above bound | >=6 | 30,369 | 9,994,651,362 | censored |
| 7,11 | planted | 42,240 | 56 / 0 | 18 / 19 | timeout | >=6 | 173,538 | partial; telemetry only | censored |
| 7,11 | unplanted | 42,240 | 56 / 0 | 18 / 19 | timeout | >=6 | 173,533 | partial; telemetry only | censored |
| 9,14 | planted | 1,025,024 | 1,896 / 0 | 23 / 24 | timeout | >=6 | 373,916 | partial; telemetry only | censored |
| 9,14 | unplanted | 1,025,024 | 1,746 / 0 | 23 / 24 | timeout | >=6 | 373,917 | partial; telemetry only | censored |

The hash-selected `m=5` unplanted target happens to be the negation of the
planted target, so both have the same abscissa, relation set, and F4 system.
The preregistered target was retained rather than replaced.  The negative
degree verdict does not depend on treating those two systems as independent.

## Gates

| gate | status | evidence |
|:--|:--:|:--|
| exact common Dickson chain | pass | 164/164 rows; complete 136+28 inventories; zero overlap |
| planted coefficient witnesses | pass | four of four; every coefficient equation evaluates to zero |
| exact relation replay | pass | 3,776/3,776 relations; zero replay failures |
| zero FP/FN on complete instances | pass | complete planted/unplanted `m=3` controls agree exactly |
| three complete sizes | **fail** | one complete size; six of eight targets censored |
| structured residual degree at most 5 | **fail** | productive degree 6 already at 13 variables |
| non-increasing degree slope | **fail / unset** | no admissible three-size fit |
| projected storage below `2^50` | **fail / unset** | no specialized P-256 solver measurement |
| usable relation below `2^103` | **fail / unset** | no P-256 relation oracle cost |
| 138,031-row collection below `2^120` | **fail / unset** | no collection or rank projection |
| promoted | **no** | degree and projection gates fail |

The exact common-chain circuit and divisor encoding are retained as useful
representations.  Generic F4 on their coefficient system is rejected as the
global selector.  Any continuation must exploit polynomial factorization,
subresultants, or evaluation-code decoding directly; merely feeding these
quadratics to a generic Gröbner engine has now been measured and fails before
the 23-variable toy cell completes.

## Resources and reproducibility

The corrected canonical run was isolated and uncontended on reserved CPU 4 of
the AMD EPYC 9V74 host.  It used 145.883047 s wall, 144.085594 s user,
1.775407 s system, and 2,208,204 KiB peak RSS.  The nominal 20-second F4
deadline is cooperative; a dense reduction already in progress can return
after it, which explains the roughly 39--40-second `m=9` target times.

Deadline-censored field-operation counts depend on where cancellation is
observed.  They are preserved in `telemetry.json` but set to null in the
canonical scientific artifact.  Degree, matrix dimensions, bounded-pair
counts, references, witnesses, and gate decisions replay byte-identically.
The first pre-correction isolated run is also preserved as the first JSONL
row; the corrected canonical run is the second.

- canonical result: 14,332 bytes, SHA-256
  `52cbea224df1e1d4d6e8406f368ba413fd4e508f5309fcde65f67a98bc91c5ba`;
- semantic evidence SHA-256:
  `7bae79706405d65d3582429675f44a9ba664f5a1da6113503f944eaa6137022c`;
- canonical telemetry: 1,366 bytes, SHA-256
  `93aa37e82144e815b23059a6565cf0dc1b279885e9a449fbdcf7ea15db85b18a`;
- two-run isolation receipt: 4,851 bytes, SHA-256
  `5655ed521e0f6ec42de28084e28260ae4998a46d70818e8e1bf5133cc527b275`;
- independent result replay: byte-identical, SHA-256
  `52cbea224df1e1d4d6e8406f368ba413fd4e508f5309fcde65f67a98bc91c5ba`.

```bash
cargo test --release --bin p256_rr_selector
cargo clippy --release --bin p256_rr_selector -- -D warnings
cargo build --release --bin p256_rr_selector --bin isolated_bench
target/release/isolated_bench run --wait --cpus 4 \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_rr_selector_round296_20261007/isolation.jsonl \
  --label p256-rr-selector-round296-canonical-v2 -- \
  target/release/p256_rr_selector \
  --round295 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_dickson_union_round295_20261007/dickson-union-result.json \
  --factor-base research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_dickson_union_round295_20261007/selected-factor-base.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_rr_selector_round296_20261007/rr-selector-result.json \
  --telemetry-out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_rr_selector_round296_20261007/telemetry.json
```
