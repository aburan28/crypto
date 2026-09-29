# Initial affine restriction reduces work, but the dramatic gate is rejected

All three new treatments fail the universal >2x strongest-reference timing gate.
Each passes **0/18** dramatic comparisons. None passes every >1x incremental
comparison either. The current reference includes all 22 predecessor methods,
including all three SIMD treatments, plus the three new scalar/list controls.
The original objective remains unmet.

| Treatment | >2x groups passed | >1x groups passed | Matched-policy >1.05x gate |
|---|---:|---:|---|
| `gray_delta_simd` | 0/18 | 5/18 | REJECTED |
| `initial_simd` | 0/18 | 10/18 | PASS |
| `packed_gray16_delta_simd` | 0/18 | 0/18 | REJECTED |

For transported direct blocks and leaves, the matched control is the corresponding
retained SIMD method. Initial restriction's matched control is its direct list/
scalar implementation; its paired backend ratios span 4.389–7.250. This backend
gain does not establish a complete-method gain against the strongest reference.

The full primary contains **144 distinct generated systems and 28,224 observations**:
28 methods and seven paired repetitions. All complete with verified results, with
107 SAT and 37 UNSAT fixtures per arm. The complete 120-system predecessor grid
remains, with 24 fresh holdouts. Regression tests confirm every retained method's
models, logical work and trace remain identical on all shared inputs. Old failures
and frozen verdicts are not overwritten.

## What the mechanism did

Initial affine restriction finds rank one on all 48 cross-planted fixtures and
rank zero on the other 96. The projected-rank certificate exactly certifies the
empty affine intersection on those 96 cases, avoiding a full coefficient reduction.
Every failed certificate takes complete affine-tail extraction. This is a sound
linear-algebra shortcut, not a probabilistic hash-equality assumption.

On the fresh n24 cross-planted group, initial restriction's paired ratio is
**1.991582**, with a 95% interval **[1.835488, 2.152253]**. The lower bound misses
2.0, so the gate is rejected. Its pooled median enumerated-point count falls from
8,413,584 to 4,206,792, but setup and the competing methods still matter. These
figures are finite standalone diagnostics; no exponent claim follows.

The old n24/cross-planted/seed2097153 exception remains in the regression split.
Initial restriction has a per-fixture paired median ratio of **0.160432** against
the expanded current reference. Its direct transported-block counterpart has ratio
0.080770; the transported leaf hybrid has ratio 1.164675. Reducing the enumeration
domain does not automatically beat an alternative search order that finds a model
earlier. These current-reference ratios must not be compared as if the predecessor
had used the same reference roster.

## Cold complete standalone solve costs

Units are milliseconds, including all setup, certification and failed probes,
coefficient transformation, solving, recovery, destruction and result validation.
This table pools the two fresh n24 holdout seeds and seven repetitions. The
last numeric column is the descriptive ratio of pooled planted costs relative
to retained direct SIMD. Acceptance instead uses paired pointwise fastest-reference
costs on every designated group. All rows are **engineering**, with correctness PASS.

| Method | Planted (ms) | Cross-planted (ms) | Unplanted (ms) | Direct SIMD / method, planted | Correctness |
|---|---:|---:|---:|---:|---|
| Search-only | 58.977604 | 66.921708 | 131.455542 | 0.028 | PASS |
| Degree-3 flat kernel | 211.675750 | 154.656688 | 496.847625 | 0.008 | PASS |
| Degree-3 sparse bucket | 834.953604 | 569.018396 | 2017.926459 | 0.002 | PASS |
| Degree-3 hybrid kernel | 425.680834 | 258.226729 | 1070.070500 | 0.004 | PASS |
| Selective degree-2 flat | 69.096271 | 76.453354 | 152.386209 | 0.024 | PASS |
| Selective degree-2 one-word | 60.876854 | 67.377167 | 137.764104 | 0.027 | PASS |
| Ordered-specialization search | 51.737938 | 59.208897 | 115.901396 | 0.032 | PASS |
| Fixed-quadratic state | 21.899041 | 24.354313 | 48.403396 | 0.075 | PASS |
| Retained packed state | 13.286791 | 14.862875 | 29.163438 | 0.123 | PASS |
| Full RREF, list | 263.747000 | 189.058312 | 362.392271 | 0.006 | PASS |
| Full RREF, wide | 57.391458 | 42.091729 | 80.346542 | 0.029 | PASS |
| Echelon/affine tail, list | 103.822063 | 75.383938 | 230.920729 | 0.016 | PASS |
| Echelon/affine tail, wide | 31.455313 | 22.961083 | 64.298438 | 0.052 | PASS |
| Packed without diagnostic hashing | 9.358438 | 10.563563 | 20.912500 | 0.175 | PASS |
| Recursive affine/products, list | 201.793584 | 194.776563 | 511.947521 | 0.008 | PASS |
| Recursive affine/products, compact | 21.520312 | 22.411625 | 51.417437 | 0.076 | PASS |
| Direct Gray, scalar | 10.418708 | 10.400021 | 20.718916 | 0.157 | PASS |
| Direct Gray, SIMD | 1.638146 | 1.641104 | 3.135958 | 1.000 | PASS |
| 12-variable leaves, scalar | 11.177500 | 13.295917 | 27.136333 | 0.147 | PASS |
| 12-variable leaves, SIMD | 5.024521 | 5.933875 | 12.128438 | 0.326 | PASS |
| 16-variable leaves, scalar | 9.101875 | 10.794166 | 21.356188 | 0.180 | PASS |
| 16-variable leaves, SIMD | 1.657396 | 1.985146 | 3.804146 | 0.988 | PASS |
| Transported blocks, scalar | 10.239313 | 10.295500 | 20.634563 | 0.160 | PASS |
| Transported blocks, SIMD | 1.431812 | 1.428292 | 2.839542 | 1.144 | PASS |
| Initial restriction, list/scalar | 10.395687 | 5.240375 | 20.388834 | 0.158 | PASS |
| Initial restriction, packed/SIMD | 1.423979 | 0.732729 | 2.839938 | 1.150 | PASS |
| Transported 16-variable leaves, scalar | 9.078062 | 11.037646 | 21.386854 | 0.180 | PASS |
| Transported 16-variable leaves, SIMD | 1.592000 | 1.878562 | 3.727667 | 1.029 | PASS |

The differing statistics matter: direct SIMD wins all seven repetitions of the
fresh planted seed 20261018, while the retained 16-variable SIMD hybrid wins all
seven of seed 2359297. A favorable ratio against a pooled single-method baseline
cannot replace the paired strongest-reference gate.

[run_01/results.json](run_01/results.json) records all cells, exclusive phase costs,
paired intervals and null boundaries; [run_01/RESULT.md](run_01/RESULT.md) lists every
group comparison. [SUMMARY.json](SUMMARY.json) adds individual-case ratios and
structural counts. Confidence intervals concern these fixed fixtures, not a
population-wide guarantee.

## Validation and custody

The frozen producer passes **42 Rust tests**. Twenty Python evidence/census tests
pass, including exact replay, missing/changed data, source/model/trace corruption,
projection and phase accounting, censored values, current-SIMD reference retention,
whole-cohort position balance and all shared predecessor outcomes. The mathematical
contracts and limits are in [README.md](README.md) and [CORRECTNESS.md](CORRECTNESS.md).
These are producer checks, not an independent external audit.

Seven rotations per fixture distribute each method across every position nine times
within each size's 36-fixture cohort. This is not complete position balance within
one fixture. The full campaign took **971.803 seconds**. Whole-worker peak RSS was
**29,769,728 bytes**, including all methods and reference preparation; candidate-
specific allocation remains unmeasured. There were no censored cells.

[RUN_LEDGER.json](RUN_LEDGER.json) binds the two discovery probes, completed primary,
structural census, source lineage and summary. All executed sources and manifests
remain immutable. The primary is **ineligible for confirmation**, and no confirmation
run has been launched. Production solver cost, full index-calculus cost,
calibrated-operation ratio and rho ratio remain **null**.

## A larger conditional linear block is still unmeasured

A separate exact discovery census finds maximum independent sets in the original
quadratic interaction graph of size 2–3 at n16 and 4–5 at n24. Fixing all other
variables leaves a linear system in those selected variables. The six n24 maxima
are 4, 4, 4, 5, 5 and 4. Every selected set and the exact bounded search work are
retained in [linear_fiber_census_01/results.json](linear_fiber_census_01/results.json).

This identifies blocks of 4–32 possible assignments that might be handled by one
small changing linear system. The current SIMD baseline already evaluates sixteen
points together, so the counting factor alone is insufficient. Graph selection,
coefficient updates, all rank-deficient cases, recovery and full solve cost must
be charged. [NEXT_EXPERIMENT.md](NEXT_EXPERIMENT.md) records the prospective contract.
That solver is unimplemented and no speed ratio is assigned to the census.
