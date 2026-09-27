# Construction-only optimization is not supported by this cost diagnostic

The cost split does not justify implementing a constructor-only candidate for a
universal dramatic claim. Under the measured unchanged-scan assumption, six of
nine relevant 16-point groups and three of nine dispatcher groups have guarded
optimistic ceilings below 2.0. The other groups remain **INCONCLUSIVE**, including
every 64-point group: extracting that scan changed timing beyond the declared
representativeness tolerance. The guards and thresholds were not relaxed.

| Policy | Guarded groups below 2x | Not ruled out | Inconclusive |
|---|---:|---:|---:|
| 16-point scan | 6/9 | 0/9 | 3/9 |
| 64-point scan | 0/9 | 0/9 | 9/9 |
| Fixed dispatcher | 3/9 | 0/9 | 6/9 |

This is an **accounting diagnostic**, not a measured speedup or a general
impossibility theorem. No constructor optimization is implemented. The frozen
protocol's classification string is retained as `engineering diagnostic`; this
report classifies the measurement work as accounting under the repository rule.
The original dramatic-gain objective remains open.

## What was measured

All **10,560 observations over 24 discovery systems and 55 methods** complete with
verified outcomes. Forty-nine retained methods are joined by profiled and
unprofiled rebound arms for the two compiled kernels and their dispatcher.
Profile and rebound share exactly the same non-inlined scan function and builder.
Every new result matches the original policy's model, semantic work and trace.

Encoding, construction, scan and release are separate timed phases. Complete
totals also charge wrapper and result-validation work. Under a fixed scan cost,
even free construction leaves `T_new >= T_scan`, so the maximum ratio is
`T_reference/T_scan`. The optimistic diagnostic subtracts the whole positive paired
profile/rebound total difference and the largest observed clock pair from the scan
timer. Nonpositive estimates yield null ceilings. This is a deliberately favorable
conditional estimate, not an absolute timing-error guarantee.

Both profile/rebound and rebound/original 95% intervals must fit in [0.8,1.25]
before a ceiling is interpreted. The original and extracted scans can differ in
code layout and data placement. A constructor change that also changes scan cost
would be a different mechanism, outside this unchanged-scan bound.

## Relevant conditional ceilings

Intervals below describe the optimistic ratio to the pointwise fastest complete
reference. An interval below 2.0 is usable only when its comparability guards pass.
Inconclusive rows stay inconclusive even when their displayed ceiling is small.
All n12 rows and every phase value remain in the machine-readable result.

| Variables | Family | Policy | Optimistic ceiling, 95% interval | Decision |
|---:|---|---|---|---|
| 16 | cross_planted | 16 | null | INCONCLUSIVE |
| 16 | cross_planted | 64 | [0.8067, 1.7846] | INCONCLUSIVE |
| 16 | cross_planted | dispatch | [1.4786, 2.3018] | INCONCLUSIVE |
| 16 | planted | 16 | null | INCONCLUSIVE |
| 16 | planted | 64 | null | INCONCLUSIVE |
| 16 | planted | dispatch | null | INCONCLUSIVE |
| 16 | unplanted | 16 | [0.9660, 1.2713] | INCONCLUSIVE |
| 16 | unplanted | 64 | [0.9451, 1.2210] | INCONCLUSIVE |
| 16 | unplanted | dispatch | [0.9162, 1.3146] | INCONCLUSIVE |
| 20 | cross_planted | 16 | [0.7954, 0.9187] | RULED_OUT_UNDER_MEASURED_SCAN |
| 20 | cross_planted | 64 | [1.0022, 1.1385] | INCONCLUSIVE |
| 20 | cross_planted | dispatch | [0.7760, 0.9434] | RULED_OUT_UNDER_MEASURED_SCAN |
| 20 | planted | 16 | [0.8308, 0.9663] | RULED_OUT_UNDER_MEASURED_SCAN |
| 20 | planted | 64 | [1.1094, 1.2649] | INCONCLUSIVE |
| 20 | planted | dispatch | [0.8447, 0.9729] | RULED_OUT_UNDER_MEASURED_SCAN |
| 20 | unplanted | 16 | [0.6064, 0.6193] | RULED_OUT_UNDER_MEASURED_SCAN |
| 20 | unplanted | 64 | [0.8539, 0.8947] | INCONCLUSIVE |
| 20 | unplanted | dispatch | [0.5943, 0.6366] | RULED_OUT_UNDER_MEASURED_SCAN |
| 24 | cross_planted | 16 | [0.5249, 0.5396] | RULED_OUT_UNDER_MEASURED_SCAN |
| 24 | cross_planted | 64 | [0.7837, 0.8143] | INCONCLUSIVE |
| 24 | cross_planted | dispatch | [0.7832, 0.7999] | INCONCLUSIVE |
| 24 | planted | 16 | [0.5391, 0.5516] | RULED_OUT_UNDER_MEASURED_SCAN |
| 24 | planted | 64 | [0.8016, 0.8158] | INCONCLUSIVE |
| 24 | planted | dispatch | [0.8094, 0.8208] | INCONCLUSIVE |
| 24 | unplanted | 16 | [0.5282, 0.5454] | RULED_OUT_UNDER_MEASURED_SCAN |
| 24 | unplanted | 64 | [0.7860, 0.8080] | INCONCLUSIVE |
| 24 | unplanted | dispatch | [0.8041, 0.8174] | INCONCLUSIVE |

The valid dispatcher groups are all at n20; their upper ceilings are below 0.973.
The n16 comparisons are too sensitive to profiling and extraction to support a
bound. The 64-point kernel's rebound/original intervals exceed the tolerance,
including approximately 1.21–1.33 in several groups. Those observations do not
permit declaring the original 64-point constructor unhelpful.

## Complete costs on n24 discovery fixtures

Units are milliseconds per cold solve plus result validation, pooling the two
discovery seeds and eight random/reverse repetitions. The descriptive ratio uses
the retained dispatcher's planted cost. It is not a promotion statistic. Profiled
methods include their observer work; no setup phase is deducted from these totals.
All methods are shown and correctness is PASS.

| Method | Planted (ms) | Cross-planted (ms) | Unplanted (ms) | Dispatcher / method, planted | Correctness |
|---|---:|---:|---:|---:|---|
| search | 92.869396 | 93.909291 | 119.293480 | 0.007 | PASS |
| flat | 359.675958 | 121.436750 | 464.239875 | 0.002 | PASS |
| bucket | 1192.472812 | 415.240312 | 1587.468188 | 0.001 | PASS |
| hybrid | 561.672042 | 189.090563 | 799.997146 | 0.001 | PASS |
| small_flat | 108.459709 | 109.787001 | 139.168250 | 0.006 | PASS |
| word_tail | 97.579708 | 99.022750 | 125.447687 | 0.007 | PASS |
| merge_search | 81.550938 | 82.529605 | 104.025729 | 0.008 | PASS |
| quadratic_state | 33.265979 | 33.999875 | 43.020625 | 0.020 | PASS |
| packed_state | 20.445124 | 20.981834 | 26.395771 | 0.033 | PASS |
| basis_list | 266.550125 | 124.032333 | 414.505187 | 0.003 | PASS |
| basis_wide | 58.842500 | 27.874146 | 90.740146 | 0.012 | PASS |
| tail_list | 182.814979 | 56.310396 | 233.922708 | 0.004 | PASS |
| tail_wide | 53.833166 | 17.131605 | 67.274230 | 0.013 | PASS |
| packed_untraced | 14.541896 | 14.836000 | 18.862250 | 0.047 | PASS |
| affine_sl_basis_list | 195.860125 | 142.035521 | 488.362396 | 0.003 | PASS |
| affine_sl_basis_fast | 21.255646 | 13.536500 | 51.402083 | 0.032 | PASS |
| gray_scalar | 17.070687 | 17.033062 | 20.704833 | 0.040 | PASS |
| gray_simd | 2.581895 | 2.660458 | 3.165625 | 0.262 | PASS |
| packed_gray12_scalar | 18.688542 | 19.668750 | 25.630000 | 0.036 | PASS |
| packed_gray12_simd | 8.194771 | 8.753479 | 11.484187 | 0.083 | PASS |
| packed_gray16_scalar | 14.625250 | 16.118833 | 21.349541 | 0.046 | PASS |
| packed_gray16_simd | 2.590396 | 2.906417 | 3.843958 | 0.262 | PASS |
| gray_delta_scalar | 16.773521 | 16.542021 | 20.425604 | 0.040 | PASS |
| gray_delta_simd | 2.327104 | 2.326604 | 2.846458 | 0.291 | PASS |
| initial_list | 16.637687 | 8.244188 | 20.450084 | 0.041 | PASS |
| initial_simd | 2.322583 | 1.188708 | 2.813687 | 0.292 | PASS |
| packed_gray16_delta_scalar | 14.789146 | 16.062771 | 21.265729 | 0.046 | PASS |
| packed_gray16_delta_simd | 2.528708 | 2.798480 | 3.713583 | 0.268 | PASS |
| fiber_rows | 4.229188 | 3.758208 | 8.149166 | 0.160 | PASS |
| fiber_columns | 3.526687 | 3.272333 | 7.715667 | 0.192 | PASS |
| fiber_simd | 1.329166 | 1.182208 | 3.093291 | 0.510 | PASS |
| fiber_zero_simd | 1.788125 | 1.797125 | 4.555396 | 0.379 | PASS |
| gray_quiet | 2.212813 | 2.296750 | 2.633729 | 0.306 | PASS |
| wide64_quiet | 1.446834 | 1.450521 | 1.774126 | 0.468 | PASS |
| byte_scalar | 5.831583 | 6.254875 | 7.116667 | 0.116 | PASS |
| byte_simd | 1.524188 | 1.431355 | 1.916980 | 0.445 | PASS |
| leaf16_quiet | 2.284250 | 2.644229 | 3.345354 | 0.297 | PASS |
| leaf16_byte_scalar | 5.504105 | 5.863208 | 8.037125 | 0.123 | PASS |
| leaf16_byte_simd | 1.768625 | 1.666062 | 2.664292 | 0.383 | PASS |
| byte_single_quiet | 1.819042 | 1.706916 | 2.322625 | 0.373 | PASS |
| leaf16_single_quiet | 2.050604 | 1.940000 | 3.063292 | 0.330 | PASS |
| byte_planes | 1.605292 | 1.545000 | 2.002667 | 0.422 | PASS |
| leaf16_byte_planes | 1.832104 | 1.832229 | 2.836500 | 0.370 | PASS |
| byte_unrolled | 1.246729 | 1.149334 | 1.514771 | 0.544 | PASS |
| wide64_unrolled | 0.692105 | 0.684729 | 0.831125 | 0.979 | PASS |
| leaf16_byte_unrolled | 1.488333 | 1.439625 | 2.257500 | 0.455 | PASS |
| word16_unrolled | 1.073771 | 1.105145 | 1.307125 | 0.631 | PASS |
| word_dispatch | 0.677625 | 0.680750 | 0.843646 | 1.000 | PASS |
| leaf16_word_unrolled | 1.381709 | 1.554812 | 2.051083 | 0.490 | PASS |
| construction16_rebound | 1.261729 | 1.288312 | 1.555959 | 0.537 | PASS |
| construction16_profile | 1.278917 | 1.298083 | 1.554917 | 0.530 | PASS |
| construction64_rebound | 0.863166 | 0.850355 | 1.060958 | 0.785 | PASS |
| construction64_profile | 0.854437 | 0.857854 | 1.044750 | 0.793 | PASS |
| construction_dispatch_rebound | 0.867812 | 0.866667 | 1.055271 | 0.781 | PASS |
| construction_dispatch_profile | 0.858792 | 0.911833 | 1.045521 | 0.789 | PASS |

The campaign took **201.581 seconds**. Whole-worker peak RSS was **19,873,792 bytes**,
including every method, reference preparation and clock calibration. Candidate-
specific memory and calibrated-operation costs remain unmeasured.

## Preserved analysis failure and validation

The original analyzer stopped on a dispatch-arm naming mismatch after every worker
completed. Its source and failure remain in the sealed `run_01` bundle. The additive
`analysis_01` correction binds that input manifest and fixes only the dispatch
name and output location. It changes no timed Rust source, observations, formulas,
guards, thresholds or cohorts. [ANALYSIS_CORRECTION.md](ANALYSIS_CORRECTION.md)
documents the execution and replay paths.

The producer passes **63 Rust tests**. Twenty Python evidence/census tests pass,
including exact replay, original-source retention, source/work corruption, null
ceilings, failed comparability guards, censored outcomes and the preserved original
analysis failure. These are producer checks, not external review. The experiment
uses discovery inputs only and has no holdout or performance promotion.

## A different mathematical degree of freedom

An exact census projects equation-coordinate words along the span of the quadratic
coefficients internal to a chosen low-variable block. The projected system is affine
in that block for every outside assignment, even when the original interaction graph
has edges inside it. Every original solution survives the projection.

The converse is false: projecting `x*y+1` can produce the zero system. A complete
method must enumerate any remaining affine freedom and verify every candidate on
all original equations. The tests explicitly preserve this necessary/sufficient
distinction. The census contains 72 configurations on the same discovery inputs:

| Variables | Low-block size | Observed annihilated ranks | Equation-quotient dimensions |
|---:|---:|---|---|
| 12 | 4 | 5, 6 | 8, 9 |
| 12 | 5 | 9, 10 | 4, 5 |
| 12 | 6 | 11, 12, 13, 14 | 0, 1, 2, 3 |
| 16 | 4 | 5, 6 | 12, 13 |
| 16 | 5 | 7, 8, 9 | 9, 10, 11 |
| 16 | 6 | 10, 11, 12, 13 | 5, 6, 7, 8 |
| 20 | 4 | 4, 5, 6 | 16, 17, 18 |
| 20 | 5 | 8, 10 | 12, 14 |
| 20 | 6 | 10, 12, 13 | 9, 10, 12 |
| 24 | 4 | 3, 4, 5 | 21, 22, 23 |
| 24 | 5 | 6, 8 | 18, 20 |
| 24 | 6 | 9, 12 | 14, 17 |

These are exact quotient dimensions, not independent-relation counts or predicted
success rates. The complete projected-fiber solver and its cost are unimplemented.
[NEXT_EXPERIMENT.md](NEXT_EXPERIMENT.md) records the required original-equation
checks, rank-deficient cases, setup accounting and matched controls.

[RUN_LEDGER.json](RUN_LEDGER.json) binds the original execution, corrected analysis,
structural census and report. Production solver, full index-calculus, calibrated-
operation and rho costs remain **null**. Neither this cost bound nor the structural
census establishes an asymptotic or cryptanalytic result.
