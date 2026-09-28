# Specialized Boolean linear-tail workload

Correctness **PASS**: 192 cells, 13824 batch-arm samples, 59904 oracle-verified complete affine-tail outputs.

Dramatic gate: **{'stream_plain': 'REJECTED', 'stream_exchange': 'REJECTED', 'stream_census': 'REJECTED', 'hybrid_census': 'REJECTED'}**. Gates passed: **{'stream_plain': 3, 'stream_exchange': 3, 'stream_census': 9, 'hybrid_census': 9}**, out of twelve per candidate.

Cold batch8 milliseconds include all setup, product generation, cache guards/misses, high elimination, affine canonicalization, validation and destruction. Values are medians over two holdout seeds and 12 balanced repetitions. Ratios below use pooled medians against the cached flat control; acceptance uses paired minima over both retained controls.

| Variables | Family | Variant | Cold batch (ms) | Flat / arm | Layout hits / 8 | Nonempty outputs / 8 | Max sparse terms |
|---:|---|---|---:|---:|---:|---:|---:|
| 12 | quadratic | flat_cached | 0.251041 | 1.000 | 0 | 0 | null |
| 12 | quadratic | sparse_bucket | 0.502833 | 0.499 | 0 | 0 | 98 |
| 12 | quadratic | stream_plain | 0.384667 | 0.653 | 0 | 0 | 98 |
| 12 | quadratic | stream_exchange | 0.329354 | 0.762 | 0 | 0 | 94 |
| 12 | quadratic | stream_census | 0.274479 | 0.915 | 0 | 0 | 94 |
| 12 | quadratic | hybrid_census | 0.255146 | 0.984 | 0 | 0 | null |
| 12 | linear_drop | flat_cached | 0.297729 | 1.000 | 0 | 2 | null |
| 12 | linear_drop | sparse_bucket | 0.631771 | 0.471 | 0 | 2 | 98 |
| 12 | linear_drop | stream_plain | 0.529312 | 0.562 | 0 | 2 | 100 |
| 12 | linear_drop | stream_exchange | 0.437792 | 0.680 | 0 | 2 | 100 |
| 12 | linear_drop | stream_census | 0.378375 | 0.787 | 0 | 2 | 100 |
| 12 | linear_drop | hybrid_census | 0.319416 | 0.932 | 0 | 2 | null |
| 12 | restricted_cycle | flat_cached | 0.103521 | 1.000 | 7 | 4 | null |
| 12 | restricted_cycle | sparse_bucket | 0.702583 | 0.147 | 0 | 4 | 45 |
| 12 | restricted_cycle | stream_plain | 0.722354 | 0.143 | 0 | 4 | 44 |
| 12 | restricted_cycle | stream_exchange | 0.578750 | 0.179 | 0 | 4 | 46 |
| 12 | restricted_cycle | stream_census | 0.544229 | 0.190 | 0 | 4 | 46 |
| 12 | restricted_cycle | hybrid_census | 0.133958 | 0.773 | 7 | 4 | null |
| 12 | cross_cancel | flat_cached | 0.259396 | 1.000 | 0 | 8 | null |
| 12 | cross_cancel | sparse_bucket | 0.434833 | 0.597 | 0 | 8 | 84 |
| 12 | cross_cancel | stream_plain | 0.302853 | 0.857 | 0 | 8 | 94 |
| 12 | cross_cancel | stream_exchange | 0.267625 | 0.969 | 0 | 8 | 92 |
| 12 | cross_cancel | stream_census | 0.212959 | 1.218 | 0 | 8 | 92 |
| 12 | cross_cancel | hybrid_census | 0.270230 | 0.960 | 0 | 8 | null |
| 20 | quadratic | flat_cached | 0.634542 | 1.000 | 0 | 0 | null |
| 20 | quadratic | sparse_bucket | 0.607187 | 1.045 | 0 | 0 | 68 |
| 20 | quadratic | stream_plain | 0.323354 | 1.962 | 0 | 0 | 79 |
| 20 | quadratic | stream_exchange | 0.292480 | 2.170 | 0 | 0 | 72 |
| 20 | quadratic | stream_census | 0.171604 | 3.698 | 0 | 0 | 72 |
| 20 | quadratic | hybrid_census | 0.161750 | 3.923 | 0 | 0 | 72 |
| 20 | linear_drop | flat_cached | 0.834646 | 1.000 | 0 | 2 | null |
| 20 | linear_drop | sparse_bucket | 0.721541 | 1.157 | 0 | 2 | 82 |
| 20 | linear_drop | stream_plain | 0.402062 | 2.076 | 0 | 2 | 98 |
| 20 | linear_drop | stream_exchange | 0.370480 | 2.253 | 0 | 2 | 79 |
| 20 | linear_drop | stream_census | 0.234584 | 3.558 | 0 | 2 | 79 |
| 20 | linear_drop | hybrid_census | 0.214500 | 3.891 | 0 | 2 | 79 |
| 20 | restricted_cycle | flat_cached | 0.100125 | 1.000 | 7 | 4 | null |
| 20 | restricted_cycle | sparse_bucket | 0.712896 | 0.140 | 0 | 4 | 45 |
| 20 | restricted_cycle | stream_plain | 0.715187 | 0.140 | 0 | 4 | 44 |
| 20 | restricted_cycle | stream_exchange | 0.588646 | 0.170 | 0 | 4 | 46 |
| 20 | restricted_cycle | stream_census | 0.545833 | 0.183 | 0 | 4 | 46 |
| 20 | restricted_cycle | hybrid_census | 0.133376 | 0.751 | 7 | 4 | null |
| 20 | cross_cancel | flat_cached | 0.697145 | 1.000 | 0 | 8 | null |
| 20 | cross_cancel | sparse_bucket | 0.635375 | 1.097 | 0 | 8 | 66 |
| 20 | cross_cancel | stream_plain | 0.293021 | 2.379 | 0 | 8 | 70 |
| 20 | cross_cancel | stream_exchange | 0.263271 | 2.648 | 0 | 8 | 70 |
| 20 | cross_cancel | stream_census | 0.171646 | 4.062 | 0 | 8 | 70 |
| 20 | cross_cancel | hybrid_census | 0.158521 | 4.398 | 0 | 8 | 70 |
| 28 | quadratic | flat_cached | 1.049667 | 1.000 | 0 | 0 | null |
| 28 | quadratic | sparse_bucket | 0.839562 | 1.250 | 0 | 0 | 58 |
| 28 | quadratic | stream_plain | 0.391458 | 2.681 | 0 | 0 | 58 |
| 28 | quadratic | stream_exchange | 0.366687 | 2.863 | 0 | 0 | 58 |
| 28 | quadratic | stream_census | 0.204104 | 5.143 | 0 | 0 | 58 |
| 28 | quadratic | hybrid_census | 0.188917 | 5.556 | 0 | 0 | 58 |
| 28 | linear_drop | flat_cached | 1.819605 | 1.000 | 0 | 2 | null |
| 28 | linear_drop | sparse_bucket | 1.031354 | 1.764 | 0 | 2 | 100 |
| 28 | linear_drop | stream_plain | 0.536854 | 3.389 | 0 | 2 | 90 |
| 28 | linear_drop | stream_exchange | 0.499000 | 3.647 | 0 | 2 | 90 |
| 28 | linear_drop | stream_census | 0.295480 | 6.158 | 0 | 2 | 90 |
| 28 | linear_drop | hybrid_census | 0.278375 | 6.537 | 0 | 2 | 90 |
| 28 | restricted_cycle | flat_cached | 0.109291 | 1.000 | 7 | 4 | null |
| 28 | restricted_cycle | sparse_bucket | 0.706792 | 0.155 | 0 | 4 | 45 |
| 28 | restricted_cycle | stream_plain | 0.713209 | 0.153 | 0 | 4 | 44 |
| 28 | restricted_cycle | stream_exchange | 0.595437 | 0.184 | 0 | 4 | 46 |
| 28 | restricted_cycle | stream_census | 0.546646 | 0.200 | 0 | 4 | 46 |
| 28 | restricted_cycle | hybrid_census | 0.129000 | 0.847 | 7 | 4 | null |
| 28 | cross_cancel | flat_cached | 1.107209 | 1.000 | 0 | 8 | null |
| 28 | cross_cancel | sparse_bucket | 0.819646 | 1.351 | 0 | 8 | 52 |
| 28 | cross_cancel | stream_plain | 0.387854 | 2.855 | 0 | 8 | 52 |
| 28 | cross_cancel | stream_exchange | 0.364542 | 3.037 | 0 | 8 | 52 |
| 28 | cross_cancel | stream_census | 0.210854 | 5.251 | 0 | 8 | 52 |
| 28 | cross_cancel | hybrid_census | 0.197166 | 5.616 | 0 | 8 | 52 |
| 36 | quadratic | flat_cached | 1.555937 | 1.000 | 0 | 0 | null |
| 36 | quadratic | sparse_bucket | 1.023771 | 1.520 | 0 | 0 | 50 |
| 36 | quadratic | stream_plain | 0.539709 | 2.883 | 0 | 0 | 50 |
| 36 | quadratic | stream_exchange | 0.517896 | 3.004 | 0 | 0 | 50 |
| 36 | quadratic | stream_census | 0.249104 | 6.246 | 0 | 0 | 50 |
| 36 | quadratic | hybrid_census | 0.239813 | 6.488 | 0 | 0 | 50 |
| 36 | linear_drop | flat_cached | 3.227667 | 1.000 | 0 | 2 | null |
| 36 | linear_drop | sparse_bucket | 1.307083 | 2.469 | 0 | 2 | 85 |
| 36 | linear_drop | stream_plain | 0.714854 | 4.515 | 0 | 2 | 76 |
| 36 | linear_drop | stream_exchange | 0.681021 | 4.739 | 0 | 2 | 83 |
| 36 | linear_drop | stream_census | 0.355833 | 9.071 | 0 | 2 | 83 |
| 36 | linear_drop | hybrid_census | 0.340126 | 9.490 | 0 | 2 | 83 |
| 36 | restricted_cycle | flat_cached | 0.102708 | 1.000 | 7 | 4 | null |
| 36 | restricted_cycle | sparse_bucket | 0.701896 | 0.146 | 0 | 4 | 45 |
| 36 | restricted_cycle | stream_plain | 0.710875 | 0.144 | 0 | 4 | 44 |
| 36 | restricted_cycle | stream_exchange | 0.577792 | 0.178 | 0 | 4 | 46 |
| 36 | restricted_cycle | stream_census | 0.533187 | 0.193 | 0 | 4 | 46 |
| 36 | restricted_cycle | hybrid_census | 0.133750 | 0.768 | 7 | 4 | null |
| 36 | cross_cancel | flat_cached | 1.625479 | 1.000 | 0 | 8 | null |
| 36 | cross_cancel | sparse_bucket | 1.018854 | 1.595 | 0 | 8 | 50 |
| 36 | cross_cancel | stream_plain | 0.544749 | 2.984 | 0 | 8 | 50 |
| 36 | cross_cancel | stream_exchange | 0.520479 | 3.123 | 0 | 8 | 50 |
| 36 | cross_cancel | stream_census | 0.264646 | 6.142 | 0 | 8 | 50 |
| 36 | cross_cancel | hybrid_census | 0.251146 | 6.472 | 0 | 8 | 50 |

The active mask is recomputed from each input. Restricted-cycle fixtures use only eight embedded variables; cross-cancel fixtures guarantee an affine consequence from two nonlinear generators. Exact source dimensions, high rank and canonical tail are checked. Empty and nonempty tail workloads are both retained.

Schedule/layout retention is reported as counts, not byte estimates. Whole-worker RSS includes all arms and references. Packed XORs and sparse merge items are different diagnostics, not a common complete operation unit. Independent oracle preparation is outside arm timing and inside process receipts.

The controls retain pre-existing kernel strategies; no novelty claim attaches to sparse high-column elimination or linear-tail restriction. These standalone measurements do not execute the production solver, enumerate roots, recover scalars or measure full index calculus. Those costs remain null.

This source-order run uses descending DegRevLex input terms and degree-layered multiplier enumeration, matching the retained source strategies. The earlier numeric-order run is preserved separately and is not used for a source-order performance claim.

The original per-cell 2x gates above are retained. The separately declared complete-mixture portfolio gate is **REJECTED**. It sums all twelve batch8 size/family cells at n>=20 for each paired seed/repetition and requires both cumulative improvement thresholds plus every non-regression guardrail.

{"all_prior": {"ci95_paired_median": [1.696901379111896, 1.7443256828625133], "paired_ratio_median": 1.723365851772635, "pass": true, "threshold": 1.05}, "retained_controls": {"ci95_paired_median": [3.34430785256344, 3.5329204921019244], "paired_ratio_median": 3.476217243777708, "pass": true, "threshold": 2.0}}

The hybrid keeps the exact flat control when the degree-bounded active monomial universe has at most 512 columns; otherwise it uses pivot-exchange streaming with an exact bitmap column census. The affine output, original source dimensions and caps are unchanged. Portfolio weights are declared benchmark weights, not measured solver-call frequencies.
