# Specialized Boolean linear-tail workload

Correctness **PASS**: 192 cells, 56448 batch-arm samples, 244608 oracle-verified complete affine-tail outputs.

Dramatic gate: **{'stream_plain': 'REJECTED', 'stream_exchange': 'REJECTED', 'stream_census': 'REJECTED', 'hybrid_census': 'REJECTED'}**. Gates passed: **{'stream_plain': 3, 'stream_exchange': 4, 'stream_census': 9, 'hybrid_census': 9}**, out of twelve per candidate.

Cold batch8 milliseconds include all setup, product generation, cache guards/misses, high elimination, affine canonicalization, validation and destruction. Values are medians over two holdout seeds and 42 balanced repetitions. Ratios below use pooled medians against the cached flat control; acceptance uses paired minima over both retained controls.

| Variables | Family | Variant | Cold batch (ms) | Flat / arm | Layout hits / 8 | Nonempty outputs / 8 | Max sparse terms |
|---:|---|---|---:|---:|---:|---:|---:|
| 12 | quadratic | flat_cached | 0.244958 | 1.000 | 0 | 0 | null |
| 12 | quadratic | sparse_bucket | 0.526208 | 0.466 | 0 | 0 | 104 |
| 12 | quadratic | stream_plain | 0.474958 | 0.516 | 0 | 0 | 108 |
| 12 | quadratic | stream_exchange | 0.387541 | 0.632 | 0 | 0 | 103 |
| 12 | quadratic | stream_census | 0.332396 | 0.737 | 0 | 0 | 103 |
| 12 | quadratic | hybrid_census | 0.261145 | 0.938 | 0 | 0 | null |
| 12 | quadratic | flat_alias | 0.255625 | 0.958 | 0 | 0 | null |
| 12 | linear_drop | flat_cached | 0.283938 | 1.000 | 0 | 2 | null |
| 12 | linear_drop | sparse_bucket | 0.665250 | 0.427 | 0 | 2 | 105 |
| 12 | linear_drop | stream_plain | 0.608291 | 0.467 | 0 | 2 | 106 |
| 12 | linear_drop | stream_exchange | 0.473355 | 0.600 | 0 | 2 | 105 |
| 12 | linear_drop | stream_census | 0.424438 | 0.669 | 0 | 2 | 105 |
| 12 | linear_drop | hybrid_census | 0.291959 | 0.973 | 0 | 2 | null |
| 12 | linear_drop | flat_alias | 0.287625 | 0.987 | 0 | 2 | null |
| 12 | restricted_cycle | flat_cached | 0.088124 | 1.000 | 7 | 4 | null |
| 12 | restricted_cycle | sparse_bucket | 0.702833 | 0.125 | 0 | 4 | 42 |
| 12 | restricted_cycle | stream_plain | 0.716375 | 0.123 | 0 | 4 | 47 |
| 12 | restricted_cycle | stream_exchange | 0.582855 | 0.151 | 0 | 4 | 44 |
| 12 | restricted_cycle | stream_census | 0.555958 | 0.159 | 0 | 4 | 44 |
| 12 | restricted_cycle | hybrid_census | 0.115501 | 0.763 | 7 | 4 | null |
| 12 | restricted_cycle | flat_alias | 0.114770 | 0.768 | 7 | 4 | null |
| 12 | cross_cancel | flat_cached | 0.244229 | 1.000 | 0 | 8 | null |
| 12 | cross_cancel | sparse_bucket | 0.450166 | 0.543 | 0 | 8 | 97 |
| 12 | cross_cancel | stream_plain | 0.314229 | 0.777 | 0 | 8 | 100 |
| 12 | cross_cancel | stream_exchange | 0.266834 | 0.915 | 0 | 8 | 97 |
| 12 | cross_cancel | stream_census | 0.224062 | 1.090 | 0 | 8 | 97 |
| 12 | cross_cancel | hybrid_census | 0.257604 | 0.948 | 0 | 8 | null |
| 12 | cross_cancel | flat_alias | 0.258584 | 0.944 | 0 | 8 | null |
| 20 | quadratic | flat_cached | 0.600750 | 1.000 | 0 | 0 | null |
| 20 | quadratic | sparse_bucket | 0.591542 | 1.016 | 0 | 0 | 100 |
| 20 | quadratic | stream_plain | 0.317084 | 1.895 | 0 | 0 | 84 |
| 20 | quadratic | stream_exchange | 0.297104 | 2.022 | 0 | 0 | 88 |
| 20 | quadratic | stream_census | 0.172625 | 3.480 | 0 | 0 | 88 |
| 20 | quadratic | hybrid_census | 0.165041 | 3.640 | 0 | 0 | 88 |
| 20 | quadratic | flat_alias | 0.620917 | 0.968 | 0 | 0 | null |
| 20 | linear_drop | flat_cached | 0.812729 | 1.000 | 0 | 2 | null |
| 20 | linear_drop | sparse_bucket | 0.724812 | 1.121 | 0 | 2 | 102 |
| 20 | linear_drop | stream_plain | 0.429208 | 1.894 | 0 | 2 | 134 |
| 20 | linear_drop | stream_exchange | 0.394229 | 2.062 | 0 | 2 | 110 |
| 20 | linear_drop | stream_census | 0.257458 | 3.157 | 0 | 2 | 110 |
| 20 | linear_drop | hybrid_census | 0.242708 | 3.349 | 0 | 2 | 110 |
| 20 | linear_drop | flat_alias | 0.820812 | 0.990 | 0 | 2 | null |
| 20 | restricted_cycle | flat_cached | 0.092146 | 1.000 | 7 | 4 | null |
| 20 | restricted_cycle | sparse_bucket | 0.715521 | 0.129 | 0 | 4 | 42 |
| 20 | restricted_cycle | stream_plain | 0.710041 | 0.130 | 0 | 4 | 47 |
| 20 | restricted_cycle | stream_exchange | 0.566209 | 0.163 | 0 | 4 | 44 |
| 20 | restricted_cycle | stream_census | 0.547687 | 0.168 | 0 | 4 | 44 |
| 20 | restricted_cycle | hybrid_census | 0.111334 | 0.828 | 7 | 4 | null |
| 20 | restricted_cycle | flat_alias | 0.109313 | 0.843 | 7 | 4 | null |
| 20 | cross_cancel | flat_cached | 0.629479 | 1.000 | 0 | 8 | null |
| 20 | cross_cancel | sparse_bucket | 0.606000 | 1.039 | 0 | 8 | 82 |
| 20 | cross_cancel | stream_plain | 0.286292 | 2.199 | 0 | 8 | 81 |
| 20 | cross_cancel | stream_exchange | 0.264209 | 2.383 | 0 | 8 | 81 |
| 20 | cross_cancel | stream_census | 0.167313 | 3.762 | 0 | 8 | 81 |
| 20 | cross_cancel | hybrid_census | 0.163062 | 3.860 | 0 | 8 | 81 |
| 20 | cross_cancel | flat_alias | 0.662041 | 0.951 | 0 | 8 | null |
| 28 | quadratic | flat_cached | 1.050875 | 1.000 | 0 | 0 | null |
| 28 | quadratic | sparse_bucket | 0.857812 | 1.225 | 0 | 0 | 52 |
| 28 | quadratic | stream_plain | 0.386208 | 2.721 | 0 | 0 | 51 |
| 28 | quadratic | stream_exchange | 0.358187 | 2.934 | 0 | 0 | 52 |
| 28 | quadratic | stream_census | 0.205812 | 5.106 | 0 | 0 | 52 |
| 28 | quadratic | hybrid_census | 0.193562 | 5.429 | 0 | 0 | 52 |
| 28 | quadratic | flat_alias | 1.066354 | 0.985 | 0 | 0 | null |
| 28 | linear_drop | flat_cached | 1.818312 | 1.000 | 0 | 2 | null |
| 28 | linear_drop | sparse_bucket | 1.070042 | 1.699 | 0 | 2 | 116 |
| 28 | linear_drop | stream_plain | 0.542959 | 3.349 | 0 | 2 | 120 |
| 28 | linear_drop | stream_exchange | 0.519979 | 3.497 | 0 | 2 | 114 |
| 28 | linear_drop | stream_census | 0.313833 | 5.794 | 0 | 2 | 114 |
| 28 | linear_drop | hybrid_census | 0.306521 | 5.932 | 0 | 2 | 114 |
| 28 | linear_drop | flat_alias | 1.843604 | 0.986 | 0 | 2 | null |
| 28 | restricted_cycle | flat_cached | 0.100709 | 1.000 | 7 | 4 | null |
| 28 | restricted_cycle | sparse_bucket | 0.700000 | 0.144 | 0 | 4 | 42 |
| 28 | restricted_cycle | stream_plain | 0.708042 | 0.142 | 0 | 4 | 47 |
| 28 | restricted_cycle | stream_exchange | 0.566396 | 0.178 | 0 | 4 | 44 |
| 28 | restricted_cycle | stream_census | 0.542354 | 0.186 | 0 | 4 | 44 |
| 28 | restricted_cycle | hybrid_census | 0.112813 | 0.893 | 7 | 4 | null |
| 28 | restricted_cycle | flat_alias | 0.100729 | 1.000 | 7 | 4 | null |
| 28 | cross_cancel | flat_cached | 1.108604 | 1.000 | 0 | 8 | null |
| 28 | cross_cancel | sparse_bucket | 0.820896 | 1.350 | 0 | 8 | 52 |
| 28 | cross_cancel | stream_plain | 0.389813 | 2.844 | 0 | 8 | 51 |
| 28 | cross_cancel | stream_exchange | 0.372125 | 2.979 | 0 | 8 | 52 |
| 28 | cross_cancel | stream_census | 0.215334 | 5.148 | 0 | 8 | 52 |
| 28 | cross_cancel | hybrid_census | 0.211042 | 5.253 | 0 | 8 | 52 |
| 28 | cross_cancel | flat_alias | 1.123250 | 0.987 | 0 | 8 | null |
| 36 | quadratic | flat_cached | 1.489813 | 1.000 | 0 | 0 | null |
| 36 | quadratic | sparse_bucket | 0.998437 | 1.492 | 0 | 0 | 50 |
| 36 | quadratic | stream_plain | 0.509563 | 2.924 | 0 | 0 | 50 |
| 36 | quadratic | stream_exchange | 0.506270 | 2.943 | 0 | 0 | 50 |
| 36 | quadratic | stream_census | 0.245167 | 6.077 | 0 | 0 | 50 |
| 36 | quadratic | hybrid_census | 0.234438 | 6.355 | 0 | 0 | 50 |
| 36 | quadratic | flat_alias | 1.503875 | 0.991 | 0 | 0 | null |
| 36 | linear_drop | flat_cached | 3.650084 | 1.000 | 0 | 2 | null |
| 36 | linear_drop | sparse_bucket | 1.382771 | 2.640 | 0 | 2 | 98 |
| 36 | linear_drop | stream_plain | 0.739062 | 4.939 | 0 | 2 | 116 |
| 36 | linear_drop | stream_exchange | 0.706604 | 5.166 | 0 | 2 | 98 |
| 36 | linear_drop | stream_census | 0.375688 | 9.716 | 0 | 2 | 98 |
| 36 | linear_drop | hybrid_census | 0.369104 | 9.889 | 0 | 2 | 98 |
| 36 | linear_drop | flat_alias | 3.661188 | 0.997 | 0 | 2 | null |
| 36 | restricted_cycle | flat_cached | 0.095604 | 1.000 | 7 | 4 | null |
| 36 | restricted_cycle | sparse_bucket | 0.710125 | 0.135 | 0 | 4 | 42 |
| 36 | restricted_cycle | stream_plain | 0.714105 | 0.134 | 0 | 4 | 47 |
| 36 | restricted_cycle | stream_exchange | 0.567771 | 0.168 | 0 | 4 | 44 |
| 36 | restricted_cycle | stream_census | 0.537729 | 0.178 | 0 | 4 | 44 |
| 36 | restricted_cycle | hybrid_census | 0.114167 | 0.837 | 7 | 4 | null |
| 36 | restricted_cycle | flat_alias | 0.104000 | 0.919 | 7 | 4 | null |
| 36 | cross_cancel | flat_cached | 1.626999 | 1.000 | 0 | 8 | null |
| 36 | cross_cancel | sparse_bucket | 0.991062 | 1.642 | 0 | 8 | 48 |
| 36 | cross_cancel | stream_plain | 0.509771 | 3.192 | 0 | 8 | 48 |
| 36 | cross_cancel | stream_exchange | 0.492604 | 3.303 | 0 | 8 | 48 |
| 36 | cross_cancel | stream_census | 0.249896 | 6.511 | 0 | 8 | 48 |
| 36 | cross_cancel | hybrid_census | 0.244855 | 6.645 | 0 | 8 | 48 |
| 36 | cross_cancel | flat_alias | 1.638354 | 0.993 | 0 | 8 | null |

The active mask is recomputed from each input. Restricted-cycle fixtures use only eight embedded variables; cross-cancel fixtures guarantee an affine consequence from two nonlinear generators. Exact source dimensions, high rank and canonical tail are checked. Empty and nonempty tail workloads are both retained.

Schedule/layout retention is reported as counts, not byte estimates. Whole-worker RSS includes all arms and references. Packed XORs and sparse merge items are different diagnostics, not a common complete operation unit. Independent oracle preparation is outside arm timing and inside process receipts.

The controls retain pre-existing kernel strategies; no novelty claim attaches to sparse high-column elimination or linear-tail restriction. These standalone measurements do not execute the production solver, enumerate roots, recover scalars or measure full index calculus. Those costs remain null.

This source-order run uses descending DegRevLex input terms and degree-layered multiplier enumeration, matching the retained source strategies. The earlier numeric-order run is preserved separately and is not used for a source-order performance claim.

The original per-cell 2x gates above are retained. The separately declared complete-mixture portfolio gate is **REJECTED**. It sums all twelve batch8 size/family cells at n>=20 for each paired seed/repetition and requires both cumulative improvement thresholds plus every non-regression guardrail.

{"all_prior": {"ci95_paired_median": [1.6207156528404738, 1.6775102489603113], "paired_ratio_median": 1.6432284542315772, "pass": true, "threshold": 1.05}, "retained_controls": {"ci95_paired_median": [3.280318965869103, 3.3879544292876353], "paired_ratio_median": 3.31507082393846, "pass": true, "threshold": 2.0}}

The hybrid keeps the exact flat control when the degree-bounded active monomial universe has at most 512 columns; otherwise it uses pivot-exchange streaming with an exact bitmap column census. The affine output, original source dimensions and caps are unchanged. Portfolio weights are declared benchmark weights, not measured solver-call frequencies.

Identical-control alias gate: **REJECTED**. Every directed predecessor/method pair, including self-pairs, occurs equally often. One unmeasured validated primer supplies the first predecessor; its cost is in process receipts. Every measured sample still creates and destroys a fresh context. Repetition indices enumerate each arm occurrence for pairing, not fixed-neighbor rotations.
