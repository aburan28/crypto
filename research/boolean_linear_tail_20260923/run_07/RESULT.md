# Specialized Boolean linear-tail workload

Correctness **PASS**: 192 cells, 56448 batch-arm samples, 244608 oracle-verified complete affine-tail outputs.

Dramatic gate: **{'stream_plain': 'REJECTED', 'stream_exchange': 'REJECTED', 'stream_census': 'REJECTED', 'hybrid_census': 'REJECTED'}**. Gates passed: **{'stream_plain': 3, 'stream_exchange': 5, 'stream_census': 9, 'hybrid_census': 9}**, out of twelve per candidate.

Cold batch8 milliseconds include all setup, product generation, cache guards/misses, high elimination, affine canonicalization, validation and destruction. Values are medians over two holdout seeds and 42 balanced repetitions. Ratios below use pooled medians against the cached flat control; acceptance uses paired minima over both retained controls.

| Variables | Family | Variant | Cold batch (ms) | Flat / arm | Layout hits / 8 | Nonempty outputs / 8 | Max sparse terms |
|---:|---|---|---:|---:|---:|---:|---:|
| 12 | quadratic | flat_cached | 0.237313 | 1.000 | 0 | 0 | null |
| 12 | quadratic | sparse_bucket | 0.545354 | 0.435 | 0 | 0 | 103 |
| 12 | quadratic | stream_plain | 0.498000 | 0.477 | 0 | 0 | 106 |
| 12 | quadratic | stream_exchange | 0.400292 | 0.593 | 0 | 0 | 102 |
| 12 | quadratic | stream_census | 0.352980 | 0.672 | 0 | 0 | 102 |
| 12 | quadratic | hybrid_census | 0.249541 | 0.951 | 0 | 0 | null |
| 12 | quadratic | flat_alias | 0.246917 | 0.961 | 0 | 0 | null |
| 12 | linear_drop | flat_cached | 0.281854 | 1.000 | 0 | 2 | null |
| 12 | linear_drop | sparse_bucket | 0.624854 | 0.451 | 0 | 2 | 102 |
| 12 | linear_drop | stream_plain | 0.592146 | 0.476 | 0 | 2 | 104 |
| 12 | linear_drop | stream_exchange | 0.465250 | 0.606 | 0 | 2 | 105 |
| 12 | linear_drop | stream_census | 0.411438 | 0.685 | 0 | 2 | 105 |
| 12 | linear_drop | hybrid_census | 0.292895 | 0.962 | 0 | 2 | null |
| 12 | linear_drop | flat_alias | 0.284417 | 0.991 | 0 | 2 | null |
| 12 | restricted_cycle | flat_cached | 0.101646 | 1.000 | 6 | 4 | null |
| 12 | restricted_cycle | sparse_bucket | 0.694209 | 0.146 | 0 | 4 | 44 |
| 12 | restricted_cycle | stream_plain | 0.664438 | 0.153 | 0 | 4 | 45 |
| 12 | restricted_cycle | stream_exchange | 0.565855 | 0.180 | 0 | 4 | 44 |
| 12 | restricted_cycle | stream_census | 0.536645 | 0.189 | 0 | 4 | 44 |
| 12 | restricted_cycle | hybrid_census | 0.108396 | 0.938 | 6 | 4 | null |
| 12 | restricted_cycle | flat_alias | 0.103188 | 0.985 | 6 | 4 | null |
| 12 | cross_cancel | flat_cached | 0.243355 | 1.000 | 0 | 8 | null |
| 12 | cross_cancel | sparse_bucket | 0.457145 | 0.532 | 0 | 8 | 96 |
| 12 | cross_cancel | stream_plain | 0.368833 | 0.660 | 0 | 8 | 111 |
| 12 | cross_cancel | stream_exchange | 0.296688 | 0.820 | 0 | 8 | 107 |
| 12 | cross_cancel | stream_census | 0.239500 | 1.016 | 0 | 8 | 107 |
| 12 | cross_cancel | hybrid_census | 0.254146 | 0.958 | 0 | 8 | null |
| 12 | cross_cancel | flat_alias | 0.252105 | 0.965 | 0 | 8 | null |
| 20 | quadratic | flat_cached | 0.596584 | 1.000 | 0 | 0 | null |
| 20 | quadratic | sparse_bucket | 0.572563 | 1.042 | 0 | 0 | 71 |
| 20 | quadratic | stream_plain | 0.297771 | 2.003 | 0 | 0 | 71 |
| 20 | quadratic | stream_exchange | 0.269813 | 2.211 | 0 | 0 | 71 |
| 20 | quadratic | stream_census | 0.157542 | 3.787 | 0 | 0 | 71 |
| 20 | quadratic | hybrid_census | 0.151958 | 3.926 | 0 | 0 | 71 |
| 20 | quadratic | flat_alias | 0.595563 | 1.002 | 0 | 0 | null |
| 20 | linear_drop | flat_cached | 0.786250 | 1.000 | 0 | 2 | null |
| 20 | linear_drop | sparse_bucket | 0.667146 | 1.179 | 0 | 2 | 82 |
| 20 | linear_drop | stream_plain | 0.355813 | 2.210 | 0 | 2 | 81 |
| 20 | linear_drop | stream_exchange | 0.323584 | 2.430 | 0 | 2 | 76 |
| 20 | linear_drop | stream_census | 0.206292 | 3.811 | 0 | 2 | 76 |
| 20 | linear_drop | hybrid_census | 0.202937 | 3.874 | 0 | 2 | 76 |
| 20 | linear_drop | flat_alias | 0.794479 | 0.990 | 0 | 2 | null |
| 20 | restricted_cycle | flat_cached | 0.100897 | 1.000 | 6 | 4 | null |
| 20 | restricted_cycle | sparse_bucket | 0.705979 | 0.143 | 0 | 4 | 44 |
| 20 | restricted_cycle | stream_plain | 0.671292 | 0.150 | 0 | 4 | 45 |
| 20 | restricted_cycle | stream_exchange | 0.558542 | 0.181 | 0 | 4 | 44 |
| 20 | restricted_cycle | stream_census | 0.537917 | 0.188 | 0 | 4 | 44 |
| 20 | restricted_cycle | hybrid_census | 0.112229 | 0.899 | 6 | 4 | null |
| 20 | restricted_cycle | flat_alias | 0.106833 | 0.944 | 6 | 4 | null |
| 20 | cross_cancel | flat_cached | 0.632292 | 1.000 | 0 | 8 | null |
| 20 | cross_cancel | sparse_bucket | 0.574625 | 1.100 | 0 | 8 | 61 |
| 20 | cross_cancel | stream_plain | 0.271792 | 2.326 | 0 | 8 | 61 |
| 20 | cross_cancel | stream_exchange | 0.254833 | 2.481 | 0 | 8 | 61 |
| 20 | cross_cancel | stream_census | 0.160166 | 3.948 | 0 | 8 | 61 |
| 20 | cross_cancel | hybrid_census | 0.156459 | 4.041 | 0 | 8 | 61 |
| 20 | cross_cancel | flat_alias | 0.658875 | 0.960 | 0 | 8 | null |
| 28 | quadratic | flat_cached | 1.017105 | 1.000 | 0 | 0 | null |
| 28 | quadratic | sparse_bucket | 0.813230 | 1.251 | 0 | 0 | 51 |
| 28 | quadratic | stream_plain | 0.368208 | 2.762 | 0 | 0 | 50 |
| 28 | quadratic | stream_exchange | 0.363855 | 2.795 | 0 | 0 | 51 |
| 28 | quadratic | stream_census | 0.199688 | 5.093 | 0 | 0 | 51 |
| 28 | quadratic | hybrid_census | 0.198063 | 5.135 | 0 | 0 | 51 |
| 28 | quadratic | flat_alias | 1.035292 | 0.982 | 0 | 0 | null |
| 28 | linear_drop | flat_cached | 1.544208 | 1.000 | 0 | 2 | null |
| 28 | linear_drop | sparse_bucket | 0.928250 | 1.664 | 0 | 2 | 66 |
| 28 | linear_drop | stream_plain | 0.486125 | 3.177 | 0 | 2 | 68 |
| 28 | linear_drop | stream_exchange | 0.452021 | 3.416 | 0 | 2 | 66 |
| 28 | linear_drop | stream_census | 0.267521 | 5.772 | 0 | 2 | 66 |
| 28 | linear_drop | hybrid_census | 0.259771 | 5.944 | 0 | 2 | 66 |
| 28 | linear_drop | flat_alias | 1.553792 | 0.994 | 0 | 2 | null |
| 28 | restricted_cycle | flat_cached | 0.095937 | 1.000 | 6 | 4 | null |
| 28 | restricted_cycle | sparse_bucket | 0.703146 | 0.136 | 0 | 4 | 44 |
| 28 | restricted_cycle | stream_plain | 0.673438 | 0.142 | 0 | 4 | 45 |
| 28 | restricted_cycle | stream_exchange | 0.556937 | 0.172 | 0 | 4 | 44 |
| 28 | restricted_cycle | stream_census | 0.534813 | 0.179 | 0 | 4 | 44 |
| 28 | restricted_cycle | hybrid_census | 0.108750 | 0.882 | 6 | 4 | null |
| 28 | restricted_cycle | flat_alias | 0.101479 | 0.945 | 6 | 4 | null |
| 28 | cross_cancel | flat_cached | 1.074270 | 1.000 | 0 | 8 | null |
| 28 | cross_cancel | sparse_bucket | 0.770583 | 1.394 | 0 | 8 | 48 |
| 28 | cross_cancel | stream_plain | 0.366395 | 2.932 | 0 | 8 | 47 |
| 28 | cross_cancel | stream_exchange | 0.340229 | 3.157 | 0 | 8 | 48 |
| 28 | cross_cancel | stream_census | 0.201375 | 5.335 | 0 | 8 | 48 |
| 28 | cross_cancel | hybrid_census | 0.200520 | 5.357 | 0 | 8 | 48 |
| 28 | cross_cancel | flat_alias | 1.104521 | 0.973 | 0 | 8 | null |
| 36 | quadratic | flat_cached | 1.454771 | 1.000 | 0 | 0 | null |
| 36 | quadratic | sparse_bucket | 0.947000 | 1.536 | 0 | 0 | 50 |
| 36 | quadratic | stream_plain | 0.491229 | 2.961 | 0 | 0 | 50 |
| 36 | quadratic | stream_exchange | 0.472959 | 3.076 | 0 | 0 | 50 |
| 36 | quadratic | stream_census | 0.228104 | 6.378 | 0 | 0 | 50 |
| 36 | quadratic | hybrid_census | 0.221709 | 6.562 | 0 | 0 | 50 |
| 36 | quadratic | flat_alias | 1.462625 | 0.995 | 0 | 0 | null |
| 36 | linear_drop | flat_cached | 3.119229 | 1.000 | 0 | 2 | null |
| 36 | linear_drop | sparse_bucket | 1.214292 | 2.569 | 0 | 2 | 78 |
| 36 | linear_drop | stream_plain | 0.649958 | 4.799 | 0 | 2 | 73 |
| 36 | linear_drop | stream_exchange | 0.635563 | 4.908 | 0 | 2 | 78 |
| 36 | linear_drop | stream_census | 0.330500 | 9.438 | 0 | 2 | 78 |
| 36 | linear_drop | hybrid_census | 0.324458 | 9.614 | 0 | 2 | 78 |
| 36 | linear_drop | flat_alias | 3.134375 | 0.995 | 0 | 2 | null |
| 36 | restricted_cycle | flat_cached | 0.103105 | 1.000 | 6 | 4 | null |
| 36 | restricted_cycle | sparse_bucket | 0.695000 | 0.148 | 0 | 4 | 44 |
| 36 | restricted_cycle | stream_plain | 0.671667 | 0.154 | 0 | 4 | 45 |
| 36 | restricted_cycle | stream_exchange | 0.562021 | 0.183 | 0 | 4 | 44 |
| 36 | restricted_cycle | stream_census | 0.527541 | 0.195 | 0 | 4 | 44 |
| 36 | restricted_cycle | hybrid_census | 0.111500 | 0.925 | 6 | 4 | null |
| 36 | restricted_cycle | flat_alias | 0.104480 | 0.987 | 6 | 4 | null |
| 36 | cross_cancel | flat_cached | 1.540229 | 1.000 | 0 | 8 | null |
| 36 | cross_cancel | sparse_bucket | 0.945271 | 1.629 | 0 | 8 | 44 |
| 36 | cross_cancel | stream_plain | 0.472813 | 3.258 | 0 | 8 | 44 |
| 36 | cross_cancel | stream_exchange | 0.454729 | 3.387 | 0 | 8 | 44 |
| 36 | cross_cancel | stream_census | 0.240521 | 6.404 | 0 | 8 | 44 |
| 36 | cross_cancel | hybrid_census | 0.236250 | 6.519 | 0 | 8 | 44 |
| 36 | cross_cancel | flat_alias | 1.555749 | 0.990 | 0 | 8 | null |

The active mask is recomputed from each input. Restricted-cycle fixtures use only eight embedded variables; cross-cancel fixtures guarantee an affine consequence from two nonlinear generators. Exact source dimensions, high rank and canonical tail are checked. Empty and nonempty tail workloads are both retained.

Schedule/layout retention is reported as counts, not byte estimates. Whole-worker RSS includes all arms and references. Packed XORs and sparse merge items are different diagnostics, not a common complete operation unit. Independent oracle preparation is outside arm timing and inside process receipts.

The controls retain pre-existing kernel strategies; no novelty claim attaches to sparse high-column elimination or linear-tail restriction. These standalone measurements do not execute the production solver, enumerate roots, recover scalars or measure full index calculus. Those costs remain null.

This source-order run uses descending DegRevLex input terms and degree-layered multiplier enumeration, matching the retained source strategies. The earlier numeric-order run is preserved separately and is not used for a source-order performance claim.

The original per-cell 2x gates above are retained. The separately declared complete-mixture portfolio gate is **REJECTED**. It sums all twelve batch8 size/family cells at n>=20 for each paired seed/repetition and requires both cumulative improvement thresholds plus every non-regression guardrail.

{"all_prior": {"ci95_paired_median": [1.6288862465695337, 1.6735277264229809], "paired_ratio_median": 1.64375224089435, "pass": true, "threshold": 1.05}, "retained_controls": {"ci95_paired_median": [3.312836794546807, 3.397353957450502], "paired_ratio_median": 3.346282505713514, "pass": true, "threshold": 2.0}}

The hybrid keeps the exact flat control when the degree-bounded active monomial universe has at most 512 columns; otherwise it uses pivot-exchange streaming with an exact bitmap column census. The affine output, original source dimensions and caps are unchanged. Portfolio weights are declared benchmark weights, not measured solver-call frequencies.

Identical-control alias gate: **REJECTED**. Every directed predecessor/method pair, including self-pairs, occurs equally often. One unmeasured validated primer supplies the first predecessor; its cost is in process receipts. Every measured sample still creates and destroys a fresh context. Repetition pairing follows the frozen measurement-order protocol, not fixed-neighbor rotations.
