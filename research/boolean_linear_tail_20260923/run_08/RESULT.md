# Specialized Boolean linear-tail workload

Correctness **PASS**: 192 cells, 56448 batch-arm samples, 3913728 oracle-verified complete affine-tail outputs.

Dramatic gate: **{'stream_plain': 'REJECTED', 'stream_exchange': 'REJECTED', 'stream_census': 'REJECTED', 'hybrid_census': 'REJECTED'}**. Gates passed: **{'stream_plain': 7, 'stream_exchange': 9, 'stream_census': 9, 'hybrid_census': 9}**, out of twelve per candidate.

Cold batch8 milliseconds include all setup, product generation, cache guards/misses, high elimination, affine canonicalization, validation and destruction. Values are medians over two holdout seeds and 42 balanced repetitions. Ratios below use pooled medians against the cached flat control; acceptance uses paired minima over both retained controls.

| Variables | Family | Variant | Cold batch (ms) | Flat / arm | Layout hits / 8 | Nonempty outputs / 8 | Max sparse terms |
|---:|---|---|---:|---:|---:|---:|---:|
| 12 | quadratic | flat_cached | 0.240820 | 1.000 | 0 | 0 | null |
| 12 | quadratic | sparse_bucket | 0.608349 | 0.396 | 0 | 0 | 113 |
| 12 | quadratic | stream_plain | 0.536059 | 0.449 | 0 | 0 | 114 |
| 12 | quadratic | stream_exchange | 0.425566 | 0.566 | 0 | 0 | 115 |
| 12 | quadratic | stream_census | 0.384154 | 0.627 | 0 | 0 | 115 |
| 12 | quadratic | hybrid_census | 0.245326 | 0.982 | 0 | 0 | null |
| 12 | quadratic | flat_alias | 0.243707 | 0.988 | 0 | 0 | null |
| 12 | linear_drop | flat_cached | 0.273579 | 1.000 | 0 | 2 | null |
| 12 | linear_drop | sparse_bucket | 0.653019 | 0.419 | 0 | 2 | 107 |
| 12 | linear_drop | stream_plain | 0.617508 | 0.443 | 0 | 2 | 104 |
| 12 | linear_drop | stream_exchange | 0.467350 | 0.585 | 0 | 2 | 110 |
| 12 | linear_drop | stream_census | 0.412801 | 0.663 | 0 | 2 | 110 |
| 12 | linear_drop | hybrid_census | 0.274898 | 0.995 | 0 | 2 | null |
| 12 | linear_drop | flat_alias | 0.273992 | 0.998 | 0 | 2 | null |
| 12 | restricted_cycle | flat_cached | 0.084224 | 1.000 | 6 | 4 | null |
| 12 | restricted_cycle | sparse_bucket | 0.677819 | 0.124 | 0 | 4 | 45 |
| 12 | restricted_cycle | stream_plain | 0.662729 | 0.127 | 0 | 4 | 46 |
| 12 | restricted_cycle | stream_exchange | 0.543107 | 0.155 | 0 | 4 | 44 |
| 12 | restricted_cycle | stream_census | 0.511134 | 0.165 | 0 | 4 | 44 |
| 12 | restricted_cycle | hybrid_census | 0.090065 | 0.935 | 6 | 4 | null |
| 12 | restricted_cycle | flat_alias | 0.088937 | 0.947 | 6 | 4 | null |
| 12 | cross_cancel | flat_cached | 0.247341 | 1.000 | 0 | 8 | null |
| 12 | cross_cancel | sparse_bucket | 0.455346 | 0.543 | 0 | 8 | 100 |
| 12 | cross_cancel | stream_plain | 0.318193 | 0.777 | 0 | 8 | 104 |
| 12 | cross_cancel | stream_exchange | 0.275796 | 0.897 | 0 | 8 | 107 |
| 12 | cross_cancel | stream_census | 0.221337 | 1.117 | 0 | 8 | 107 |
| 12 | cross_cancel | hybrid_census | 0.251200 | 0.985 | 0 | 8 | null |
| 12 | cross_cancel | flat_alias | 0.249514 | 0.991 | 0 | 8 | null |
| 20 | quadratic | flat_cached | 0.581645 | 1.000 | 0 | 0 | null |
| 20 | quadratic | sparse_bucket | 0.565697 | 1.028 | 0 | 0 | 62 |
| 20 | quadratic | stream_plain | 0.263189 | 2.210 | 0 | 0 | 60 |
| 20 | quadratic | stream_exchange | 0.254023 | 2.290 | 0 | 0 | 60 |
| 20 | quadratic | stream_census | 0.144514 | 4.025 | 0 | 0 | 60 |
| 20 | quadratic | hybrid_census | 0.140173 | 4.149 | 0 | 0 | 60 |
| 20 | quadratic | flat_alias | 0.588168 | 0.989 | 0 | 0 | null |
| 20 | linear_drop | flat_cached | 0.783143 | 1.000 | 0 | 2 | null |
| 20 | linear_drop | sparse_bucket | 0.676833 | 1.157 | 0 | 2 | 78 |
| 20 | linear_drop | stream_plain | 0.329221 | 2.379 | 0 | 2 | 84 |
| 20 | linear_drop | stream_exchange | 0.318380 | 2.460 | 0 | 2 | 78 |
| 20 | linear_drop | stream_census | 0.197846 | 3.958 | 0 | 2 | 78 |
| 20 | linear_drop | hybrid_census | 0.199758 | 3.920 | 0 | 2 | 78 |
| 20 | linear_drop | flat_alias | 0.790331 | 0.991 | 0 | 2 | null |
| 20 | restricted_cycle | flat_cached | 0.084784 | 1.000 | 6 | 4 | null |
| 20 | restricted_cycle | sparse_bucket | 0.675732 | 0.125 | 0 | 4 | 45 |
| 20 | restricted_cycle | stream_plain | 0.657181 | 0.129 | 0 | 4 | 46 |
| 20 | restricted_cycle | stream_exchange | 0.524600 | 0.162 | 0 | 4 | 44 |
| 20 | restricted_cycle | stream_census | 0.504885 | 0.168 | 0 | 4 | 44 |
| 20 | restricted_cycle | hybrid_census | 0.087539 | 0.969 | 6 | 4 | null |
| 20 | restricted_cycle | flat_alias | 0.087534 | 0.969 | 6 | 4 | null |
| 20 | cross_cancel | flat_cached | 0.619186 | 1.000 | 0 | 8 | null |
| 20 | cross_cancel | sparse_bucket | 0.581904 | 1.064 | 0 | 8 | 52 |
| 20 | cross_cancel | stream_plain | 0.236604 | 2.617 | 0 | 8 | 52 |
| 20 | cross_cancel | stream_exchange | 0.227673 | 2.720 | 0 | 8 | 50 |
| 20 | cross_cancel | stream_census | 0.145686 | 4.250 | 0 | 8 | 50 |
| 20 | cross_cancel | hybrid_census | 0.143190 | 4.324 | 0 | 8 | 50 |
| 20 | cross_cancel | flat_alias | 0.638159 | 0.970 | 0 | 8 | null |
| 28 | quadratic | flat_cached | 1.014458 | 1.000 | 0 | 0 | null |
| 28 | quadratic | sparse_bucket | 0.800138 | 1.268 | 0 | 0 | 54 |
| 28 | quadratic | stream_plain | 0.328095 | 3.092 | 0 | 0 | 54 |
| 28 | quadratic | stream_exchange | 0.320652 | 3.164 | 0 | 0 | 54 |
| 28 | quadratic | stream_census | 0.181290 | 5.596 | 0 | 0 | 54 |
| 28 | quadratic | hybrid_census | 0.177836 | 5.704 | 0 | 0 | 54 |
| 28 | quadratic | flat_alias | 1.029568 | 0.985 | 0 | 0 | null |
| 28 | linear_drop | flat_cached | 1.670395 | 1.000 | 0 | 2 | null |
| 28 | linear_drop | sparse_bucket | 0.974275 | 1.715 | 0 | 2 | 79 |
| 28 | linear_drop | stream_plain | 0.437792 | 3.816 | 0 | 2 | 98 |
| 28 | linear_drop | stream_exchange | 0.439249 | 3.803 | 0 | 2 | 77 |
| 28 | linear_drop | stream_census | 0.256724 | 6.507 | 0 | 2 | 77 |
| 28 | linear_drop | hybrid_census | 0.254964 | 6.552 | 0 | 2 | 77 |
| 28 | linear_drop | flat_alias | 1.676142 | 0.997 | 0 | 2 | null |
| 28 | restricted_cycle | flat_cached | 0.084549 | 1.000 | 6 | 4 | null |
| 28 | restricted_cycle | sparse_bucket | 0.650270 | 0.130 | 0 | 4 | 45 |
| 28 | restricted_cycle | stream_plain | 0.650654 | 0.130 | 0 | 4 | 46 |
| 28 | restricted_cycle | stream_exchange | 0.515115 | 0.164 | 0 | 4 | 44 |
| 28 | restricted_cycle | stream_census | 0.494354 | 0.171 | 0 | 4 | 44 |
| 28 | restricted_cycle | hybrid_census | 0.089932 | 0.940 | 6 | 4 | null |
| 28 | restricted_cycle | flat_alias | 0.086661 | 0.976 | 6 | 4 | null |
| 28 | cross_cancel | flat_cached | 1.088480 | 1.000 | 0 | 8 | null |
| 28 | cross_cancel | sparse_bucket | 0.779837 | 1.396 | 0 | 8 | 54 |
| 28 | cross_cancel | stream_plain | 0.330874 | 3.290 | 0 | 8 | 54 |
| 28 | cross_cancel | stream_exchange | 0.328082 | 3.318 | 0 | 8 | 54 |
| 28 | cross_cancel | stream_census | 0.189421 | 5.746 | 0 | 8 | 54 |
| 28 | cross_cancel | hybrid_census | 0.186760 | 5.828 | 0 | 8 | 54 |
| 28 | cross_cancel | flat_alias | 1.110040 | 0.981 | 0 | 8 | null |
| 36 | quadratic | flat_cached | 1.529358 | 1.000 | 0 | 0 | null |
| 36 | quadratic | sparse_bucket | 0.972990 | 1.572 | 0 | 0 | 52 |
| 36 | quadratic | stream_plain | 0.468305 | 3.266 | 0 | 0 | 52 |
| 36 | quadratic | stream_exchange | 0.453310 | 3.374 | 0 | 0 | 52 |
| 36 | quadratic | stream_census | 0.221577 | 6.902 | 0 | 0 | 52 |
| 36 | quadratic | hybrid_census | 0.216579 | 7.061 | 0 | 0 | 52 |
| 36 | quadratic | flat_alias | 1.540154 | 0.993 | 0 | 0 | null |
| 36 | linear_drop | flat_cached | 2.971131 | 1.000 | 0 | 2 | null |
| 36 | linear_drop | sparse_bucket | 1.251026 | 2.375 | 0 | 2 | 90 |
| 36 | linear_drop | stream_plain | 0.615853 | 4.824 | 0 | 2 | 98 |
| 36 | linear_drop | stream_exchange | 0.599562 | 4.955 | 0 | 2 | 95 |
| 36 | linear_drop | stream_census | 0.319836 | 9.290 | 0 | 2 | 95 |
| 36 | linear_drop | hybrid_census | 0.318159 | 9.339 | 0 | 2 | 95 |
| 36 | linear_drop | flat_alias | 2.981753 | 0.996 | 0 | 2 | null |
| 36 | restricted_cycle | flat_cached | 0.082891 | 1.000 | 6 | 4 | null |
| 36 | restricted_cycle | sparse_bucket | 0.645996 | 0.128 | 0 | 4 | 45 |
| 36 | restricted_cycle | stream_plain | 0.635158 | 0.131 | 0 | 4 | 46 |
| 36 | restricted_cycle | stream_exchange | 0.512918 | 0.162 | 0 | 4 | 44 |
| 36 | restricted_cycle | stream_census | 0.470521 | 0.176 | 0 | 4 | 44 |
| 36 | restricted_cycle | hybrid_census | 0.086267 | 0.961 | 6 | 4 | null |
| 36 | restricted_cycle | flat_alias | 0.083320 | 0.995 | 6 | 4 | null |
| 36 | cross_cancel | flat_cached | 1.628077 | 1.000 | 0 | 8 | null |
| 36 | cross_cancel | sparse_bucket | 0.961986 | 1.692 | 0 | 8 | 51 |
| 36 | cross_cancel | stream_plain | 0.443960 | 3.667 | 0 | 8 | 52 |
| 36 | cross_cancel | stream_exchange | 0.439884 | 3.701 | 0 | 8 | 51 |
| 36 | cross_cancel | stream_census | 0.232690 | 6.997 | 0 | 8 | 51 |
| 36 | cross_cancel | hybrid_census | 0.229171 | 7.104 | 0 | 8 | 51 |
| 36 | cross_cancel | flat_alias | 1.641911 | 0.992 | 0 | 8 | null |

The active mask is recomputed from each input. Restricted-cycle fixtures use only eight embedded variables; cross-cancel fixtures guarantee an affine consequence from two nonlinear generators. Exact source dimensions, high rank and canonical tail are checked. Empty and nonempty tail workloads are both retained.

Schedule/layout retention is reported as counts, not byte estimates. Whole-worker RSS includes all arms and references. Packed XORs and sparse merge items are different diagnostics, not a common complete operation unit. Independent oracle preparation is outside arm timing and inside process receipts.

The controls retain pre-existing kernel strategies; no novelty claim attaches to sparse high-column elimination or linear-tail restriction. These standalone measurements do not execute the production solver, enumerate roots, recover scalars or measure full index calculus. Those costs remain null.

This source-order run uses descending DegRevLex input terms and degree-layered multiplier enumeration, matching the retained source strategies. The earlier numeric-order run is preserved separately and is not used for a source-order performance claim.

The original per-cell 2x gates above are retained. The separately declared complete-mixture portfolio gate is **REJECTED**. It sums all twelve batch8 size/family cells at n>=20 for each paired seed/repetition and requires both cumulative improvement thresholds plus every non-regression guardrail.

{"all_prior": {"ci95_paired_median": [1.6646127428793296, 1.685798143883083], "paired_ratio_median": 1.6745189063893426, "pass": true, "threshold": 1.05}, "retained_controls": {"ci95_paired_median": [3.637541341529706, 3.7082763923793935], "paired_ratio_median": 3.6634050963167057, "pass": true, "threshold": 2.0}}

The hybrid keeps the exact flat control when the degree-bounded active monomial universe has at most 512 columns; otherwise it uses pivot-exchange streaming with an exact bitmap column census. The affine output, original source dimensions and caps are unchanged. Portfolio weights are declared benchmark weights, not measured solver-call frequencies.

Identical-control alias gate: **PASS**. Every directed predecessor/method pair, including self-pairs, occurs equally often. One unmeasured validated primer supplies the first predecessor; its cost is in process receipts. Every measured sample still creates and destroys a fresh context. Repetition pairing follows the frozen measurement-order protocol, not fixed-neighbor rotations.

Each timing observation contains 16 independent cold batches, including the primer observation. Contexts reset between them and every output is validated. Raw times/counts retain the observation totals; displayed times and per-batch counts divide by this fixed factor. Ratio gates use equal-work raw totals. The flat labels share one non-inlined leaf function.
