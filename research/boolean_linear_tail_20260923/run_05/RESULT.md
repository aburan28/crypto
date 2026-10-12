# Specialized Boolean linear-tail workload

Correctness **PASS**: 192 cells, 13824 batch-arm samples, 59904 oracle-verified complete affine-tail outputs.

Dramatic gate: **{'stream_plain': 'REJECTED', 'stream_exchange': 'REJECTED', 'stream_census': 'REJECTED', 'hybrid_census': 'REJECTED'}**. Gates passed: **{'stream_plain': 2, 'stream_exchange': 5, 'stream_census': 9, 'hybrid_census': 9}**, out of twelve per candidate.

Cold batch8 milliseconds include all setup, product generation, cache guards/misses, high elimination, affine canonicalization, validation and destruction. Values are medians over two holdout seeds and 12 balanced repetitions. Ratios below use pooled medians against the cached flat control; acceptance uses paired minima over both retained controls.

| Variables | Family | Variant | Cold batch (ms) | Flat / arm | Layout hits / 8 | Nonempty outputs / 8 | Max sparse terms |
|---:|---|---|---:|---:|---:|---:|---:|
| 12 | quadratic | flat_cached | 0.275251 | 1.000 | 0 | 0 | null |
| 12 | quadratic | sparse_bucket | 0.643958 | 0.427 | 0 | 0 | 112 |
| 12 | quadratic | stream_plain | 0.636063 | 0.433 | 0 | 0 | 116 |
| 12 | quadratic | stream_exchange | 0.476250 | 0.578 | 0 | 0 | 111 |
| 12 | quadratic | stream_census | 0.421959 | 0.652 | 0 | 0 | 111 |
| 12 | quadratic | hybrid_census | 0.292355 | 0.941 | 0 | 0 | null |
| 12 | linear_drop | flat_cached | 0.319375 | 1.000 | 0 | 2 | null |
| 12 | linear_drop | sparse_bucket | 0.774042 | 0.413 | 0 | 2 | 102 |
| 12 | linear_drop | stream_plain | 0.767875 | 0.416 | 0 | 2 | 112 |
| 12 | linear_drop | stream_exchange | 0.553292 | 0.577 | 0 | 2 | 108 |
| 12 | linear_drop | stream_census | 0.537542 | 0.594 | 0 | 2 | 108 |
| 12 | linear_drop | hybrid_census | 0.335020 | 0.953 | 0 | 2 | null |
| 12 | restricted_cycle | flat_cached | 0.108500 | 1.000 | 7 | 4 | null |
| 12 | restricted_cycle | sparse_bucket | 0.787896 | 0.138 | 0 | 4 | 43 |
| 12 | restricted_cycle | stream_plain | 0.763541 | 0.142 | 0 | 4 | 48 |
| 12 | restricted_cycle | stream_exchange | 0.640792 | 0.169 | 0 | 4 | 46 |
| 12 | restricted_cycle | stream_census | 0.597375 | 0.182 | 0 | 4 | 46 |
| 12 | restricted_cycle | hybrid_census | 0.139604 | 0.777 | 7 | 4 | null |
| 12 | cross_cancel | flat_cached | 0.280334 | 1.000 | 0 | 8 | null |
| 12 | cross_cancel | sparse_bucket | 0.565604 | 0.496 | 0 | 8 | 104 |
| 12 | cross_cancel | stream_plain | 0.435937 | 0.643 | 0 | 8 | 108 |
| 12 | cross_cancel | stream_exchange | 0.343209 | 0.817 | 0 | 8 | 106 |
| 12 | cross_cancel | stream_census | 0.288500 | 0.972 | 0 | 8 | 106 |
| 12 | cross_cancel | hybrid_census | 0.295063 | 0.950 | 0 | 8 | null |
| 20 | quadratic | flat_cached | 0.632459 | 1.000 | 0 | 0 | null |
| 20 | quadratic | sparse_bucket | 0.626084 | 1.010 | 0 | 0 | 93 |
| 20 | quadratic | stream_plain | 0.321145 | 1.969 | 0 | 0 | 93 |
| 20 | quadratic | stream_exchange | 0.298542 | 2.118 | 0 | 0 | 93 |
| 20 | quadratic | stream_census | 0.174313 | 3.628 | 0 | 0 | 93 |
| 20 | quadratic | hybrid_census | 0.163521 | 3.868 | 0 | 0 | 93 |
| 20 | linear_drop | flat_cached | 0.868291 | 1.000 | 0 | 2 | null |
| 20 | linear_drop | sparse_bucket | 0.791438 | 1.097 | 0 | 2 | 94 |
| 20 | linear_drop | stream_plain | 0.447958 | 1.938 | 0 | 2 | 102 |
| 20 | linear_drop | stream_exchange | 0.410042 | 2.118 | 0 | 2 | 94 |
| 20 | linear_drop | stream_census | 0.264833 | 3.279 | 0 | 2 | 94 |
| 20 | linear_drop | hybrid_census | 0.236292 | 3.675 | 0 | 2 | 94 |
| 20 | restricted_cycle | flat_cached | 0.108458 | 1.000 | 7 | 4 | null |
| 20 | restricted_cycle | sparse_bucket | 0.748479 | 0.145 | 0 | 4 | 43 |
| 20 | restricted_cycle | stream_plain | 0.748542 | 0.145 | 0 | 4 | 48 |
| 20 | restricted_cycle | stream_exchange | 0.625916 | 0.173 | 0 | 4 | 46 |
| 20 | restricted_cycle | stream_census | 0.587333 | 0.185 | 0 | 4 | 46 |
| 20 | restricted_cycle | hybrid_census | 0.133333 | 0.813 | 7 | 4 | null |
| 20 | cross_cancel | flat_cached | 0.691438 | 1.000 | 0 | 8 | null |
| 20 | cross_cancel | sparse_bucket | 0.649687 | 1.064 | 0 | 8 | 72 |
| 20 | cross_cancel | stream_plain | 0.305229 | 2.265 | 0 | 8 | 74 |
| 20 | cross_cancel | stream_exchange | 0.273250 | 2.530 | 0 | 8 | 72 |
| 20 | cross_cancel | stream_census | 0.179584 | 3.850 | 0 | 8 | 72 |
| 20 | cross_cancel | hybrid_census | 0.158688 | 4.357 | 0 | 8 | 72 |
| 28 | quadratic | flat_cached | 1.087083 | 1.000 | 0 | 0 | null |
| 28 | quadratic | sparse_bucket | 0.889749 | 1.222 | 0 | 0 | 66 |
| 28 | quadratic | stream_plain | 0.418355 | 2.598 | 0 | 0 | 64 |
| 28 | quadratic | stream_exchange | 0.388813 | 2.796 | 0 | 0 | 66 |
| 28 | quadratic | stream_census | 0.217062 | 5.008 | 0 | 0 | 66 |
| 28 | quadratic | hybrid_census | 0.206875 | 5.255 | 0 | 0 | 66 |
| 28 | linear_drop | flat_cached | 1.684792 | 1.000 | 0 | 2 | null |
| 28 | linear_drop | sparse_bucket | 1.076583 | 1.565 | 0 | 2 | 116 |
| 28 | linear_drop | stream_plain | 0.577229 | 2.919 | 0 | 2 | 104 |
| 28 | linear_drop | stream_exchange | 0.534937 | 3.150 | 0 | 2 | 90 |
| 28 | linear_drop | stream_census | 0.336792 | 5.002 | 0 | 2 | 90 |
| 28 | linear_drop | hybrid_census | 0.320958 | 5.249 | 0 | 2 | 90 |
| 28 | restricted_cycle | flat_cached | 0.104625 | 1.000 | 7 | 4 | null |
| 28 | restricted_cycle | sparse_bucket | 0.747563 | 0.140 | 0 | 4 | 43 |
| 28 | restricted_cycle | stream_plain | 0.742354 | 0.141 | 0 | 4 | 48 |
| 28 | restricted_cycle | stream_exchange | 0.630355 | 0.166 | 0 | 4 | 46 |
| 28 | restricted_cycle | stream_census | 0.587480 | 0.178 | 0 | 4 | 46 |
| 28 | restricted_cycle | hybrid_census | 0.131646 | 0.795 | 7 | 4 | null |
| 28 | cross_cancel | flat_cached | 1.141458 | 1.000 | 0 | 8 | null |
| 28 | cross_cancel | sparse_bucket | 0.841604 | 1.356 | 0 | 8 | 66 |
| 28 | cross_cancel | stream_plain | 0.419604 | 2.720 | 0 | 8 | 64 |
| 28 | cross_cancel | stream_exchange | 0.400770 | 2.848 | 0 | 8 | 66 |
| 28 | cross_cancel | stream_census | 0.223833 | 5.100 | 0 | 8 | 66 |
| 28 | cross_cancel | hybrid_census | 0.216542 | 5.271 | 0 | 8 | 66 |
| 36 | quadratic | flat_cached | 1.512146 | 1.000 | 0 | 0 | null |
| 36 | quadratic | sparse_bucket | 1.018292 | 1.485 | 0 | 0 | 51 |
| 36 | quadratic | stream_plain | 0.530271 | 2.852 | 0 | 0 | 52 |
| 36 | quadratic | stream_exchange | 0.506104 | 2.988 | 0 | 0 | 51 |
| 36 | quadratic | stream_census | 0.238584 | 6.338 | 0 | 0 | 51 |
| 36 | quadratic | hybrid_census | 0.230499 | 6.560 | 0 | 0 | 51 |
| 36 | linear_drop | flat_cached | 3.010021 | 1.000 | 0 | 2 | null |
| 36 | linear_drop | sparse_bucket | 1.304312 | 2.308 | 0 | 2 | 79 |
| 36 | linear_drop | stream_plain | 0.703937 | 4.276 | 0 | 2 | 78 |
| 36 | linear_drop | stream_exchange | 0.670166 | 4.491 | 0 | 2 | 79 |
| 36 | linear_drop | stream_census | 0.361500 | 8.326 | 0 | 2 | 79 |
| 36 | linear_drop | hybrid_census | 0.337167 | 8.927 | 0 | 2 | 79 |
| 36 | restricted_cycle | flat_cached | 0.106375 | 1.000 | 7 | 4 | null |
| 36 | restricted_cycle | sparse_bucket | 0.749750 | 0.142 | 0 | 4 | 43 |
| 36 | restricted_cycle | stream_plain | 0.739459 | 0.144 | 0 | 4 | 48 |
| 36 | restricted_cycle | stream_exchange | 0.617896 | 0.172 | 0 | 4 | 46 |
| 36 | restricted_cycle | stream_census | 0.585916 | 0.182 | 0 | 4 | 46 |
| 36 | restricted_cycle | hybrid_census | 0.133083 | 0.799 | 7 | 4 | null |
| 36 | cross_cancel | flat_cached | 1.594396 | 1.000 | 0 | 8 | null |
| 36 | cross_cancel | sparse_bucket | 1.040104 | 1.533 | 0 | 8 | 46 |
| 36 | cross_cancel | stream_plain | 0.518563 | 3.075 | 0 | 8 | 46 |
| 36 | cross_cancel | stream_exchange | 0.482834 | 3.302 | 0 | 8 | 47 |
| 36 | cross_cancel | stream_census | 0.255875 | 6.231 | 0 | 8 | 47 |
| 36 | cross_cancel | hybrid_census | 0.239063 | 6.669 | 0 | 8 | 47 |

The active mask is recomputed from each input. Restricted-cycle fixtures use only eight embedded variables; cross-cancel fixtures guarantee an affine consequence from two nonlinear generators. Exact source dimensions, high rank and canonical tail are checked. Empty and nonempty tail workloads are both retained.

Schedule/layout retention is reported as counts, not byte estimates. Whole-worker RSS includes all arms and references. Packed XORs and sparse merge items are different diagnostics, not a common complete operation unit. Independent oracle preparation is outside arm timing and inside process receipts.

The controls retain pre-existing kernel strategies; no novelty claim attaches to sparse high-column elimination or linear-tail restriction. These standalone measurements do not execute the production solver, enumerate roots, recover scalars or measure full index calculus. Those costs remain null.

This source-order run uses descending DegRevLex input terms and degree-layered multiplier enumeration, matching the retained source strategies. The earlier numeric-order run is preserved separately and is not used for a source-order performance claim.

The original per-cell 2x gates above are retained. The separately declared complete-mixture portfolio gate is **REJECTED**. It sums all twelve batch8 size/family cells at n>=20 for each paired seed/repetition and requires both cumulative improvement thresholds plus every non-regression guardrail.

{"all_prior": {"ci95_paired_median": [1.6648211412894882, 1.7486637972688686], "paired_ratio_median": 1.713883651023833, "pass": true, "threshold": 1.05}, "retained_controls": {"ci95_paired_median": [3.374842263789146, 3.4806103999249705], "paired_ratio_median": 3.4113749041078694, "pass": true, "threshold": 2.0}}

The hybrid keeps the exact flat control when the degree-bounded active monomial universe has at most 512 columns; otherwise it uses pivot-exchange streaming with an exact bitmap column census. The affine output, original source dimensions and caps are unchanged. Portfolio weights are declared benchmark weights, not measured solver-call frequencies.
