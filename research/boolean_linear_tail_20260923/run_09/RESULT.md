# Specialized Boolean linear-tail workload

Correctness **PASS**: 192 cells, 56448 batch-arm samples, 3913728 oracle-verified complete affine-tail outputs.

Dramatic gate: **{'stream_plain': 'REJECTED', 'stream_exchange': 'REJECTED', 'stream_census': 'REJECTED', 'hybrid_census': 'REJECTED'}**. Gates passed: **{'stream_plain': 8, 'stream_exchange': 8, 'stream_census': 9, 'hybrid_census': 9}**, out of twelve per candidate.

Cold batch8 milliseconds include all setup, product generation, cache guards/misses, high elimination, affine canonicalization, validation and destruction. Values are medians over two holdout seeds and 42 balanced repetitions. Ratios below use pooled medians against the cached flat control; acceptance uses paired minima over both retained controls.

| Variables | Family | Variant | Cold batch (ms) | Flat / arm | Layout hits / 8 | Nonempty outputs / 8 | Max sparse terms |
|---:|---|---|---:|---:|---:|---:|---:|
| 12 | quadratic | flat_cached | 0.228578 | 1.000 | 0 | 0 | null |
| 12 | quadratic | sparse_bucket | 0.598030 | 0.382 | 0 | 0 | 111 |
| 12 | quadratic | stream_plain | 0.566839 | 0.403 | 0 | 0 | 118 |
| 12 | quadratic | stream_exchange | 0.460211 | 0.497 | 0 | 0 | 110 |
| 12 | quadratic | stream_census | 0.407005 | 0.562 | 0 | 0 | 110 |
| 12 | quadratic | hybrid_census | 0.241355 | 0.947 | 0 | 0 | null |
| 12 | quadratic | flat_alias | 0.233410 | 0.979 | 0 | 0 | null |
| 12 | linear_drop | flat_cached | 0.257267 | 1.000 | 0 | 2 | null |
| 12 | linear_drop | sparse_bucket | 0.700512 | 0.367 | 0 | 2 | 113 |
| 12 | linear_drop | stream_plain | 0.658973 | 0.390 | 0 | 2 | 118 |
| 12 | linear_drop | stream_exchange | 0.508839 | 0.506 | 0 | 2 | 109 |
| 12 | linear_drop | stream_census | 0.467608 | 0.550 | 0 | 2 | 109 |
| 12 | linear_drop | hybrid_census | 0.261874 | 0.982 | 0 | 2 | null |
| 12 | linear_drop | flat_alias | 0.262023 | 0.982 | 0 | 2 | null |
| 12 | restricted_cycle | flat_cached | 0.075185 | 1.000 | 7 | 4 | null |
| 12 | restricted_cycle | sparse_bucket | 0.593138 | 0.127 | 0 | 4 | 40 |
| 12 | restricted_cycle | stream_plain | 0.580668 | 0.129 | 0 | 4 | 42 |
| 12 | restricted_cycle | stream_exchange | 0.464599 | 0.162 | 0 | 4 | 42 |
| 12 | restricted_cycle | stream_census | 0.440953 | 0.171 | 0 | 4 | 42 |
| 12 | restricted_cycle | hybrid_census | 0.082223 | 0.914 | 7 | 4 | null |
| 12 | restricted_cycle | flat_alias | 0.081256 | 0.925 | 7 | 4 | null |
| 12 | cross_cancel | flat_cached | 0.234991 | 1.000 | 0 | 8 | null |
| 12 | cross_cancel | sparse_bucket | 0.446462 | 0.526 | 0 | 8 | 102 |
| 12 | cross_cancel | stream_plain | 0.315112 | 0.746 | 0 | 8 | 101 |
| 12 | cross_cancel | stream_exchange | 0.276802 | 0.849 | 0 | 8 | 99 |
| 12 | cross_cancel | stream_census | 0.231100 | 1.017 | 0 | 8 | 99 |
| 12 | cross_cancel | hybrid_census | 0.237767 | 0.988 | 0 | 8 | null |
| 12 | cross_cancel | flat_alias | 0.236352 | 0.994 | 0 | 8 | null |
| 20 | quadratic | flat_cached | 0.552969 | 1.000 | 0 | 0 | null |
| 20 | quadratic | sparse_bucket | 0.536820 | 1.030 | 0 | 0 | 73 |
| 20 | quadratic | stream_plain | 0.250552 | 2.207 | 0 | 0 | 96 |
| 20 | quadratic | stream_exchange | 0.243914 | 2.267 | 0 | 0 | 90 |
| 20 | quadratic | stream_census | 0.145125 | 3.810 | 0 | 0 | 90 |
| 20 | quadratic | hybrid_census | 0.144091 | 3.838 | 0 | 0 | 90 |
| 20 | quadratic | flat_alias | 0.564451 | 0.980 | 0 | 0 | null |
| 20 | linear_drop | flat_cached | 0.749260 | 1.000 | 0 | 2 | null |
| 20 | linear_drop | sparse_bucket | 0.650163 | 1.152 | 0 | 2 | 105 |
| 20 | linear_drop | stream_plain | 0.336633 | 2.226 | 0 | 2 | 136 |
| 20 | linear_drop | stream_exchange | 0.324495 | 2.309 | 0 | 2 | 111 |
| 20 | linear_drop | stream_census | 0.207599 | 3.609 | 0 | 2 | 111 |
| 20 | linear_drop | hybrid_census | 0.204034 | 3.672 | 0 | 2 | 111 |
| 20 | linear_drop | flat_alias | 0.761225 | 0.984 | 0 | 2 | null |
| 20 | restricted_cycle | flat_cached | 0.073323 | 1.000 | 7 | 4 | null |
| 20 | restricted_cycle | sparse_bucket | 0.580780 | 0.126 | 0 | 4 | 40 |
| 20 | restricted_cycle | stream_plain | 0.570756 | 0.128 | 0 | 4 | 42 |
| 20 | restricted_cycle | stream_exchange | 0.464109 | 0.158 | 0 | 4 | 42 |
| 20 | restricted_cycle | stream_census | 0.441824 | 0.166 | 0 | 4 | 42 |
| 20 | restricted_cycle | hybrid_census | 0.081504 | 0.900 | 7 | 4 | null |
| 20 | restricted_cycle | flat_alias | 0.076009 | 0.965 | 7 | 4 | null |
| 20 | cross_cancel | flat_cached | 0.586633 | 1.000 | 0 | 8 | null |
| 20 | cross_cancel | sparse_bucket | 0.541616 | 1.083 | 0 | 8 | 66 |
| 20 | cross_cancel | stream_plain | 0.228884 | 2.563 | 0 | 8 | 74 |
| 20 | cross_cancel | stream_exchange | 0.223340 | 2.627 | 0 | 8 | 74 |
| 20 | cross_cancel | stream_census | 0.146277 | 4.010 | 0 | 8 | 74 |
| 20 | cross_cancel | hybrid_census | 0.143034 | 4.101 | 0 | 8 | 74 |
| 20 | cross_cancel | flat_alias | 0.599242 | 0.979 | 0 | 8 | null |
| 28 | quadratic | flat_cached | 0.957387 | 1.000 | 0 | 0 | null |
| 28 | quadratic | sparse_bucket | 0.754822 | 1.268 | 0 | 0 | 66 |
| 28 | quadratic | stream_plain | 0.311385 | 3.075 | 0 | 0 | 61 |
| 28 | quadratic | stream_exchange | 0.300509 | 3.186 | 0 | 0 | 61 |
| 28 | quadratic | stream_census | 0.175995 | 5.440 | 0 | 0 | 61 |
| 28 | quadratic | hybrid_census | 0.175741 | 5.448 | 0 | 0 | 61 |
| 28 | quadratic | flat_alias | 0.971100 | 0.986 | 0 | 0 | null |
| 28 | linear_drop | flat_cached | 1.588445 | 1.000 | 0 | 2 | null |
| 28 | linear_drop | sparse_bucket | 0.912363 | 1.741 | 0 | 2 | 78 |
| 28 | linear_drop | stream_plain | 0.412078 | 3.855 | 0 | 2 | 85 |
| 28 | linear_drop | stream_exchange | 0.413155 | 3.845 | 0 | 2 | 78 |
| 28 | linear_drop | stream_census | 0.246250 | 6.451 | 0 | 2 | 78 |
| 28 | linear_drop | hybrid_census | 0.247077 | 6.429 | 0 | 2 | 78 |
| 28 | linear_drop | flat_alias | 1.584027 | 1.003 | 0 | 2 | null |
| 28 | restricted_cycle | flat_cached | 0.075000 | 1.000 | 7 | 4 | null |
| 28 | restricted_cycle | sparse_bucket | 0.582980 | 0.129 | 0 | 4 | 40 |
| 28 | restricted_cycle | stream_plain | 0.572007 | 0.131 | 0 | 4 | 42 |
| 28 | restricted_cycle | stream_exchange | 0.469159 | 0.160 | 0 | 4 | 42 |
| 28 | restricted_cycle | stream_census | 0.443302 | 0.169 | 0 | 4 | 42 |
| 28 | restricted_cycle | hybrid_census | 0.080671 | 0.930 | 7 | 4 | null |
| 28 | restricted_cycle | flat_alias | 0.078835 | 0.951 | 7 | 4 | null |
| 28 | cross_cancel | flat_cached | 1.022499 | 1.000 | 0 | 8 | null |
| 28 | cross_cancel | sparse_bucket | 0.738253 | 1.385 | 0 | 8 | 66 |
| 28 | cross_cancel | stream_plain | 0.310419 | 3.294 | 0 | 8 | 61 |
| 28 | cross_cancel | stream_exchange | 0.313622 | 3.260 | 0 | 8 | 61 |
| 28 | cross_cancel | stream_census | 0.185891 | 5.501 | 0 | 8 | 61 |
| 28 | cross_cancel | hybrid_census | 0.183576 | 5.570 | 0 | 8 | 61 |
| 28 | cross_cancel | flat_alias | 1.031056 | 0.992 | 0 | 8 | null |
| 36 | quadratic | flat_cached | 1.481716 | 1.000 | 0 | 0 | null |
| 36 | quadratic | sparse_bucket | 0.943720 | 1.570 | 0 | 0 | 48 |
| 36 | quadratic | stream_plain | 0.441522 | 3.356 | 0 | 0 | 48 |
| 36 | quadratic | stream_exchange | 0.437848 | 3.384 | 0 | 0 | 48 |
| 36 | quadratic | stream_census | 0.213036 | 6.955 | 0 | 0 | 48 |
| 36 | quadratic | hybrid_census | 0.213638 | 6.936 | 0 | 0 | 48 |
| 36 | quadratic | flat_alias | 1.484120 | 0.998 | 0 | 0 | null |
| 36 | linear_drop | flat_cached | 2.955561 | 1.000 | 0 | 2 | null |
| 36 | linear_drop | sparse_bucket | 1.195530 | 2.472 | 0 | 2 | 66 |
| 36 | linear_drop | stream_plain | 0.564499 | 5.236 | 0 | 2 | 64 |
| 36 | linear_drop | stream_exchange | 0.565529 | 5.226 | 0 | 2 | 65 |
| 36 | linear_drop | stream_census | 0.302470 | 9.771 | 0 | 2 | 65 |
| 36 | linear_drop | hybrid_census | 0.301923 | 9.789 | 0 | 2 | 65 |
| 36 | linear_drop | flat_alias | 2.994671 | 0.987 | 0 | 2 | null |
| 36 | restricted_cycle | flat_cached | 0.074004 | 1.000 | 7 | 4 | null |
| 36 | restricted_cycle | sparse_bucket | 0.591923 | 0.125 | 0 | 4 | 40 |
| 36 | restricted_cycle | stream_plain | 0.579945 | 0.128 | 0 | 4 | 42 |
| 36 | restricted_cycle | stream_exchange | 0.466715 | 0.159 | 0 | 4 | 42 |
| 36 | restricted_cycle | stream_census | 0.444639 | 0.166 | 0 | 4 | 42 |
| 36 | restricted_cycle | hybrid_census | 0.082290 | 0.899 | 7 | 4 | null |
| 36 | restricted_cycle | flat_alias | 0.076327 | 0.970 | 7 | 4 | null |
| 36 | cross_cancel | flat_cached | 1.590747 | 1.000 | 0 | 8 | null |
| 36 | cross_cancel | sparse_bucket | 0.943421 | 1.686 | 0 | 8 | 48 |
| 36 | cross_cancel | stream_plain | 0.435589 | 3.652 | 0 | 8 | 48 |
| 36 | cross_cancel | stream_exchange | 0.423594 | 3.755 | 0 | 8 | 48 |
| 36 | cross_cancel | stream_census | 0.226766 | 7.015 | 0 | 8 | 48 |
| 36 | cross_cancel | hybrid_census | 0.223602 | 7.114 | 0 | 8 | 48 |
| 36 | cross_cancel | flat_alias | 1.607216 | 0.990 | 0 | 8 | null |

The active mask is recomputed from each input. Restricted-cycle fixtures use only eight embedded variables; cross-cancel fixtures guarantee an affine consequence from two nonlinear generators. Exact source dimensions, high rank and canonical tail are checked. Empty and nonempty tail workloads are both retained.

Schedule/layout retention is reported as counts, not byte estimates. Whole-worker RSS includes all arms and references. Packed XORs and sparse merge items are different diagnostics, not a common complete operation unit. Independent oracle preparation is outside arm timing and inside process receipts.

The controls retain pre-existing kernel strategies; no novelty claim attaches to sparse high-column elimination or linear-tail restriction. These standalone measurements do not execute the production solver, enumerate roots, recover scalars or measure full index calculus. Those costs remain null.

This source-order run uses descending DegRevLex input terms and degree-layered multiplier enumeration, matching the retained source strategies. The earlier numeric-order run is preserved separately and is not used for a source-order performance claim.

The original per-cell 2x gates above are retained. The separately declared complete-mixture portfolio gate is **REJECTED**. It sums all twelve batch8 size/family cells at n>=20 for each paired seed/repetition and requires both cumulative improvement thresholds plus every non-regression guardrail.

{"all_prior": {"ci95_paired_median": [1.621558677108875, 1.6573163229222105], "paired_ratio_median": 1.6403829799163987, "pass": true, "threshold": 1.05}, "retained_controls": {"ci95_paired_median": [3.5743848700754173, 3.6193571555847392], "paired_ratio_median": 3.590953663178877, "pass": true, "threshold": 2.0}}

The hybrid keeps the exact flat control when the degree-bounded active monomial universe has at most 512 columns; otherwise it uses pivot-exchange streaming with an exact bitmap column census. The affine output, original source dimensions and caps are unchanged. Portfolio weights are declared benchmark weights, not measured solver-call frequencies.

Identical-control alias gate: **PASS**. Every directed predecessor/method pair, including self-pairs, occurs equally often. One unmeasured validated primer supplies the first predecessor; its cost is in process receipts. Every measured sample still creates and destroys a fresh context. Repetition pairing follows the frozen measurement-order protocol, not fixed-neighbor rotations.

Each timing observation contains 16 independent cold batches, including the primer observation. Contexts reset between them and every output is validated. Raw times/counts retain the observation totals; displayed times and per-batch counts divide by this fixed factor. Ratio gates use equal-work raw totals. The flat labels share one non-inlined leaf function.
