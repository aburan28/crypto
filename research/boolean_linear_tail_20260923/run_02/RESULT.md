# Specialized Boolean linear-tail workload

Correctness **PASS**: 192 cells, 9216 batch-arm samples, 39936 oracle-verified complete affine-tail outputs.

Dramatic gate: **{'stream_plain': 'REJECTED', 'stream_exchange': 'REJECTED'}**. Gates passed: **{'stream_plain': 2, 'stream_exchange': 7}**, out of twelve per candidate.

Cold batch8 milliseconds include all setup, product generation, cache guards/misses, high elimination, affine canonicalization, validation and destruction. Values are medians over two holdout seeds and 12 balanced repetitions. Ratios below use pooled medians against the cached flat control; acceptance uses paired minima over both retained controls.

| Variables | Family | Variant | Cold batch (ms) | Flat / arm | Layout hits / 8 | Nonempty outputs / 8 | Max sparse terms |
|---:|---|---|---:|---:|---:|---:|---:|
| 12 | quadratic | flat_cached | 0.293938 | 1.000 | 0 | 0 | null |
| 12 | quadratic | sparse_bucket | 0.717646 | 0.410 | 0 | 0 | 112 |
| 12 | quadratic | stream_plain | 0.674354 | 0.436 | 0 | 0 | 110 |
| 12 | quadratic | stream_exchange | 0.555041 | 0.530 | 0 | 0 | 112 |
| 12 | linear_drop | flat_cached | 0.318063 | 1.000 | 0 | 2 | null |
| 12 | linear_drop | sparse_bucket | 0.787083 | 0.404 | 0 | 2 | 112 |
| 12 | linear_drop | stream_plain | 0.752937 | 0.422 | 0 | 2 | 110 |
| 12 | linear_drop | stream_exchange | 0.607750 | 0.523 | 0 | 2 | 110 |
| 12 | restricted_cycle | flat_cached | 0.134458 | 1.000 | 7 | 4 | null |
| 12 | restricted_cycle | sparse_bucket | 0.716083 | 0.188 | 0 | 4 | 42 |
| 12 | restricted_cycle | stream_plain | 0.709167 | 0.190 | 0 | 4 | 46 |
| 12 | restricted_cycle | stream_exchange | 0.593000 | 0.227 | 0 | 4 | 48 |
| 12 | cross_cancel | flat_cached | 0.290687 | 1.000 | 0 | 8 | null |
| 12 | cross_cancel | sparse_bucket | 0.545438 | 0.533 | 0 | 8 | 114 |
| 12 | cross_cancel | stream_plain | 0.419021 | 0.694 | 0 | 8 | 104 |
| 12 | cross_cancel | stream_exchange | 0.351062 | 0.828 | 0 | 8 | 100 |
| 20 | quadratic | flat_cached | 0.603667 | 1.000 | 0 | 0 | null |
| 20 | quadratic | sparse_bucket | 0.604812 | 0.998 | 0 | 0 | 78 |
| 20 | quadratic | stream_plain | 0.295896 | 2.040 | 0 | 0 | 70 |
| 20 | quadratic | stream_exchange | 0.272834 | 2.213 | 0 | 0 | 68 |
| 20 | linear_drop | flat_cached | 0.852229 | 1.000 | 0 | 2 | null |
| 20 | linear_drop | sparse_bucket | 0.762563 | 1.118 | 0 | 2 | 100 |
| 20 | linear_drop | stream_plain | 0.458646 | 1.858 | 0 | 2 | 104 |
| 20 | linear_drop | stream_exchange | 0.411875 | 2.069 | 0 | 2 | 103 |
| 20 | restricted_cycle | flat_cached | 0.136896 | 1.000 | 7 | 4 | null |
| 20 | restricted_cycle | sparse_bucket | 0.717417 | 0.191 | 0 | 4 | 42 |
| 20 | restricted_cycle | stream_plain | 0.711063 | 0.193 | 0 | 4 | 46 |
| 20 | restricted_cycle | stream_exchange | 0.587375 | 0.233 | 0 | 4 | 48 |
| 20 | cross_cancel | flat_cached | 0.685563 | 1.000 | 0 | 8 | null |
| 20 | cross_cancel | sparse_bucket | 0.643583 | 1.065 | 0 | 8 | 67 |
| 20 | cross_cancel | stream_plain | 0.303063 | 2.262 | 0 | 8 | 67 |
| 20 | cross_cancel | stream_exchange | 0.274562 | 2.497 | 0 | 8 | 67 |
| 28 | quadratic | flat_cached | 1.099687 | 1.000 | 0 | 0 | null |
| 28 | quadratic | sparse_bucket | 0.874458 | 1.258 | 0 | 0 | 59 |
| 28 | quadratic | stream_plain | 0.419938 | 2.619 | 0 | 0 | 59 |
| 28 | quadratic | stream_exchange | 0.377937 | 2.910 | 0 | 0 | 59 |
| 28 | linear_drop | flat_cached | 1.768750 | 1.000 | 0 | 2 | null |
| 28 | linear_drop | sparse_bucket | 1.082791 | 1.634 | 0 | 2 | 102 |
| 28 | linear_drop | stream_plain | 0.542750 | 3.259 | 0 | 2 | 101 |
| 28 | linear_drop | stream_exchange | 0.517917 | 3.415 | 0 | 2 | 102 |
| 28 | restricted_cycle | flat_cached | 0.128312 | 1.000 | 7 | 4 | null |
| 28 | restricted_cycle | sparse_bucket | 0.714896 | 0.179 | 0 | 4 | 42 |
| 28 | restricted_cycle | stream_plain | 0.724000 | 0.177 | 0 | 4 | 46 |
| 28 | restricted_cycle | stream_exchange | 0.578626 | 0.222 | 0 | 4 | 48 |
| 28 | cross_cancel | flat_cached | 1.147958 | 1.000 | 0 | 8 | null |
| 28 | cross_cancel | sparse_bucket | 0.844979 | 1.359 | 0 | 8 | 59 |
| 28 | cross_cancel | stream_plain | 0.397375 | 2.889 | 0 | 8 | 59 |
| 28 | cross_cancel | stream_exchange | 0.389688 | 2.946 | 0 | 8 | 59 |
| 36 | quadratic | flat_cached | 1.617730 | 1.000 | 0 | 0 | null |
| 36 | quadratic | sparse_bucket | 1.070396 | 1.511 | 0 | 0 | 50 |
| 36 | quadratic | stream_plain | 0.538791 | 3.003 | 0 | 0 | 50 |
| 36 | quadratic | stream_exchange | 0.516667 | 3.131 | 0 | 0 | 50 |
| 36 | linear_drop | flat_cached | 3.613958 | 1.000 | 0 | 2 | null |
| 36 | linear_drop | sparse_bucket | 1.402125 | 2.577 | 0 | 2 | 102 |
| 36 | linear_drop | stream_plain | 0.761667 | 4.745 | 0 | 2 | 100 |
| 36 | linear_drop | stream_exchange | 0.714355 | 5.059 | 0 | 2 | 98 |
| 36 | restricted_cycle | flat_cached | 0.137292 | 1.000 | 7 | 4 | null |
| 36 | restricted_cycle | sparse_bucket | 0.717916 | 0.191 | 0 | 4 | 42 |
| 36 | restricted_cycle | stream_plain | 0.711938 | 0.193 | 0 | 4 | 46 |
| 36 | restricted_cycle | stream_exchange | 0.601417 | 0.228 | 0 | 4 | 48 |
| 36 | cross_cancel | flat_cached | 1.700250 | 1.000 | 0 | 8 | null |
| 36 | cross_cancel | sparse_bucket | 1.035333 | 1.642 | 0 | 8 | 48 |
| 36 | cross_cancel | stream_plain | 0.500396 | 3.398 | 0 | 8 | 48 |
| 36 | cross_cancel | stream_exchange | 0.475604 | 3.575 | 0 | 8 | 48 |

The active mask is recomputed from each input. Restricted-cycle fixtures use only eight embedded variables; cross-cancel fixtures guarantee an affine consequence from two nonlinear generators. Exact source dimensions, high rank and canonical tail are checked. Empty and nonempty tail workloads are both retained.

Schedule/layout retention is reported as counts, not byte estimates. Whole-worker RSS includes all arms and references. Packed XORs and sparse merge items are different diagnostics, not a common complete operation unit. Independent oracle preparation is outside arm timing and inside process receipts.

The controls retain pre-existing kernel strategies; no novelty claim attaches to sparse high-column elimination or linear-tail restriction. These standalone measurements do not execute the production solver, enumerate roots, recover scalars or measure full index calculus. Those costs remain null.

This source-order run uses descending DegRevLex input terms and degree-layered multiplier enumeration, matching the retained source strategies. The earlier numeric-order run is preserved separately and is not used for a source-order performance claim.
