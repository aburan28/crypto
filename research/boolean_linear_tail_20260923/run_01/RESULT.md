# Specialized Boolean linear-tail workload

Correctness **PASS**: 192 cells, 9216 batch-arm samples, 39936 oracle-verified complete affine-tail outputs.

Dramatic gate: **{'stream_plain': 'REJECTED', 'stream_exchange': 'REJECTED'}**. Gates passed: **{'stream_plain': 0, 'stream_exchange': 1}**, out of twelve per candidate.

Cold batch8 milliseconds include all setup, product generation, cache guards/misses, high elimination, affine canonicalization, validation and destruction. Values are medians over two holdout seeds and 12 balanced repetitions. Ratios below use pooled medians against the cached flat control; acceptance uses paired minima over both retained controls.

| Variables | Family | Variant | Cold batch (ms) | Flat / arm | Layout hits / 8 | Nonempty outputs / 8 | Max sparse terms |
|---:|---|---|---:|---:|---:|---:|---:|
| 12 | quadratic | flat_cached | 0.278021 | 1.000 | 0 | 0 | null |
| 12 | quadratic | sparse_bucket | 0.646479 | 0.430 | 0 | 0 | 105 |
| 12 | quadratic | stream_plain | 0.628896 | 0.442 | 0 | 0 | 104 |
| 12 | quadratic | stream_exchange | 0.524813 | 0.530 | 0 | 0 | 106 |
| 12 | linear_drop | flat_cached | 0.308458 | 1.000 | 0 | 2 | null |
| 12 | linear_drop | sparse_bucket | 0.742896 | 0.415 | 0 | 2 | 105 |
| 12 | linear_drop | stream_plain | 0.710167 | 0.434 | 0 | 2 | 104 |
| 12 | linear_drop | stream_exchange | 0.590792 | 0.522 | 0 | 2 | 102 |
| 12 | restricted_cycle | flat_cached | 0.122521 | 1.000 | 7 | 4 | null |
| 12 | restricted_cycle | sparse_bucket | 0.765208 | 0.160 | 0 | 4 | 42 |
| 12 | restricted_cycle | stream_plain | 0.739917 | 0.166 | 0 | 4 | 46 |
| 12 | restricted_cycle | stream_exchange | 0.624687 | 0.196 | 0 | 4 | 45 |
| 12 | cross_cancel | flat_cached | 0.267688 | 1.000 | 0 | 8 | null |
| 12 | cross_cancel | sparse_bucket | 0.507354 | 0.528 | 0 | 8 | 96 |
| 12 | cross_cancel | stream_plain | 0.391730 | 0.683 | 0 | 8 | 96 |
| 12 | cross_cancel | stream_exchange | 0.349417 | 0.766 | 0 | 8 | 98 |
| 20 | quadratic | flat_cached | 0.610770 | 1.000 | 0 | 0 | null |
| 20 | quadratic | sparse_bucket | 0.599938 | 1.018 | 0 | 0 | 92 |
| 20 | quadratic | stream_plain | 0.351729 | 1.736 | 0 | 0 | 92 |
| 20 | quadratic | stream_exchange | 0.335041 | 1.823 | 0 | 0 | 90 |
| 20 | linear_drop | flat_cached | 0.826416 | 1.000 | 0 | 2 | null |
| 20 | linear_drop | sparse_bucket | 0.730271 | 1.132 | 0 | 2 | 102 |
| 20 | linear_drop | stream_plain | 0.462500 | 1.787 | 0 | 2 | 112 |
| 20 | linear_drop | stream_exchange | 0.410250 | 2.014 | 0 | 2 | 99 |
| 20 | restricted_cycle | flat_cached | 0.124771 | 1.000 | 7 | 4 | null |
| 20 | restricted_cycle | sparse_bucket | 0.752625 | 0.166 | 0 | 4 | 42 |
| 20 | restricted_cycle | stream_plain | 0.745645 | 0.167 | 0 | 4 | 46 |
| 20 | restricted_cycle | stream_exchange | 0.620708 | 0.201 | 0 | 4 | 45 |
| 20 | cross_cancel | flat_cached | 0.669792 | 1.000 | 0 | 8 | null |
| 20 | cross_cancel | sparse_bucket | 0.618417 | 1.083 | 0 | 8 | 79 |
| 20 | cross_cancel | stream_plain | 0.326771 | 2.050 | 0 | 8 | 79 |
| 20 | cross_cancel | stream_exchange | 0.293646 | 2.281 | 0 | 8 | 79 |
| 28 | quadratic | flat_cached | 1.089479 | 1.000 | 0 | 0 | null |
| 28 | quadratic | sparse_bucket | 0.835208 | 1.304 | 0 | 0 | 58 |
| 28 | quadratic | stream_plain | 0.423271 | 2.574 | 0 | 0 | 58 |
| 28 | quadratic | stream_exchange | 0.419874 | 2.595 | 0 | 0 | 58 |
| 28 | linear_drop | flat_cached | 1.830958 | 1.000 | 0 | 2 | null |
| 28 | linear_drop | sparse_bucket | 1.038563 | 1.763 | 0 | 2 | 82 |
| 28 | linear_drop | stream_plain | 0.585792 | 3.126 | 0 | 2 | 86 |
| 28 | linear_drop | stream_exchange | 0.556083 | 3.293 | 0 | 2 | 82 |
| 28 | restricted_cycle | flat_cached | 0.121750 | 1.000 | 7 | 4 | null |
| 28 | restricted_cycle | sparse_bucket | 0.752896 | 0.162 | 0 | 4 | 42 |
| 28 | restricted_cycle | stream_plain | 0.735375 | 0.166 | 0 | 4 | 46 |
| 28 | restricted_cycle | stream_exchange | 0.611437 | 0.199 | 0 | 4 | 45 |
| 28 | cross_cancel | flat_cached | 1.130667 | 1.000 | 0 | 8 | null |
| 28 | cross_cancel | sparse_bucket | 0.815875 | 1.386 | 0 | 8 | 58 |
| 28 | cross_cancel | stream_plain | 0.443833 | 2.548 | 0 | 8 | 58 |
| 28 | cross_cancel | stream_exchange | 0.412125 | 2.744 | 0 | 8 | 58 |
| 36 | quadratic | flat_cached | 1.516354 | 1.000 | 0 | 0 | null |
| 36 | quadratic | sparse_bucket | 1.018063 | 1.489 | 0 | 0 | 45 |
| 36 | quadratic | stream_plain | 0.586021 | 2.588 | 0 | 0 | 45 |
| 36 | quadratic | stream_exchange | 0.575917 | 2.633 | 0 | 0 | 45 |
| 36 | linear_drop | flat_cached | 3.433854 | 1.000 | 0 | 2 | null |
| 36 | linear_drop | sparse_bucket | 1.335230 | 2.572 | 0 | 2 | 84 |
| 36 | linear_drop | stream_plain | 0.801480 | 4.284 | 0 | 2 | 90 |
| 36 | linear_drop | stream_exchange | 0.771437 | 4.451 | 0 | 2 | 82 |
| 36 | restricted_cycle | flat_cached | 0.123813 | 1.000 | 7 | 4 | null |
| 36 | restricted_cycle | sparse_bucket | 0.757020 | 0.164 | 0 | 4 | 42 |
| 36 | restricted_cycle | stream_plain | 0.729688 | 0.170 | 0 | 4 | 46 |
| 36 | restricted_cycle | stream_exchange | 0.627000 | 0.197 | 0 | 4 | 45 |
| 36 | cross_cancel | flat_cached | 1.620354 | 1.000 | 0 | 8 | null |
| 36 | cross_cancel | sparse_bucket | 0.994000 | 1.630 | 0 | 8 | 44 |
| 36 | cross_cancel | stream_plain | 0.566229 | 2.862 | 0 | 8 | 44 |
| 36 | cross_cancel | stream_exchange | 0.551833 | 2.936 | 0 | 8 | 44 |

The active mask is recomputed from each input. Restricted-cycle fixtures use only eight embedded variables; cross-cancel fixtures guarantee an affine consequence from two nonlinear generators. Exact source dimensions, high rank and canonical tail are checked. Empty and nonempty tail workloads are both retained.

Schedule/layout retention is reported as counts, not byte estimates. Whole-worker RSS includes all arms and references. Packed XORs and sparse merge items are different diagnostics, not a common complete operation unit. Independent oracle preparation is outside arm timing and inside process receipts.

The controls retain pre-existing kernel strategies; no novelty claim attaches to sparse high-column elimination or linear-tail restriction. These standalone measurements do not execute the production solver, enumerate roots, recover scalars or measure full index calculus. Those costs remain null.
