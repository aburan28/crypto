# Complete generic Boolean solve comparison

Evidence integrity: **PASS** across 48 cells and 2352 full-solve observations.

All results completed and independently verified: **True**. Full-solve gates against flat: **{'search': 'REJECTED', 'bucket': 'REJECTED', 'hybrid': 'REJECTED', 'small_flat': 'REJECTED', 'word_tail': 'REJECTED', 'merge_search': 'REJECTED'}**. Hybrid versus fastest control: **REJECTED**.

The table reports cold solve plus result-validation milliseconds over two holdout seeds and the declared rotated repetitions. A completion cost is null if any sample in the group is censored or lacks the required verification. Observed capped work is retained separately. Ratios use pooled medians; acceptance uses paired intervals.

| Variables | Family | Arm | Complete cost (ms) | Observed cost (ms) | Flat / arm | Median nodes | Median kernel calls | Median kernel fraction |
|---:|---|---|---:|---:|---:|---:|---:|---:|
| 12 | planted | search | 0.317875 | 0.317875 | 1.789 | 152 | 0 | 0.000 |
| 12 | planted | flat | 0.568771 | 0.568771 | 1.000 | 11 | 11 | 0.951 |
| 12 | planted | bucket | 5.964001 | 5.964001 | 0.095 | 11 | 11 | 0.994 |
| 12 | planted | hybrid | 0.536063 | 0.536063 | 1.061 | 11 | 11 | 0.945 |
| 12 | planted | small_flat | 0.499896 | 0.499896 | 1.138 | 82 | 72 | 0.585 |
| 12 | planted | word_tail | 0.419291 | 0.419291 | 1.357 | 82 | 72 | 0.546 |
| 12 | planted | merge_search | 0.291521 | 0.291521 | 1.951 | 152 | 0 | 0.000 |
| 12 | cross_planted | search | 0.297479 | 0.297479 | 1.241 | 144 | 0 | 0.000 |
| 12 | cross_planted | flat | 0.369188 | 0.369188 | 1.000 | 6 | 7 | 0.943 |
| 12 | cross_planted | bucket | 3.449500 | 3.449500 | 0.107 | 6 | 7 | 0.993 |
| 12 | cross_planted | hybrid | 0.360895 | 0.360895 | 1.023 | 6 | 7 | 0.942 |
| 12 | cross_planted | small_flat | 0.306875 | 0.306875 | 1.203 | 52 | 43 | 0.528 |
| 12 | cross_planted | word_tail | 0.239250 | 0.239250 | 1.543 | 52 | 43 | 0.519 |
| 12 | cross_planted | merge_search | 0.276209 | 0.276209 | 1.337 | 144 | 0 | 0.000 |
| 12 | unplanted | search | 0.671958 | 0.671958 | 1.011 | 353 | 0 | 0.000 |
| 12 | unplanted | flat | 0.679604 | 0.679604 | 1.000 | 15 | 15 | 0.942 |
| 12 | unplanted | bucket | 6.547292 | 6.547292 | 0.104 | 15 | 15 | 0.993 |
| 12 | unplanted | hybrid | 0.640459 | 0.640459 | 1.061 | 15 | 15 | 0.935 |
| 12 | unplanted | small_flat | 0.908062 | 0.908062 | 0.748 | 167 | 143 | 0.511 |
| 12 | unplanted | word_tail | 0.869312 | 0.869312 | 0.782 | 167 | 143 | 0.520 |
| 12 | unplanted | merge_search | 0.639333 | 0.639333 | 1.063 | 353 | 0 | 0.000 |
| 16 | planted | search | 0.961875 | 0.961875 | 3.601 | 342 | 0 | 0.000 |
| 16 | planted | flat | 3.463896 | 3.463896 | 1.000 | 52 | 52 | 0.950 |
| 16 | planted | bucket | 26.954313 | 26.954313 | 0.129 | 52 | 52 | 0.994 |
| 16 | planted | hybrid | 8.979271 | 8.979271 | 0.386 | 52 | 52 | 0.983 |
| 16 | planted | small_flat | 1.574666 | 1.574666 | 2.200 | 209 | 142 | 0.564 |
| 16 | planted | word_tail | 1.309437 | 1.309437 | 2.645 | 209 | 142 | 0.479 |
| 16 | planted | merge_search | 0.864667 | 0.864667 | 4.006 | 342 | 0 | 0.000 |
| 16 | cross_planted | search | 1.344542 | 1.344542 | 1.895 | 432 | 0 | 0.000 |
| 16 | cross_planted | flat | 2.548500 | 2.548500 | 1.000 | 30 | 31 | 0.948 |
| 16 | cross_planted | bucket | 17.600104 | 17.600104 | 0.145 | 30 | 31 | 0.993 |
| 16 | cross_planted | hybrid | 6.427041 | 6.427041 | 0.397 | 30 | 31 | 0.982 |
| 16 | cross_planted | small_flat | 1.397146 | 1.397146 | 1.824 | 179 | 138 | 0.533 |
| 16 | cross_planted | word_tail | 1.187542 | 1.187542 | 2.146 | 179 | 138 | 0.476 |
| 16 | cross_planted | merge_search | 1.075271 | 1.075271 | 2.370 | 432 | 0 | 0.000 |
| 16 | unplanted | search | 3.165896 | 3.165896 | 2.298 | 1120 | 0 | 0.000 |
| 16 | unplanted | flat | 7.274730 | 7.274730 | 1.000 | 126 | 126 | 0.939 |
| 16 | unplanted | bucket | 53.470292 | 53.470292 | 0.136 | 126 | 126 | 0.990 |
| 16 | unplanted | hybrid | 14.086937 | 14.086937 | 0.516 | 126 | 126 | 0.967 |
| 16 | unplanted | small_flat | 4.496812 | 4.496812 | 1.618 | 580 | 435 | 0.558 |
| 16 | unplanted | word_tail | 3.896959 | 3.896959 | 1.867 | 580 | 435 | 0.501 |
| 16 | unplanted | merge_search | 2.888146 | 2.888146 | 2.519 | 1120 | 0 | 0.000 |
| 20 | planted | search | 13.569542 | 13.569542 | 1.894 | 4014 | 0 | 0.000 |
| 20 | planted | flat | 25.698166 | 25.698166 | 1.000 | 304 | 316 | 0.949 |
| 20 | planted | bucket | 162.206979 | 162.206979 | 0.158 | 304 | 316 | 0.991 |
| 20 | planted | hybrid | 69.864937 | 69.864937 | 0.368 | 304 | 316 | 0.980 |
| 20 | planted | small_flat | 19.200751 | 19.200751 | 1.338 | 2192 | 1299 | 0.546 |
| 20 | planted | word_tail | 16.101937 | 16.101937 | 1.596 | 2192 | 1299 | 0.471 |
| 20 | planted | merge_search | 12.142021 | 12.142021 | 2.116 | 4014 | 0 | 0.000 |
| 20 | cross_planted | search | 11.906125 | 11.906125 | 1.833 | 3626 | 0 | 0.000 |
| 20 | cross_planted | flat | 21.819626 | 21.819626 | 1.000 | 306 | 308 | 0.949 |
| 20 | cross_planted | bucket | 105.830438 | 105.830438 | 0.206 | 306 | 308 | 0.989 |
| 20 | cross_planted | hybrid | 40.167667 | 40.167667 | 0.543 | 306 | 308 | 0.973 |
| 20 | cross_planted | small_flat | 14.346271 | 14.346271 | 1.521 | 1616 | 1001 | 0.551 |
| 20 | cross_planted | word_tail | 11.870708 | 11.870708 | 1.838 | 1616 | 1001 | 0.473 |
| 20 | cross_planted | merge_search | 10.589313 | 10.589313 | 2.061 | 3626 | 0 | 0.000 |
| 20 | unplanted | search | 8.174479 | 8.174479 | 3.275 | 2358 | 0 | 0.000 |
| 20 | unplanted | flat | 26.767750 | 26.767750 | 1.000 | 282 | 304 | 0.956 |
| 20 | unplanted | bucket | 144.901979 | 144.901979 | 0.185 | 282 | 304 | 0.991 |
| 20 | unplanted | hybrid | 59.572313 | 59.572313 | 0.449 | 282 | 304 | 0.979 |
| 20 | unplanted | small_flat | 11.551500 | 11.551500 | 2.317 | 1348 | 727 | 0.533 |
| 20 | unplanted | word_tail | 9.486937 | 9.486937 | 2.822 | 1348 | 727 | 0.446 |
| 20 | unplanted | merge_search | 7.203854 | 7.203854 | 3.716 | 2358 | 0 | 0.000 |
| 24 | planted | search | 98.194229 | 98.194229 | 3.033 | 24182 | 0 | 0.000 |
| 24 | planted | flat | 297.822917 | 297.822917 | 1.000 | 2056 | 2128 | 0.968 |
| 24 | planted | bucket | 1068.500479 | 1068.500479 | 0.279 | 2056 | 2128 | 0.990 |
| 24 | planted | hybrid | 429.624542 | 429.624542 | 0.693 | 2056 | 2128 | 0.978 |
| 24 | planted | small_flat | 114.024709 | 114.024709 | 2.612 | 14553 | 4456 | 0.371 |
| 24 | planted | word_tail | 99.382437 | 99.382437 | 2.997 | 14553 | 4456 | 0.292 |
| 24 | planted | merge_search | 87.016291 | 87.016291 | 3.423 | 24182 | 0 | 0.000 |
| 24 | cross_planted | search | 86.760063 | 86.760063 | 1.441 | 22274 | 0 | 0.000 |
| 24 | cross_planted | flat | 125.040250 | 125.040250 | 1.000 | 891 | 932 | 0.969 |
| 24 | cross_planted | bucket | 466.614896 | 466.614896 | 0.268 | 891 | 932 | 0.991 |
| 24 | cross_planted | hybrid | 176.332562 | 176.332562 | 0.709 | 891 | 932 | 0.980 |
| 24 | cross_planted | small_flat | 99.104333 | 99.104333 | 1.262 | 13376 | 3837 | 0.372 |
| 24 | cross_planted | word_tail | 85.962625 | 85.962625 | 1.455 | 13376 | 3837 | 0.294 |
| 24 | cross_planted | merge_search | 76.643687 | 76.643687 | 1.631 | 22274 | 0 | 0.000 |
| 24 | unplanted | search | 146.671084 | 146.671084 | 4.071 | 35097 | 0 | 0.000 |
| 24 | unplanted | flat | 597.041521 | 597.041521 | 1.000 | 3935 | 3955 | 0.969 |
| 24 | unplanted | bucket | 2137.548020 | 2137.548020 | 0.279 | 3935 | 3955 | 0.990 |
| 24 | unplanted | hybrid | 928.377542 | 928.377542 | 0.643 | 3935 | 3955 | 0.979 |
| 24 | unplanted | small_flat | 171.460291 | 171.460291 | 3.482 | 24762 | 4804 | 0.343 |
| 24 | unplanted | word_tail | 151.653458 | 151.653458 | 3.937 | 24762 | 4804 | 0.267 |
| 24 | unplanted | merge_search | 128.637125 | 128.637125 | 4.641 | 35097 | 0 | 0.000 |

Kernel-backed solvers have equal outcomes/models, logical counters and trace digests within each matched policy group. The search-only control may explore a different tree. SAT models are checked against original equations; UNSAT verification requires a completed independent search reference. UNKNOWN is censored, not UNSAT.

This driver solves bounded generated Boolean systems. It does not benchmark the repository inherited-F4 implementation, accept curve targets, recover scalars or establish index-calculus performance. Production and cryptanalytic costs stay null.

Selective one-word degree-2 full-solve gate: **REJECTED**. The small_flat and word_tail methods share a policy and must match exact traces/counters. They are compared separately with the always-degree-3 methods and search-only.

Ordered-specialization complete-solve gate: **REJECTED**. The search and merge_search arms have identical outcomes/models, logical work and trace digests. Only representation work in specialization changes. The gate compares merge_search with the fastest of all six prior full-solve methods.
