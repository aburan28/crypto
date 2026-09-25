# Complete generic Boolean solve comparison

Evidence integrity: **PASS** across 48 cells and 3072 full-solve observations.

All results completed and independently verified: **True**. Full-solve gates against flat: **{'search': 'REJECTED', 'bucket': 'REJECTED', 'hybrid': 'REJECTED', 'small_flat': 'REJECTED', 'word_tail': 'REJECTED', 'merge_search': 'REJECTED', 'quadratic_state': 'PASS'}**. Hybrid versus fastest control: **REJECTED**.

The table reports cold solve plus result-validation milliseconds over two holdout seeds and the declared rotated repetitions. A completion cost is null if any sample in the group is censored or lacks the required verification. Observed capped work is retained separately. Ratios use pooled medians; acceptance uses paired intervals.

| Variables | Family | Arm | Complete cost (ms) | Observed cost (ms) | Flat / arm | Median nodes | Median kernel calls | Median kernel fraction |
|---:|---|---|---:|---:|---:|---:|---:|---:|
| 12 | planted | search | 0.334417 | 0.334417 | 1.598 | 158 | 0 | 0.000 |
| 12 | planted | flat | 0.534500 | 0.534500 | 1.000 | 12 | 12 | 0.941 |
| 12 | planted | bucket | 5.493083 | 5.493083 | 0.097 | 12 | 12 | 0.993 |
| 12 | planted | hybrid | 0.520396 | 0.520396 | 1.027 | 12 | 12 | 0.939 |
| 12 | planted | small_flat | 0.518938 | 0.518938 | 1.030 | 86 | 73 | 0.550 |
| 12 | planted | word_tail | 0.446458 | 0.446458 | 1.197 | 86 | 73 | 0.529 |
| 12 | planted | merge_search | 0.317167 | 0.317167 | 1.685 | 158 | 0 | 0.000 |
| 12 | planted | quadratic_state | 0.158209 | 0.158209 | 3.378 | 158 | 0 | 0.000 |
| 12 | cross_planted | search | 0.344021 | 0.344021 | 1.016 | 168 | 0 | 0.000 |
| 12 | cross_planted | flat | 0.349542 | 0.349542 | 1.000 | 6 | 7 | 0.945 |
| 12 | cross_planted | bucket | 3.016938 | 3.016938 | 0.116 | 6 | 7 | 0.994 |
| 12 | cross_planted | hybrid | 0.332062 | 0.332062 | 1.053 | 6 | 7 | 0.944 |
| 12 | cross_planted | small_flat | 0.310500 | 0.310500 | 1.126 | 57 | 47 | 0.554 |
| 12 | cross_planted | word_tail | 0.263041 | 0.263041 | 1.329 | 57 | 47 | 0.537 |
| 12 | cross_planted | merge_search | 0.321229 | 0.321229 | 1.088 | 168 | 0 | 0.000 |
| 12 | cross_planted | quadratic_state | 0.124917 | 0.124917 | 2.798 | 168 | 0 | 0.000 |
| 12 | unplanted | search | 0.454729 | 0.454729 | 1.557 | 206 | 0 | 0.000 |
| 12 | unplanted | flat | 0.708187 | 0.708187 | 1.000 | 15 | 15 | 0.946 |
| 12 | unplanted | bucket | 7.221292 | 7.221292 | 0.098 | 15 | 15 | 0.993 |
| 12 | unplanted | hybrid | 0.708562 | 0.708562 | 0.999 | 15 | 15 | 0.937 |
| 12 | unplanted | small_flat | 0.737604 | 0.737604 | 0.960 | 127 | 94 | 0.519 |
| 12 | unplanted | word_tail | 0.611437 | 0.611437 | 1.158 | 127 | 94 | 0.494 |
| 12 | unplanted | merge_search | 0.430312 | 0.430312 | 1.646 | 206 | 0 | 0.000 |
| 12 | unplanted | quadratic_state | 0.181417 | 0.181417 | 3.904 | 206 | 0 | 0.000 |
| 16 | planted | search | 0.897188 | 0.897188 | 2.558 | 322 | 0 | 0.000 |
| 16 | planted | flat | 2.295396 | 2.295396 | 1.000 | 34 | 34 | 0.955 |
| 16 | planted | bucket | 21.711063 | 21.711063 | 0.106 | 34 | 34 | 0.994 |
| 16 | planted | hybrid | 7.475604 | 7.475604 | 0.307 | 34 | 34 | 0.985 |
| 16 | planted | small_flat | 1.237208 | 1.237208 | 1.855 | 152 | 124 | 0.543 |
| 16 | planted | word_tail | 1.020417 | 1.020417 | 2.249 | 152 | 124 | 0.510 |
| 16 | planted | merge_search | 0.784500 | 0.784500 | 2.926 | 322 | 0 | 0.000 |
| 16 | planted | quadratic_state | 0.298354 | 0.298354 | 7.694 | 322 | 0 | 0.000 |
| 16 | cross_planted | search | 0.944416 | 0.944416 | 1.568 | 386 | 0 | 0.000 |
| 16 | cross_planted | flat | 1.481250 | 1.481250 | 1.000 | 20 | 21 | 0.954 |
| 16 | cross_planted | bucket | 12.334146 | 12.334146 | 0.120 | 20 | 21 | 0.994 |
| 16 | cross_planted | hybrid | 5.498583 | 5.498583 | 0.269 | 20 | 21 | 0.987 |
| 16 | cross_planted | small_flat | 0.907979 | 0.907979 | 1.631 | 118 | 100 | 0.547 |
| 16 | cross_planted | word_tail | 0.749958 | 0.749958 | 1.975 | 118 | 100 | 0.517 |
| 16 | cross_planted | merge_search | 0.847271 | 0.847271 | 1.748 | 386 | 0 | 0.000 |
| 16 | cross_planted | quadratic_state | 0.346646 | 0.346646 | 4.273 | 386 | 0 | 0.000 |
| 16 | unplanted | search | 3.929458 | 3.929458 | 1.875 | 1409 | 0 | 0.000 |
| 16 | unplanted | flat | 7.368563 | 7.368563 | 1.000 | 127 | 127 | 0.938 |
| 16 | unplanted | bucket | 56.798271 | 56.798271 | 0.130 | 127 | 127 | 0.991 |
| 16 | unplanted | hybrid | 13.249146 | 13.249146 | 0.556 | 127 | 127 | 0.965 |
| 16 | unplanted | small_flat | 5.709708 | 5.709708 | 1.291 | 773 | 584 | 0.555 |
| 16 | unplanted | word_tail | 4.951625 | 4.951625 | 1.488 | 773 | 584 | 0.495 |
| 16 | unplanted | merge_search | 3.498041 | 3.498041 | 2.106 | 1409 | 0 | 0.000 |
| 16 | unplanted | quadratic_state | 1.401750 | 1.401750 | 5.257 | 1409 | 0 | 0.000 |
| 20 | planted | search | 13.752230 | 13.752230 | 1.622 | 4232 | 0 | 0.000 |
| 20 | planted | flat | 22.312188 | 22.312188 | 1.000 | 278 | 299 | 0.948 |
| 20 | planted | bucket | 132.581063 | 132.581063 | 0.168 | 278 | 299 | 0.990 |
| 20 | planted | hybrid | 56.793021 | 56.793021 | 0.393 | 278 | 299 | 0.978 |
| 20 | planted | small_flat | 19.062020 | 19.062020 | 1.171 | 2196 | 1351 | 0.525 |
| 20 | planted | word_tail | 15.604791 | 15.604791 | 1.430 | 2196 | 1351 | 0.457 |
| 20 | planted | merge_search | 12.127271 | 12.127271 | 1.840 | 4232 | 0 | 0.000 |
| 20 | planted | quadratic_state | 5.012083 | 5.012083 | 4.452 | 4232 | 0 | 0.000 |
| 20 | cross_planted | search | 15.339063 | 15.339063 | 1.602 | 4728 | 0 | 0.000 |
| 20 | cross_planted | flat | 24.565667 | 24.565667 | 1.000 | 322 | 327 | 0.950 |
| 20 | cross_planted | bucket | 129.697354 | 129.697354 | 0.189 | 322 | 327 | 0.990 |
| 20 | cross_planted | hybrid | 45.366667 | 45.366667 | 0.541 | 322 | 327 | 0.972 |
| 20 | cross_planted | small_flat | 17.356375 | 17.356375 | 1.415 | 1736 | 1136 | 0.547 |
| 20 | cross_planted | word_tail | 13.731375 | 13.731375 | 1.789 | 1736 | 1136 | 0.465 |
| 20 | cross_planted | merge_search | 13.652479 | 13.652479 | 1.799 | 4728 | 0 | 0.000 |
| 20 | cross_planted | quadratic_state | 5.509021 | 5.509021 | 4.459 | 4728 | 0 | 0.000 |
| 20 | unplanted | search | 16.625062 | 16.625062 | 2.177 | 4716 | 0 | 0.000 |
| 20 | unplanted | flat | 36.199792 | 36.199792 | 1.000 | 332 | 348 | 0.955 |
| 20 | unplanted | bucket | 212.895895 | 212.895895 | 0.170 | 332 | 348 | 0.990 |
| 20 | unplanted | hybrid | 76.601104 | 76.601104 | 0.473 | 332 | 348 | 0.971 |
| 20 | unplanted | small_flat | 21.998500 | 21.998500 | 1.646 | 2264 | 1420 | 0.579 |
| 20 | unplanted | word_tail | 18.003084 | 18.003084 | 2.011 | 2264 | 1420 | 0.504 |
| 20 | unplanted | merge_search | 14.555333 | 14.555333 | 2.487 | 4716 | 0 | 0.000 |
| 20 | unplanted | quadratic_state | 5.958083 | 5.958083 | 6.076 | 4716 | 0 | 0.000 |
| 24 | planted | search | 83.948353 | 83.948353 | 3.866 | 20110 | 0 | 0.000 |
| 24 | planted | flat | 324.505584 | 324.505584 | 1.000 | 2204 | 2235 | 0.968 |
| 24 | planted | bucket | 1039.781000 | 1039.781000 | 0.312 | 2204 | 2235 | 0.989 |
| 24 | planted | hybrid | 432.203792 | 432.203792 | 0.751 | 2204 | 2235 | 0.976 |
| 24 | planted | small_flat | 98.751625 | 98.751625 | 3.286 | 13684 | 2952 | 0.370 |
| 24 | planted | word_tail | 86.087291 | 86.087291 | 3.769 | 13684 | 2952 | 0.285 |
| 24 | planted | merge_search | 73.976250 | 73.976250 | 4.387 | 20110 | 0 | 0.000 |
| 24 | planted | quadratic_state | 30.557917 | 30.557917 | 10.619 | 20110 | 0 | 0.000 |
| 24 | cross_planted | search | 50.172584 | 50.172584 | 2.200 | 12320 | 0 | 0.000 |
| 24 | cross_planted | flat | 110.399500 | 110.399500 | 1.000 | 792 | 806 | 0.966 |
| 24 | cross_planted | bucket | 340.670584 | 340.670584 | 0.324 | 792 | 806 | 0.989 |
| 24 | cross_planted | hybrid | 156.305584 | 156.305584 | 0.706 | 792 | 806 | 0.977 |
| 24 | cross_planted | small_flat | 57.306479 | 57.306479 | 1.926 | 8448 | 1683 | 0.348 |
| 24 | cross_planted | word_tail | 49.644729 | 49.644729 | 2.224 | 8448 | 1683 | 0.270 |
| 24 | cross_planted | merge_search | 43.469708 | 43.469708 | 2.540 | 12320 | 0 | 0.000 |
| 24 | cross_planted | quadratic_state | 17.770312 | 17.770312 | 6.213 | 12320 | 0 | 0.000 |
| 24 | unplanted | search | 172.670625 | 172.670625 | 3.256 | 41017 | 0 | 0.000 |
| 24 | unplanted | flat | 562.174501 | 562.174501 | 1.000 | 3908 | 3950 | 0.967 |
| 24 | unplanted | bucket | 2074.648312 | 2074.648312 | 0.271 | 3908 | 3950 | 0.990 |
| 24 | unplanted | hybrid | 920.137104 | 920.137104 | 0.611 | 3908 | 3950 | 0.979 |
| 24 | unplanted | small_flat | 201.832708 | 201.832708 | 2.785 | 26658 | 6652 | 0.400 |
| 24 | unplanted | word_tail | 175.701604 | 175.701604 | 3.200 | 26658 | 6652 | 0.315 |
| 24 | unplanted | merge_search | 151.862917 | 151.862917 | 3.702 | 41017 | 0 | 0.000 |
| 24 | unplanted | quadratic_state | 62.343125 | 62.343125 | 9.017 | 41017 | 0 | 0.000 |

Kernel-backed solvers have equal outcomes/models, logical counters and trace digests within each matched policy group. The search-only control may explore a different tree. SAT models are checked against original equations; UNSAT verification requires a completed independent search reference. UNKNOWN is censored, not UNSAT.

This driver solves bounded generated Boolean systems. It does not benchmark the repository inherited-F4 implementation, accept curve targets, recover scalars or establish index-calculus performance. Production and cryptanalytic costs stay null.

Selective one-word degree-2 full-solve gate: **REJECTED**. The small_flat and word_tail methods share a policy and must match exact traces/counters. They are compared separately with the always-degree-3 methods and search-only.

Ordered-specialization complete-solve gate: **REJECTED**. The search and merge_search arms have identical outcomes/models, logical work and trace digests. Only representation work in specialization changes. The gate compares merge_search with the fastest of all six prior full-solve methods.

Fixed-quadratic complete-solve gate: **REJECTED**. Compilation is charged to each cold solve. Search, merge_search and quadratic_state must match outcomes/models, logical counters and trace digests. The reference is the pointwise fastest of all seven prior methods.

| Variables | Family | Fastest prior / quadratic | 95% paired interval | Gate |
|---:|---|---:|---|---|
| 16 | planted | 2.439 | 1.950–2.609 | REJECTED |
| 16 | cross_planted | 2.102 | 1.994–2.237 | REJECTED |
| 16 | unplanted | 2.553 | 2.517–2.611 | PASS |
| 20 | planted | 2.435 | 2.423–2.479 | PASS |
| 20 | cross_planted | 2.441 | 2.373–2.501 | PASS |
| 20 | unplanted | 2.450 | 2.417–2.523 | PASS |
| 24 | planted | 2.409 | 2.399–2.426 | PASS |
| 24 | cross_planted | 2.448 | 2.392–2.464 | PASS |
| 24 | unplanted | 2.430 | 2.412–2.451 | PASS |

These paired intervals describe repeated timings on two fixed holdout fixtures per size/family, not a confidence interval over a population of Boolean systems. Fresh-worker RSS includes all arms and is not candidate-specific memory.
