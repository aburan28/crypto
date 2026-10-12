# Complete generic Boolean solve comparison

Evidence integrity: **PASS** across 72 cells and 5832 full-solve observations.

All results completed and independently verified: **True**. Full-solve gates against flat: **{'search': 'REJECTED', 'bucket': 'REJECTED', 'hybrid': 'REJECTED', 'small_flat': 'REJECTED', 'word_tail': 'REJECTED', 'merge_search': 'REJECTED', 'quadratic_state': 'PASS', 'packed_state': 'PASS'}**. Hybrid versus fastest control: **REJECTED**.

The table reports cold solve plus result-validation milliseconds over two holdout seeds and the declared rotated repetitions. A completion cost is null if any sample in the group is censored or lacks the required verification. Observed capped work is retained separately. Ratios use pooled medians; acceptance uses paired intervals.

| Variables | Family | Arm | Complete cost (ms) | Observed cost (ms) | Flat / arm | Median nodes | Median kernel calls | Median kernel fraction |
|---:|---|---|---:|---:|---:|---:|---:|---:|
| 12 | planted | search | 0.285292 | 0.285292 | 1.835 | 143 | 0 | 0.000 |
| 12 | planted | flat | 0.523521 | 0.523521 | 1.000 | 12 | 12 | 0.944 |
| 12 | planted | bucket | 5.622312 | 5.622312 | 0.093 | 12 | 12 | 0.994 |
| 12 | planted | hybrid | 0.510375 | 0.510375 | 1.026 | 12 | 12 | 0.935 |
| 12 | planted | small_flat | 0.429396 | 0.429396 | 1.219 | 86 | 72 | 0.540 |
| 12 | planted | word_tail | 0.415646 | 0.415646 | 1.260 | 86 | 72 | 0.543 |
| 12 | planted | merge_search | 0.268854 | 0.268854 | 1.947 | 143 | 0 | 0.000 |
| 12 | planted | quadratic_state | 0.107875 | 0.107875 | 4.853 | 143 | 0 | 0.000 |
| 12 | planted | packed_state | 0.070855 | 0.070855 | 7.389 | 143 | 0 | 0.000 |
| 12 | cross_planted | search | 0.354312 | 0.354312 | 1.054 | 171 | 0 | 0.000 |
| 12 | cross_planted | flat | 0.373271 | 0.373271 | 1.000 | 6 | 8 | 0.939 |
| 12 | cross_planted | bucket | 3.586479 | 3.586479 | 0.104 | 6 | 8 | 0.994 |
| 12 | cross_planted | hybrid | 0.344396 | 0.344396 | 1.084 | 6 | 8 | 0.938 |
| 12 | cross_planted | small_flat | 0.394208 | 0.394208 | 0.947 | 70 | 56 | 0.521 |
| 12 | cross_planted | word_tail | 0.318979 | 0.318979 | 1.170 | 70 | 56 | 0.500 |
| 12 | cross_planted | merge_search | 0.329520 | 0.329520 | 1.133 | 171 | 0 | 0.000 |
| 12 | cross_planted | quadratic_state | 0.128833 | 0.128833 | 2.897 | 171 | 0 | 0.000 |
| 12 | cross_planted | packed_state | 0.082396 | 0.082396 | 4.530 | 171 | 0 | 0.000 |
| 12 | unplanted | search | 0.426354 | 0.426354 | 1.331 | 210 | 0 | 0.000 |
| 12 | unplanted | flat | 0.567625 | 0.567625 | 1.000 | 12 | 12 | 0.945 |
| 12 | unplanted | bucket | 5.341563 | 5.341563 | 0.106 | 12 | 12 | 0.993 |
| 12 | unplanted | hybrid | 0.561875 | 0.561875 | 1.010 | 12 | 12 | 0.936 |
| 12 | unplanted | small_flat | 0.600958 | 0.600958 | 0.945 | 102 | 89 | 0.557 |
| 12 | unplanted | word_tail | 0.513000 | 0.513000 | 1.106 | 102 | 89 | 0.535 |
| 12 | unplanted | merge_search | 0.387562 | 0.387562 | 1.465 | 210 | 0 | 0.000 |
| 12 | unplanted | quadratic_state | 0.156562 | 0.156562 | 3.626 | 210 | 0 | 0.000 |
| 12 | unplanted | packed_state | 0.105042 | 0.105042 | 5.404 | 210 | 0 | 0.000 |
| 16 | planted | search | 0.608750 | 0.608750 | 3.476 | 220 | 0 | 0.000 |
| 16 | planted | flat | 2.115979 | 2.115979 | 1.000 | 28 | 28 | 0.957 |
| 16 | planted | bucket | 14.914813 | 14.914813 | 0.142 | 28 | 28 | 0.993 |
| 16 | planted | hybrid | 6.570688 | 6.570688 | 0.322 | 28 | 28 | 0.985 |
| 16 | planted | small_flat | 1.037896 | 1.037896 | 2.039 | 117 | 84 | 0.586 |
| 16 | planted | word_tail | 0.805625 | 0.805625 | 2.627 | 117 | 84 | 0.502 |
| 16 | planted | merge_search | 0.551480 | 0.551480 | 3.837 | 220 | 0 | 0.000 |
| 16 | planted | quadratic_state | 0.239709 | 0.239709 | 8.827 | 220 | 0 | 0.000 |
| 16 | planted | packed_state | 0.146334 | 0.146334 | 14.460 | 220 | 0 | 0.000 |
| 16 | cross_planted | search | 0.750104 | 0.750104 | 1.529 | 279 | 0 | 0.000 |
| 16 | cross_planted | flat | 1.147063 | 1.147063 | 1.000 | 12 | 13 | 0.956 |
| 16 | cross_planted | bucket | 8.596417 | 8.596417 | 0.133 | 12 | 13 | 0.994 |
| 16 | cross_planted | hybrid | 4.356125 | 4.356125 | 0.263 | 12 | 13 | 0.989 |
| 16 | cross_planted | small_flat | 0.692729 | 0.692729 | 1.656 | 86 | 66 | 0.563 |
| 16 | cross_planted | word_tail | 0.548063 | 0.548063 | 2.093 | 86 | 66 | 0.506 |
| 16 | cross_planted | merge_search | 0.672270 | 0.672270 | 1.706 | 279 | 0 | 0.000 |
| 16 | cross_planted | quadratic_state | 0.273083 | 0.273083 | 4.200 | 279 | 0 | 0.000 |
| 16 | cross_planted | packed_state | 0.171958 | 0.171958 | 6.671 | 279 | 0 | 0.000 |
| 16 | unplanted | search | 3.766542 | 3.766542 | 1.909 | 1336 | 0 | 0.000 |
| 16 | unplanted | flat | 7.192125 | 7.192125 | 1.000 | 123 | 123 | 0.941 |
| 16 | unplanted | bucket | 56.803229 | 56.803229 | 0.127 | 123 | 123 | 0.992 |
| 16 | unplanted | hybrid | 14.728021 | 14.728021 | 0.488 | 123 | 123 | 0.971 |
| 16 | unplanted | small_flat | 5.503833 | 5.503833 | 1.307 | 761 | 561 | 0.546 |
| 16 | unplanted | word_tail | 4.835792 | 4.835792 | 1.487 | 761 | 561 | 0.492 |
| 16 | unplanted | merge_search | 3.325229 | 3.325229 | 2.163 | 1336 | 0 | 0.000 |
| 16 | unplanted | quadratic_state | 1.359229 | 1.359229 | 5.291 | 1336 | 0 | 0.000 |
| 16 | unplanted | packed_state | 0.842333 | 0.842333 | 8.538 | 1336 | 0 | 0.000 |
| 20 | planted | search | 17.167167 | 17.167167 | 2.918 | 4851 | 0 | 0.000 |
| 20 | planted | flat | 50.087958 | 50.087958 | 1.000 | 516 | 549 | 0.953 |
| 20 | planted | bucket | 267.898374 | 267.898374 | 0.187 | 516 | 549 | 0.988 |
| 20 | planted | hybrid | 114.422250 | 114.422250 | 0.438 | 516 | 549 | 0.980 |
| 20 | planted | small_flat | 24.021083 | 24.021083 | 2.085 | 2759 | 1510 | 0.537 |
| 20 | planted | word_tail | 19.972145 | 19.972145 | 2.508 | 2759 | 1510 | 0.450 |
| 20 | planted | merge_search | 15.187145 | 15.187145 | 3.298 | 4851 | 0 | 0.000 |
| 20 | planted | quadratic_state | 6.101104 | 6.101104 | 8.210 | 4851 | 0 | 0.000 |
| 20 | planted | packed_state | 3.720813 | 3.720813 | 13.462 | 4851 | 0 | 0.000 |
| 20 | cross_planted | search | 15.505959 | 15.505959 | 1.117 | 4548 | 0 | 0.000 |
| 20 | cross_planted | flat | 17.312625 | 17.312625 | 1.000 | 219 | 221 | 0.956 |
| 20 | cross_planted | bucket | 94.879959 | 94.879959 | 0.182 | 219 | 221 | 0.989 |
| 20 | cross_planted | hybrid | 39.704499 | 39.704499 | 0.436 | 219 | 221 | 0.980 |
| 20 | cross_planted | small_flat | 17.605667 | 17.605667 | 0.983 | 1910 | 1232 | 0.523 |
| 20 | cross_planted | word_tail | 14.627334 | 14.627334 | 1.184 | 1910 | 1232 | 0.430 |
| 20 | cross_planted | merge_search | 13.600708 | 13.600708 | 1.273 | 4548 | 0 | 0.000 |
| 20 | cross_planted | quadratic_state | 5.501105 | 5.501105 | 3.147 | 4548 | 0 | 0.000 |
| 20 | cross_planted | packed_state | 3.374125 | 3.374125 | 5.131 | 4548 | 0 | 0.000 |
| 20 | unplanted | search | 23.459709 | 23.459709 | 2.300 | 6504 | 0 | 0.000 |
| 20 | unplanted | flat | 53.960604 | 53.960604 | 1.000 | 533 | 578 | 0.956 |
| 20 | unplanted | bucket | 261.285979 | 261.285979 | 0.207 | 533 | 578 | 0.988 |
| 20 | unplanted | hybrid | 109.221166 | 109.221166 | 0.494 | 533 | 578 | 0.974 |
| 20 | unplanted | small_flat | 30.318103 | 30.318103 | 1.780 | 3267 | 1849 | 0.555 |
| 20 | unplanted | word_tail | 25.147771 | 25.147771 | 2.146 | 3267 | 1849 | 0.469 |
| 20 | unplanted | merge_search | 20.647749 | 20.647749 | 2.613 | 6504 | 0 | 0.000 |
| 20 | unplanted | quadratic_state | 8.308708 | 8.308708 | 6.494 | 6504 | 0 | 0.000 |
| 20 | unplanted | packed_state | 5.152166 | 5.152166 | 10.473 | 6504 | 0 | 0.000 |
| 24 | planted | search | 53.764896 | 53.764896 | 4.708 | 13046 | 0 | 0.000 |
| 24 | planted | flat | 253.129771 | 253.129771 | 1.000 | 1420 | 1464 | 0.973 |
| 24 | planted | bucket | 1050.400459 | 1050.400459 | 0.241 | 1420 | 1464 | 0.991 |
| 24 | planted | hybrid | 447.711125 | 447.711125 | 0.565 | 1420 | 1464 | 0.982 |
| 24 | planted | small_flat | 61.827375 | 61.827375 | 4.094 | 10342 | 1264 | 0.240 |
| 24 | planted | word_tail | 56.283229 | 56.283229 | 4.497 | 10342 | 1264 | 0.166 |
| 24 | planted | merge_search | 47.039854 | 47.039854 | 5.381 | 13046 | 0 | 0.000 |
| 24 | planted | quadratic_state | 19.183375 | 19.183375 | 13.195 | 13046 | 0 | 0.000 |
| 24 | planted | packed_state | 11.966916 | 11.966916 | 21.152 | 13046 | 0 | 0.000 |
| 24 | cross_planted | search | 44.279021 | 44.279021 | 3.331 | 10720 | 0 | 0.000 |
| 24 | cross_planted | flat | 147.513229 | 147.513229 | 1.000 | 885 | 894 | 0.972 |
| 24 | cross_planted | bucket | 549.571083 | 549.571083 | 0.268 | 885 | 894 | 0.990 |
| 24 | cross_planted | hybrid | 253.588729 | 253.588729 | 0.582 | 885 | 894 | 0.980 |
| 24 | cross_planted | small_flat | 51.415292 | 51.415292 | 2.869 | 8831 | 913 | 0.240 |
| 24 | cross_planted | word_tail | 46.981979 | 46.981979 | 3.140 | 8831 | 913 | 0.165 |
| 24 | cross_planted | merge_search | 38.542645 | 38.542645 | 3.827 | 10720 | 0 | 0.000 |
| 24 | cross_planted | quadratic_state | 15.664146 | 15.664146 | 9.417 | 10720 | 0 | 0.000 |
| 24 | cross_planted | packed_state | 9.869708 | 9.869708 | 14.946 | 10720 | 0 | 0.000 |
| 24 | unplanted | search | 127.432626 | 127.432626 | 4.136 | 30179 | 0 | 0.000 |
| 24 | unplanted | flat | 527.080687 | 527.080687 | 1.000 | 3360 | 3462 | 0.969 |
| 24 | unplanted | bucket | 1767.499937 | 1767.499937 | 0.298 | 3360 | 3462 | 0.990 |
| 24 | unplanted | hybrid | 884.726334 | 884.726334 | 0.596 | 3360 | 3462 | 0.981 |
| 24 | unplanted | small_flat | 149.500770 | 149.500770 | 3.526 | 22558 | 3548 | 0.307 |
| 24 | unplanted | word_tail | 131.702937 | 131.702937 | 4.002 | 22558 | 3548 | 0.221 |
| 24 | unplanted | merge_search | 111.184270 | 111.184270 | 4.741 | 30179 | 0 | 0.000 |
| 24 | unplanted | quadratic_state | 46.031709 | 46.031709 | 11.450 | 30179 | 0 | 0.000 |
| 24 | unplanted | packed_state | 28.152125 | 28.152125 | 18.723 | 30179 | 0 | 0.000 |

Kernel-backed solvers have equal outcomes/models, logical counters and trace digests within each matched policy group. The search-only control may explore a different tree. SAT models are checked against original equations; UNSAT verification requires a completed independent search reference. UNKNOWN is censored, not UNSAT.

This driver solves bounded generated Boolean systems. It does not benchmark the repository inherited-F4 implementation, accept curve targets, recover scalars or establish index-calculus performance. Production and cryptanalytic costs stay null.

Selective one-word degree-2 full-solve gate: **REJECTED**. The small_flat and word_tail methods share a policy and must match exact traces/counters. They are compared separately with the always-degree-3 methods and search-only.

Ordered-specialization complete-solve gate: **REJECTED**. The search and merge_search arms have identical outcomes/models, logical work and trace digests. Only representation work in specialization changes. The gate compares merge_search with the fastest of all six prior full-solve methods.

Fixed-quadratic complete-solve gate: **REJECTED**. Compilation is charged to each cold solve. Search, merge_search and quadratic_state must match outcomes/models, logical counters and trace digests. The reference is the pointwise fastest of all seven prior methods.

| Variables | Family | Fastest prior / quadratic | 95% paired interval | Gate |
|---:|---|---:|---|---|
| 16 | planted | 2.298 | 2.271–2.308 | PASS |
| 16 | cross_planted | 1.659 | 1.069–2.410 | REJECTED |
| 16 | unplanted | 2.489 | 2.425–2.506 | PASS |
| 20 | planted | 2.454 | 2.432–2.476 | PASS |
| 20 | cross_planted | 2.191 | 2.175–2.436 | PASS |
| 20 | unplanted | 2.472 | 2.427–2.495 | PASS |
| 24 | planted | 2.397 | 2.313–2.453 | PASS |
| 24 | cross_planted | 2.417 | 2.309–2.467 | PASS |
| 24 | unplanted | 2.412 | 2.404–2.422 | PASS |

These paired intervals describe repeated timings on two fixed holdout fixtures per size/family, not a confidence interval over a population of Boolean systems. Fresh-worker RSS includes all arms and is not candidate-specific memory.

Packed-state current-frontier 2x gate: **REJECTED**. Incremental improvement gate (>1x lower bound): **PASS**. Cumulative historical-seven 2x gate: **REJECTED**.

Both fresh holdouts and the entire prior confirmation regression grid must pass. The current frontier includes the already faster quadratic_state method; the historical-seven comparison cannot substitute for it.

| Split | Variables | Family | Reference | Paired ratio | 95% interval | >2x gate |
|---|---:|---|---|---:|---|---|
| regression | 16 | planted | current_frontier | 1.641 | 1.586–1.665 | REJECTED |
| regression | 16 | planted | historical_seven | 3.638 | 3.158–4.068 | PASS |
| regression | 16 | cross_planted | current_frontier | 1.564 | 1.531–1.635 | REJECTED |
| regression | 16 | cross_planted | historical_seven | 3.486 | 3.264–3.582 | PASS |
| regression | 16 | unplanted | current_frontier | 1.573 | 1.539–1.596 | REJECTED |
| regression | 16 | unplanted | historical_seven | 3.963 | 3.842–4.016 | PASS |
| regression | 20 | planted | current_frontier | 1.644 | 1.623–1.659 | REJECTED |
| regression | 20 | planted | historical_seven | 4.009 | 3.931–4.048 | PASS |
| regression | 20 | cross_planted | current_frontier | 1.621 | 1.601–1.635 | REJECTED |
| regression | 20 | cross_planted | historical_seven | 3.909 | 3.866–4.009 | PASS |
| regression | 20 | unplanted | current_frontier | 1.608 | 1.584–1.626 | REJECTED |
| regression | 20 | unplanted | historical_seven | 3.971 | 3.932–4.022 | PASS |
| regression | 24 | planted | current_frontier | 1.648 | 1.636–1.663 | REJECTED |
| regression | 24 | planted | historical_seven | 3.941 | 3.929–3.951 | PASS |
| regression | 24 | cross_planted | current_frontier | 1.626 | 1.607–1.679 | REJECTED |
| regression | 24 | cross_planted | historical_seven | 3.948 | 3.920–3.974 | PASS |
| regression | 24 | unplanted | current_frontier | 1.670 | 1.659–1.678 | REJECTED |
| regression | 24 | unplanted | historical_seven | 3.990 | 3.972–4.005 | PASS |
| holdout | 16 | planted | current_frontier | 1.639 | 1.603–1.665 | REJECTED |
| holdout | 16 | planted | historical_seven | 3.741 | 3.680–3.812 | PASS |
| holdout | 16 | cross_planted | current_frontier | 1.653 | 1.621–1.694 | REJECTED |
| holdout | 16 | cross_planted | historical_seven | 2.632 | 1.792–3.957 | REJECTED |
| holdout | 16 | unplanted | current_frontier | 1.620 | 1.591–1.642 | REJECTED |
| holdout | 16 | unplanted | historical_seven | 4.021 | 3.886–4.088 | PASS |
| holdout | 20 | planted | current_frontier | 1.626 | 1.579–1.653 | REJECTED |
| holdout | 20 | planted | historical_seven | 3.942 | 3.869–4.056 | PASS |
| holdout | 20 | cross_planted | current_frontier | 1.631 | 1.618–1.657 | REJECTED |
| holdout | 20 | cross_planted | historical_seven | 3.555 | 3.511–4.027 | PASS |
| holdout | 20 | unplanted | current_frontier | 1.609 | 1.592–1.629 | REJECTED |
| holdout | 20 | unplanted | historical_seven | 3.950 | 3.885–4.056 | PASS |
| holdout | 24 | planted | current_frontier | 1.630 | 1.607–1.689 | REJECTED |
| holdout | 24 | planted | historical_seven | 3.938 | 3.909–3.954 | PASS |
| holdout | 24 | cross_planted | current_frontier | 1.597 | 1.585–1.621 | REJECTED |
| holdout | 24 | cross_planted | historical_seven | 3.911 | 3.855–3.938 | PASS |
| holdout | 24 | unplanted | current_frontier | 1.630 | 1.616–1.651 | REJECTED |
| holdout | 24 | unplanted | historical_seven | 3.918 | 3.895–3.983 | PASS |
