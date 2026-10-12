# Complete generic Boolean solve comparison

Evidence integrity: **PASS** across 48 cells and 768 full-solve observations.

All results completed and independently verified: **True**. Full-solve gates against flat: **{'search': 'REJECTED', 'bucket': 'REJECTED', 'hybrid': 'REJECTED'}**. Hybrid versus fastest control: **REJECTED**.

The table reports cold solve plus result-validation milliseconds over two holdout seeds and four rotated repetitions. A completion cost is null if any sample in the group is censored or lacks the required verification. Observed capped work is retained separately. Ratios use pooled medians; acceptance uses paired intervals.

| Variables | Family | Arm | Complete cost (ms) | Observed cost (ms) | Flat / arm | Median nodes | Median kernel calls | Median kernel fraction |
|---:|---|---|---:|---:|---:|---:|---:|---:|
| 12 | planted | search | 0.276104 | 0.276104 | 1.776 | 126 | 0 | 0.000 |
| 12 | planted | flat | 0.490333 | 0.490333 | 1.000 | 12 | 12 | 0.945 |
| 12 | planted | bucket | 5.540146 | 5.540146 | 0.089 | 12 | 12 | 0.993 |
| 12 | planted | hybrid | 0.508729 | 0.508729 | 0.964 | 12 | 12 | 0.943 |
| 12 | cross_planted | search | 0.298458 | 0.298458 | 1.114 | 148 | 0 | 0.000 |
| 12 | cross_planted | flat | 0.332437 | 0.332437 | 1.000 | 6 | 8 | 0.942 |
| 12 | cross_planted | bucket | 3.093520 | 3.093520 | 0.107 | 6 | 8 | 0.993 |
| 12 | cross_planted | hybrid | 0.356729 | 0.356729 | 0.932 | 6 | 8 | 0.944 |
| 12 | unplanted | search | 0.524021 | 0.524021 | 1.113 | 244 | 0 | 0.000 |
| 12 | unplanted | flat | 0.583146 | 0.583146 | 1.000 | 15 | 15 | 0.935 |
| 12 | unplanted | bucket | 6.870563 | 6.870563 | 0.085 | 15 | 15 | 0.993 |
| 12 | unplanted | hybrid | 0.664687 | 0.664687 | 0.877 | 15 | 15 | 0.933 |
| 16 | planted | search | 2.376354 | 2.376354 | 1.836 | 830 | 0 | 0.000 |
| 16 | planted | flat | 4.361854 | 4.361854 | 1.000 | 65 | 65 | 0.948 |
| 16 | planted | bucket | 41.199959 | 41.199959 | 0.106 | 65 | 65 | 0.994 |
| 16 | planted | hybrid | 10.660562 | 10.660562 | 0.409 | 65 | 65 | 0.979 |
| 16 | cross_planted | search | 1.633584 | 1.633584 | 1.233 | 600 | 0 | 0.000 |
| 16 | cross_planted | flat | 2.013645 | 2.013645 | 1.000 | 32 | 32 | 0.952 |
| 16 | cross_planted | bucket | 18.552396 | 18.552396 | 0.109 | 32 | 32 | 0.994 |
| 16 | cross_planted | hybrid | 6.325062 | 6.325062 | 0.318 | 32 | 32 | 0.983 |
| 16 | unplanted | search | 4.096292 | 4.096292 | 1.778 | 1461 | 0 | 0.000 |
| 16 | unplanted | flat | 7.283542 | 7.283542 | 1.000 | 127 | 127 | 0.938 |
| 16 | unplanted | bucket | 71.328230 | 71.328230 | 0.102 | 127 | 127 | 0.993 |
| 16 | unplanted | hybrid | 15.603438 | 15.603438 | 0.467 | 127 | 127 | 0.971 |
| 20 | planted | search | 5.588249 | 5.588249 | 2.676 | 1726 | 0 | 0.000 |
| 20 | planted | flat | 14.953062 | 14.953062 | 1.000 | 148 | 157 | 0.958 |
| 20 | planted | bucket | 87.252583 | 87.252583 | 0.171 | 148 | 157 | 0.992 |
| 20 | planted | hybrid | 40.627104 | 40.627104 | 0.368 | 148 | 157 | 0.985 |
| 20 | cross_planted | search | 4.583417 | 4.583417 | 1.829 | 1398 | 0 | 0.000 |
| 20 | cross_planted | flat | 8.381854 | 8.381854 | 1.000 | 93 | 94 | 0.959 |
| 20 | cross_planted | bucket | 53.696562 | 53.696562 | 0.156 | 93 | 94 | 0.993 |
| 20 | cross_planted | hybrid | 30.683499 | 30.683499 | 0.273 | 93 | 94 | 0.988 |
| 20 | unplanted | search | 23.124021 | 23.124021 | 2.750 | 6714 | 0 | 0.000 |
| 20 | unplanted | flat | 63.589917 | 63.589917 | 1.000 | 550 | 566 | 0.963 |
| 20 | unplanted | bucket | 393.535542 | 393.535542 | 0.162 | 550 | 566 | 0.993 |
| 20 | unplanted | hybrid | 151.691521 | 151.691521 | 0.419 | 550 | 566 | 0.983 |
| 24 | planted | search | 67.992729 | 67.992729 | 4.484 | 16314 | 0 | 0.000 |
| 24 | planted | flat | 304.904875 | 304.904875 | 1.000 | 1896 | 1944 | 0.970 |
| 24 | planted | bucket | 933.512103 | 933.512103 | 0.327 | 1896 | 1944 | 0.989 |
| 24 | planted | hybrid | 462.987125 | 462.987125 | 0.659 | 1896 | 1944 | 0.979 |
| 24 | cross_planted | search | 40.357729 | 40.357729 | 2.782 | 10089 | 0 | 0.000 |
| 24 | cross_planted | flat | 112.278563 | 112.278563 | 1.000 | 746 | 764 | 0.970 |
| 24 | cross_planted | bucket | 336.001875 | 336.001875 | 0.334 | 746 | 764 | 0.989 |
| 24 | cross_planted | hybrid | 182.109875 | 182.109875 | 0.617 | 746 | 764 | 0.980 |
| 24 | unplanted | search | 72.694104 | 72.694104 | 4.683 | 17972 | 0 | 0.000 |
| 24 | unplanted | flat | 340.439563 | 340.439563 | 1.000 | 2139 | 2206 | 0.971 |
| 24 | unplanted | bucket | 1070.028333 | 1070.028333 | 0.318 | 2139 | 2206 | 0.990 |
| 24 | unplanted | hybrid | 516.598584 | 516.598584 | 0.659 | 2139 | 2206 | 0.981 |

Kernel-backed solvers have equal outcomes/models, logical counters and trace digests on every matched cell. The search-only control may explore a different tree. SAT models are checked against original equations; UNSAT verification requires a completed independent search reference. UNKNOWN is censored, not UNSAT.

This driver solves bounded generated Boolean systems. It does not benchmark the repository inherited-F4 implementation, accept curve targets, recover scalars or establish index-calculus performance. Production and cryptanalytic costs stay null.
