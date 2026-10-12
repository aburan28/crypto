# Complete generic Boolean solve comparison

Evidence integrity: **PASS** across 48 cells and 1728 full-solve observations.

All results completed and independently verified: **True**. Full-solve gates against flat: **{'search': 'REJECTED', 'bucket': 'REJECTED', 'hybrid': 'REJECTED', 'small_flat': 'REJECTED', 'word_tail': 'REJECTED'}**. Hybrid versus fastest control: **REJECTED**.

The table reports cold solve plus result-validation milliseconds over two holdout seeds and the declared rotated repetitions. A completion cost is null if any sample in the group is censored or lacks the required verification. Observed capped work is retained separately. Ratios use pooled medians; acceptance uses paired intervals.

| Variables | Family | Arm | Complete cost (ms) | Observed cost (ms) | Flat / arm | Median nodes | Median kernel calls | Median kernel fraction |
|---:|---|---|---:|---:|---:|---:|---:|---:|
| 12 | planted | search | 0.224104 | 0.224104 | 2.059 | 111 | 0 | 0.000 |
| 12 | planted | flat | 0.461500 | 0.461500 | 1.000 | 10 | 10 | 0.950 |
| 12 | planted | bucket | 4.878188 | 4.878188 | 0.095 | 10 | 10 | 0.994 |
| 12 | planted | hybrid | 0.427312 | 0.427312 | 1.080 | 10 | 10 | 0.944 |
| 12 | planted | small_flat | 0.386000 | 0.386000 | 1.196 | 70 | 53 | 0.557 |
| 12 | planted | word_tail | 0.335896 | 0.335896 | 1.374 | 70 | 53 | 0.514 |
| 12 | cross_planted | search | 0.230292 | 0.230292 | 1.199 | 121 | 0 | 0.000 |
| 12 | cross_planted | flat | 0.276125 | 0.276125 | 1.000 | 4 | 5 | 0.957 |
| 12 | cross_planted | bucket | 2.487563 | 2.487563 | 0.111 | 4 | 5 | 0.995 |
| 12 | cross_planted | hybrid | 0.270500 | 0.270500 | 1.021 | 4 | 5 | 0.957 |
| 12 | cross_planted | small_flat | 0.135833 | 0.135833 | 2.033 | 24 | 18 | 0.603 |
| 12 | cross_planted | word_tail | 0.095354 | 0.095354 | 2.896 | 24 | 18 | 0.591 |
| 12 | unplanted | search | 0.636271 | 0.636271 | 1.067 | 322 | 0 | 0.000 |
| 12 | unplanted | flat | 0.678771 | 0.678771 | 1.000 | 15 | 15 | 0.944 |
| 12 | unplanted | bucket | 7.158291 | 7.158291 | 0.095 | 15 | 15 | 0.994 |
| 12 | unplanted | hybrid | 0.669542 | 0.669542 | 1.014 | 15 | 15 | 0.938 |
| 12 | unplanted | small_flat | 0.899021 | 0.899021 | 0.755 | 168 | 139 | 0.515 |
| 12 | unplanted | word_tail | 0.860250 | 0.860250 | 0.789 | 168 | 139 | 0.507 |
| 16 | planted | search | 0.682625 | 0.682625 | 2.741 | 267 | 0 | 0.000 |
| 16 | planted | flat | 1.871145 | 1.871145 | 1.000 | 26 | 26 | 0.959 |
| 16 | planted | bucket | 15.113583 | 15.113583 | 0.124 | 26 | 26 | 0.994 |
| 16 | planted | hybrid | 6.041605 | 6.041605 | 0.310 | 26 | 26 | 0.986 |
| 16 | planted | small_flat | 1.057688 | 1.057688 | 1.769 | 136 | 104 | 0.557 |
| 16 | planted | word_tail | 0.826896 | 0.826896 | 2.263 | 136 | 104 | 0.495 |
| 16 | cross_planted | search | 0.660021 | 0.660021 | 2.026 | 262 | 0 | 0.000 |
| 16 | cross_planted | flat | 1.337063 | 1.337063 | 1.000 | 18 | 20 | 0.949 |
| 16 | cross_planted | bucket | 9.940646 | 9.940646 | 0.135 | 18 | 20 | 0.993 |
| 16 | cross_planted | hybrid | 4.477105 | 4.477105 | 0.299 | 18 | 20 | 0.985 |
| 16 | cross_planted | small_flat | 0.656687 | 0.656687 | 2.036 | 96 | 67 | 0.514 |
| 16 | cross_planted | word_tail | 0.574041 | 0.574041 | 2.329 | 96 | 67 | 0.495 |
| 16 | unplanted | search | 3.465876 | 3.465876 | 2.013 | 1293 | 0 | 0.000 |
| 16 | unplanted | flat | 6.975875 | 6.975875 | 1.000 | 116 | 116 | 0.942 |
| 16 | unplanted | bucket | 57.818729 | 57.818729 | 0.121 | 116 | 116 | 0.992 |
| 16 | unplanted | hybrid | 13.886333 | 13.886333 | 0.502 | 116 | 116 | 0.969 |
| 16 | unplanted | small_flat | 5.112876 | 5.112876 | 1.364 | 659 | 512 | 0.554 |
| 16 | unplanted | word_tail | 4.455646 | 4.455646 | 1.566 | 659 | 512 | 0.499 |
| 20 | planted | search | 11.692730 | 11.692730 | 2.429 | 3550 | 0 | 0.000 |
| 20 | planted | flat | 28.400562 | 28.400562 | 1.000 | 296 | 324 | 0.953 |
| 20 | planted | bucket | 152.070771 | 152.070771 | 0.187 | 296 | 324 | 0.990 |
| 20 | planted | hybrid | 72.616917 | 72.616917 | 0.391 | 296 | 324 | 0.981 |
| 20 | planted | small_flat | 16.390625 | 16.390625 | 1.733 | 1824 | 1118 | 0.562 |
| 20 | planted | word_tail | 13.559875 | 13.559875 | 2.094 | 1824 | 1118 | 0.485 |
| 20 | cross_planted | search | 7.440271 | 7.440271 | 1.455 | 2430 | 0 | 0.000 |
| 20 | cross_planted | flat | 10.824687 | 10.824687 | 1.000 | 138 | 140 | 0.955 |
| 20 | cross_planted | bucket | 52.515667 | 52.515667 | 0.206 | 138 | 140 | 0.990 |
| 20 | cross_planted | hybrid | 26.586105 | 26.586105 | 0.407 | 138 | 140 | 0.980 |
| 20 | cross_planted | small_flat | 8.561084 | 8.561084 | 1.264 | 1048 | 623 | 0.543 |
| 20 | cross_planted | word_tail | 7.159020 | 7.159020 | 1.512 | 1048 | 623 | 0.468 |
| 20 | unplanted | search | 12.441167 | 12.441167 | 2.690 | 3745 | 0 | 0.000 |
| 20 | unplanted | flat | 33.472771 | 33.472771 | 1.000 | 336 | 359 | 0.955 |
| 20 | unplanted | bucket | 213.884938 | 213.884938 | 0.156 | 336 | 359 | 0.992 |
| 20 | unplanted | hybrid | 90.568000 | 90.568000 | 0.370 | 336 | 359 | 0.982 |
| 20 | unplanted | small_flat | 17.698771 | 17.698771 | 1.891 | 2071 | 1127 | 0.548 |
| 20 | unplanted | word_tail | 14.705333 | 14.705333 | 2.276 | 2071 | 1127 | 0.469 |
| 24 | planted | search | 45.432188 | 45.432188 | 3.324 | 11376 | 0 | 0.000 |
| 24 | planted | flat | 151.023375 | 151.023375 | 1.000 | 1258 | 1276 | 0.968 |
| 24 | planted | bucket | 507.155666 | 507.155666 | 0.298 | 1258 | 1276 | 0.990 |
| 24 | planted | hybrid | 239.517979 | 239.517979 | 0.631 | 1258 | 1276 | 0.979 |
| 24 | planted | small_flat | 53.877459 | 53.877459 | 2.803 | 8420 | 1423 | 0.354 |
| 24 | planted | word_tail | 48.051041 | 48.051041 | 3.143 | 8420 | 1423 | 0.274 |
| 24 | cross_planted | search | 50.946604 | 50.946604 | 1.775 | 13044 | 0 | 0.000 |
| 24 | cross_planted | flat | 90.425938 | 90.425938 | 1.000 | 798 | 806 | 0.968 |
| 24 | cross_planted | bucket | 292.500084 | 292.500084 | 0.309 | 798 | 806 | 0.990 |
| 24 | cross_planted | hybrid | 141.813625 | 141.813625 | 0.638 | 798 | 806 | 0.978 |
| 24 | cross_planted | small_flat | 60.047146 | 60.047146 | 1.506 | 9626 | 1606 | 0.337 |
| 24 | cross_planted | word_tail | 53.416021 | 53.416021 | 1.693 | 9626 | 1606 | 0.258 |
| 24 | unplanted | search | 153.465000 | 153.465000 | 3.586 | 38351 | 0 | 0.000 |
| 24 | unplanted | flat | 550.373458 | 550.373458 | 1.000 | 3938 | 3982 | 0.968 |
| 24 | unplanted | bucket | 1912.102729 | 1912.102729 | 0.288 | 3938 | 3982 | 0.990 |
| 24 | unplanted | hybrid | 773.957688 | 773.957688 | 0.711 | 3938 | 3982 | 0.976 |
| 24 | unplanted | small_flat | 182.218105 | 182.218105 | 3.020 | 26706 | 5394 | 0.363 |
| 24 | unplanted | word_tail | 160.860834 | 160.860834 | 3.421 | 26706 | 5394 | 0.284 |

Kernel-backed solvers have equal outcomes/models, logical counters and trace digests within each matched policy group. The search-only control may explore a different tree. SAT models are checked against original equations; UNSAT verification requires a completed independent search reference. UNKNOWN is censored, not UNSAT.

This driver solves bounded generated Boolean systems. It does not benchmark the repository inherited-F4 implementation, accept curve targets, recover scalars or establish index-calculus performance. Production and cryptanalytic costs stay null.

Selective one-word degree-2 full-solve gate: **REJECTED**. The small_flat and word_tail methods share a policy and must match exact traces/counters. They are compared separately with the always-degree-3 methods and search-only.
