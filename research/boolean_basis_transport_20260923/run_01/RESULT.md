# Complete generic Boolean solve comparison

Evidence integrity: **PASS** across 96 cells and 16224 full-solve observations.

All results completed and independently verified: **True**. Full-solve gates against flat: **{'search': 'REJECTED', 'bucket': 'REJECTED', 'hybrid': 'REJECTED', 'small_flat': 'REJECTED', 'word_tail': 'REJECTED', 'merge_search': 'REJECTED', 'quadratic_state': 'PASS', 'packed_state': 'PASS', 'basis_wide': 'REJECTED', 'tail_wide': 'PASS'}**. Hybrid versus fastest control: **REJECTED**.

The table reports cold solve plus result-validation milliseconds over two holdout seeds and the declared rotated repetitions. A completion cost is null if any sample in the group is censored or lacks the required verification. Observed capped work is retained separately. Ratios use pooled medians; acceptance uses paired intervals.

| Variables | Family | Arm | Complete cost (ms) | Observed cost (ms) | Flat / arm | Median nodes | Median kernel calls | Median kernel fraction |
|---:|---|---|---:|---:|---:|---:|---:|---:|
| 12 | planted | search | 0.302833 | 0.302833 | 1.653 | 136 | 0 | 0.000 |
| 12 | planted | flat | 0.500667 | 0.500667 | 1.000 | 8 | 8 | 0.948 |
| 12 | planted | bucket | 4.400313 | 4.400313 | 0.114 | 8 | 8 | 0.994 |
| 12 | planted | hybrid | 0.453813 | 0.453813 | 1.103 | 8 | 8 | 0.946 |
| 12 | planted | small_flat | 0.453209 | 0.453209 | 1.105 | 77 | 64 | 0.552 |
| 12 | planted | word_tail | 0.373229 | 0.373229 | 1.341 | 77 | 64 | 0.548 |
| 12 | planted | merge_search | 0.258104 | 0.258104 | 1.940 | 136 | 0 | 0.000 |
| 12 | planted | quadratic_state | 0.111042 | 0.111042 | 4.509 | 136 | 0 | 0.000 |
| 12 | planted | packed_state | 0.069854 | 0.069854 | 7.167 | 136 | 0 | 0.000 |
| 12 | planted | basis_list | 0.746209 | 0.746209 | 0.671 | 158 | 228 | 0.653 |
| 12 | planted | basis_wide | 0.136542 | 0.136542 | 3.667 | 158 | 228 | 0.518 |
| 12 | planted | tail_list | 0.494667 | 0.494667 | 1.012 | 104 | 145 | 0.579 |
| 12 | planted | tail_wide | 0.113604 | 0.113604 | 4.407 | 104 | 145 | 0.385 |
| 12 | cross_planted | search | 0.223271 | 0.223271 | 1.467 | 118 | 0 | 0.000 |
| 12 | cross_planted | flat | 0.327521 | 0.327521 | 1.000 | 5 | 6 | 0.945 |
| 12 | cross_planted | bucket | 2.776041 | 2.776041 | 0.118 | 5 | 6 | 0.993 |
| 12 | cross_planted | hybrid | 0.298667 | 0.298667 | 1.097 | 5 | 6 | 0.945 |
| 12 | cross_planted | small_flat | 0.210520 | 0.210520 | 1.556 | 38 | 28 | 0.539 |
| 12 | cross_planted | word_tail | 0.171230 | 0.171230 | 1.913 | 38 | 28 | 0.507 |
| 12 | cross_planted | merge_search | 0.209354 | 0.209354 | 1.564 | 118 | 0 | 0.000 |
| 12 | cross_planted | quadratic_state | 0.082542 | 0.082542 | 3.968 | 118 | 0 | 0.000 |
| 12 | cross_planted | packed_state | 0.050875 | 0.050875 | 6.438 | 118 | 0 | 0.000 |
| 12 | cross_planted | basis_list | 0.366812 | 0.366812 | 0.893 | 74 | 109 | 0.646 |
| 12 | cross_planted | basis_wide | 0.067209 | 0.067209 | 4.873 | 74 | 109 | 0.502 |
| 12 | cross_planted | tail_list | 0.277791 | 0.277791 | 1.179 | 64 | 91 | 0.592 |
| 12 | cross_planted | tail_wide | 0.070042 | 0.070042 | 4.676 | 64 | 91 | 0.382 |
| 12 | unplanted | search | 0.531313 | 0.531313 | 1.194 | 242 | 0 | 0.000 |
| 12 | unplanted | flat | 0.634166 | 0.634166 | 1.000 | 14 | 14 | 0.940 |
| 12 | unplanted | bucket | 6.454063 | 6.454063 | 0.098 | 14 | 14 | 0.993 |
| 12 | unplanted | hybrid | 0.607458 | 0.607458 | 1.044 | 14 | 14 | 0.936 |
| 12 | unplanted | small_flat | 0.735250 | 0.735250 | 0.863 | 142 | 111 | 0.522 |
| 12 | unplanted | word_tail | 0.681250 | 0.681250 | 0.931 | 142 | 111 | 0.504 |
| 12 | unplanted | merge_search | 0.465667 | 0.465667 | 1.362 | 242 | 0 | 0.000 |
| 12 | unplanted | quadratic_state | 0.208000 | 0.208000 | 3.049 | 242 | 0 | 0.000 |
| 12 | unplanted | packed_state | 0.132292 | 0.132292 | 4.794 | 242 | 0 | 0.000 |
| 12 | unplanted | basis_list | 0.885021 | 0.885021 | 0.717 | 180 | 259 | 0.660 |
| 12 | unplanted | basis_wide | 0.175813 | 0.175813 | 3.607 | 180 | 259 | 0.528 |
| 12 | unplanted | tail_list | 0.875645 | 0.875645 | 0.724 | 174 | 240 | 0.592 |
| 12 | unplanted | tail_wide | 0.220833 | 0.220833 | 2.872 | 174 | 240 | 0.397 |
| 16 | planted | search | 2.834250 | 2.834250 | 1.944 | 964 | 0 | 0.000 |
| 16 | planted | flat | 5.508938 | 5.508938 | 1.000 | 88 | 89 | 0.942 |
| 16 | planted | bucket | 53.075333 | 53.075333 | 0.104 | 88 | 89 | 0.993 |
| 16 | planted | hybrid | 12.765583 | 12.765583 | 0.432 | 88 | 89 | 0.975 |
| 16 | planted | small_flat | 3.926271 | 3.926271 | 1.403 | 464 | 374 | 0.574 |
| 16 | planted | word_tail | 3.222771 | 3.222771 | 1.709 | 464 | 374 | 0.503 |
| 16 | planted | merge_search | 2.455312 | 2.455312 | 2.244 | 964 | 0 | 0.000 |
| 16 | planted | quadratic_state | 0.995271 | 0.995271 | 5.535 | 964 | 0 | 0.000 |
| 16 | planted | packed_state | 0.611521 | 0.611521 | 9.009 | 964 | 0 | 0.000 |
| 16 | planted | basis_list | 4.813645 | 4.813645 | 1.144 | 764 | 1063 | 0.679 |
| 16 | planted | basis_wide | 0.929895 | 0.929895 | 5.924 | 764 | 1063 | 0.531 |
| 16 | planted | tail_list | 4.878188 | 4.878188 | 1.129 | 712 | 1008 | 0.592 |
| 16 | planted | tail_wide | 1.139063 | 1.139063 | 4.836 | 712 | 1008 | 0.430 |
| 16 | cross_planted | search | 2.494791 | 2.494791 | 1.452 | 894 | 0 | 0.000 |
| 16 | cross_planted | flat | 3.622250 | 3.622250 | 1.000 | 60 | 60 | 0.940 |
| 16 | cross_planted | bucket | 33.183062 | 33.183062 | 0.109 | 60 | 60 | 0.993 |
| 16 | cross_planted | hybrid | 8.450021 | 8.450021 | 0.429 | 60 | 60 | 0.974 |
| 16 | cross_planted | small_flat | 2.246583 | 2.246583 | 1.612 | 262 | 240 | 0.563 |
| 16 | cross_planted | word_tail | 1.947896 | 1.947896 | 1.860 | 262 | 240 | 0.513 |
| 16 | cross_planted | merge_search | 2.116354 | 2.116354 | 1.712 | 894 | 0 | 0.000 |
| 16 | cross_planted | quadratic_state | 0.880709 | 0.880709 | 4.113 | 894 | 0 | 0.000 |
| 16 | cross_planted | packed_state | 0.536458 | 0.536458 | 6.752 | 894 | 0 | 0.000 |
| 16 | cross_planted | basis_list | 2.244772 | 2.244772 | 1.614 | 385 | 529 | 0.681 |
| 16 | cross_planted | basis_wide | 0.452584 | 0.452584 | 8.003 | 385 | 529 | 0.523 |
| 16 | cross_planted | tail_list | 3.261188 | 3.261188 | 1.111 | 528 | 752 | 0.591 |
| 16 | cross_planted | tail_wide | 0.802562 | 0.802562 | 4.513 | 528 | 752 | 0.427 |
| 16 | unplanted | search | 2.314062 | 2.314062 | 2.015 | 840 | 0 | 0.000 |
| 16 | unplanted | flat | 4.663021 | 4.663021 | 1.000 | 72 | 72 | 0.949 |
| 16 | unplanted | bucket | 44.527229 | 44.527229 | 0.105 | 72 | 72 | 0.993 |
| 16 | unplanted | hybrid | 10.918750 | 10.918750 | 0.427 | 72 | 72 | 0.979 |
| 16 | unplanted | small_flat | 3.303041 | 3.303041 | 1.412 | 442 | 350 | 0.577 |
| 16 | unplanted | word_tail | 2.977063 | 2.977063 | 1.566 | 442 | 350 | 0.514 |
| 16 | unplanted | merge_search | 2.012749 | 2.012749 | 2.317 | 840 | 0 | 0.000 |
| 16 | unplanted | quadratic_state | 0.820312 | 0.820312 | 5.684 | 840 | 0 | 0.000 |
| 16 | unplanted | packed_state | 0.503791 | 0.503791 | 9.256 | 840 | 0 | 0.000 |
| 16 | unplanted | basis_list | 8.852291 | 8.852291 | 0.527 | 1486 | 1996 | 0.693 |
| 16 | unplanted | basis_wide | 1.765333 | 1.765333 | 2.641 | 1486 | 1996 | 0.546 |
| 16 | unplanted | tail_list | 4.413625 | 4.413625 | 1.057 | 596 | 848 | 0.605 |
| 16 | unplanted | tail_wide | 0.962938 | 0.962938 | 4.842 | 596 | 848 | 0.434 |
| 20 | planted | search | 12.455854 | 12.455854 | 2.967 | 3408 | 0 | 0.000 |
| 20 | planted | flat | 36.951270 | 36.951270 | 1.000 | 414 | 429 | 0.953 |
| 20 | planted | bucket | 190.182229 | 190.182229 | 0.194 | 414 | 429 | 0.990 |
| 20 | planted | hybrid | 72.237646 | 72.237646 | 0.512 | 414 | 429 | 0.974 |
| 20 | planted | small_flat | 16.639563 | 16.639563 | 2.221 | 1790 | 1032 | 0.545 |
| 20 | planted | word_tail | 13.778792 | 13.778792 | 2.682 | 1790 | 1032 | 0.465 |
| 20 | planted | merge_search | 10.828729 | 10.828729 | 3.412 | 3408 | 0 | 0.000 |
| 20 | planted | quadratic_state | 4.402479 | 4.402479 | 8.393 | 3408 | 0 | 0.000 |
| 20 | planted | packed_state | 2.721291 | 2.721291 | 13.579 | 3408 | 0 | 0.000 |
| 20 | planted | basis_list | 7.600500 | 7.600500 | 4.862 | 956 | 1346 | 0.708 |
| 20 | planted | basis_wide | 1.603438 | 1.603438 | 23.045 | 956 | 1346 | 0.566 |
| 20 | planted | tail_list | 20.022313 | 20.022313 | 1.846 | 2681 | 4014 | 0.599 |
| 20 | planted | tail_wide | 6.063000 | 6.063000 | 6.095 | 2681 | 4014 | 0.465 |
| 20 | cross_planted | search | 10.320729 | 10.320729 | 2.161 | 3002 | 0 | 0.000 |
| 20 | cross_planted | flat | 22.304333 | 22.304333 | 1.000 | 272 | 276 | 0.951 |
| 20 | cross_planted | bucket | 106.741104 | 106.741104 | 0.209 | 272 | 276 | 0.989 |
| 20 | cross_planted | hybrid | 38.931562 | 38.931562 | 0.573 | 272 | 276 | 0.972 |
| 20 | cross_planted | small_flat | 12.119041 | 12.119041 | 1.840 | 1454 | 796 | 0.530 |
| 20 | cross_planted | word_tail | 10.560146 | 10.560146 | 2.112 | 1454 | 796 | 0.452 |
| 20 | cross_planted | merge_search | 8.689187 | 8.689187 | 2.567 | 3002 | 0 | 0.000 |
| 20 | cross_planted | quadratic_state | 3.708000 | 3.708000 | 6.015 | 3002 | 0 | 0.000 |
| 20 | cross_planted | packed_state | 2.209583 | 2.209583 | 10.094 | 3002 | 0 | 0.000 |
| 20 | cross_planted | basis_list | 5.400749 | 5.400749 | 4.130 | 712 | 1054 | 0.694 |
| 20 | cross_planted | basis_wide | 1.179688 | 1.179688 | 18.907 | 712 | 1054 | 0.553 |
| 20 | cross_planted | tail_list | 13.078687 | 13.078687 | 1.705 | 1862 | 2780 | 0.595 |
| 20 | cross_planted | tail_wide | 3.897771 | 3.897771 | 5.722 | 1862 | 2780 | 0.465 |
| 20 | unplanted | search | 23.312521 | 23.312521 | 2.576 | 6656 | 0 | 0.000 |
| 20 | unplanted | flat | 60.049479 | 60.049479 | 1.000 | 695 | 709 | 0.953 |
| 20 | unplanted | bucket | 292.625563 | 292.625563 | 0.205 | 695 | 709 | 0.989 |
| 20 | unplanted | hybrid | 102.147958 | 102.147958 | 0.588 | 695 | 709 | 0.972 |
| 20 | unplanted | small_flat | 31.926167 | 31.926167 | 1.881 | 3575 | 2046 | 0.546 |
| 20 | unplanted | word_tail | 27.365500 | 27.365500 | 2.194 | 3575 | 2046 | 0.476 |
| 20 | unplanted | merge_search | 20.179105 | 20.179105 | 2.976 | 6656 | 0 | 0.000 |
| 20 | unplanted | quadratic_state | 8.213458 | 8.213458 | 7.311 | 6656 | 0 | 0.000 |
| 20 | unplanted | packed_state | 4.996729 | 4.996729 | 12.018 | 6656 | 0 | 0.000 |
| 20 | unplanted | basis_list | 60.206729 | 60.206729 | 0.997 | 7919 | 11370 | 0.716 |
| 20 | unplanted | basis_wide | 13.144146 | 13.144146 | 4.569 | 7919 | 11370 | 0.576 |
| 20 | unplanted | tail_list | 36.029479 | 36.029479 | 1.667 | 5395 | 7922 | 0.608 |
| 20 | unplanted | tail_wide | 11.064771 | 11.064771 | 5.427 | 5395 | 7922 | 0.467 |
| 24 | planted | search | 66.217584 | 66.217584 | 4.624 | 15729 | 0 | 0.000 |
| 24 | planted | flat | 306.182646 | 306.182646 | 1.000 | 1630 | 1654 | 0.973 |
| 24 | planted | bucket | 1209.731875 | 1209.731875 | 0.253 | 1630 | 1654 | 0.992 |
| 24 | planted | hybrid | 564.794188 | 564.794188 | 0.542 | 1630 | 1654 | 0.984 |
| 24 | planted | small_flat | 77.067354 | 77.067354 | 3.973 | 11972 | 1768 | 0.293 |
| 24 | planted | word_tail | 68.757750 | 68.757750 | 4.453 | 11972 | 1768 | 0.216 |
| 24 | planted | merge_search | 58.147124 | 58.147124 | 5.266 | 15729 | 0 | 0.000 |
| 24 | planted | quadratic_state | 24.245251 | 24.245251 | 12.629 | 15729 | 0 | 0.000 |
| 24 | planted | packed_state | 14.673687 | 14.673687 | 20.866 | 15729 | 0 | 0.000 |
| 24 | planted | basis_list | 125.438604 | 125.438604 | 2.441 | 10814 | 16178 | 0.728 |
| 24 | planted | basis_wide | 32.977396 | 32.977396 | 9.285 | 10814 | 16178 | 0.591 |
| 24 | planted | tail_list | 115.883104 | 115.883104 | 2.642 | 13052 | 19389 | 0.623 |
| 24 | planted | tail_wide | 36.544917 | 36.544917 | 8.378 | 13052 | 19389 | 0.520 |
| 24 | cross_planted | search | 40.276479 | 40.276479 | 4.496 | 10170 | 0 | 0.000 |
| 24 | cross_planted | flat | 181.077229 | 181.077229 | 1.000 | 1252 | 1270 | 0.969 |
| 24 | cross_planted | bucket | 577.337396 | 577.337396 | 0.314 | 1252 | 1270 | 0.989 |
| 24 | cross_planted | hybrid | 258.394334 | 258.394334 | 0.701 | 1252 | 1270 | 0.977 |
| 24 | cross_planted | small_flat | 48.026376 | 48.026376 | 3.770 | 7693 | 1191 | 0.297 |
| 24 | cross_planted | word_tail | 42.481395 | 42.481395 | 4.263 | 7693 | 1191 | 0.220 |
| 24 | cross_planted | merge_search | 35.623875 | 35.623875 | 5.083 | 10170 | 0 | 0.000 |
| 24 | cross_planted | quadratic_state | 14.992396 | 14.992396 | 12.078 | 10170 | 0 | 0.000 |
| 24 | cross_planted | packed_state | 9.026291 | 9.026291 | 20.061 | 10170 | 0 | 0.000 |
| 24 | cross_planted | basis_list | 51.116333 | 51.116333 | 3.542 | 5640 | 8430 | 0.724 |
| 24 | cross_planted | basis_wide | 12.047854 | 12.047854 | 15.030 | 5640 | 8430 | 0.586 |
| 24 | cross_planted | tail_list | 74.345251 | 74.345251 | 2.436 | 8654 | 12665 | 0.621 |
| 24 | cross_planted | tail_wide | 23.088125 | 23.088125 | 7.843 | 8654 | 12665 | 0.513 |
| 24 | unplanted | search | 77.027376 | 77.027376 | 4.102 | 18809 | 0 | 0.000 |
| 24 | unplanted | flat | 315.948980 | 315.948980 | 1.000 | 1951 | 1964 | 0.971 |
| 24 | unplanted | bucket | 1029.723521 | 1029.723521 | 0.307 | 1951 | 1964 | 0.990 |
| 24 | unplanted | hybrid | 402.785188 | 402.785188 | 0.784 | 1951 | 1964 | 0.978 |
| 24 | unplanted | small_flat | 89.620125 | 89.620125 | 3.525 | 12978 | 2619 | 0.337 |
| 24 | unplanted | word_tail | 79.346271 | 79.346271 | 3.982 | 12978 | 2619 | 0.258 |
| 24 | unplanted | merge_search | 67.757313 | 67.757313 | 4.663 | 18809 | 0 | 0.000 |
| 24 | unplanted | quadratic_state | 27.964833 | 27.964833 | 11.298 | 18809 | 0 | 0.000 |
| 24 | unplanted | packed_state | 16.798105 | 16.798105 | 18.809 | 18809 | 0 | 0.000 |
| 24 | unplanted | basis_list | 268.847458 | 268.847458 | 1.175 | 27930 | 42857 | 0.726 |
| 24 | unplanted | basis_wide | 64.268354 | 64.268354 | 4.916 | 27930 | 42857 | 0.585 |
| 24 | unplanted | tail_list | 135.799709 | 135.799709 | 2.327 | 15913 | 24270 | 0.613 |
| 24 | unplanted | tail_wide | 44.602750 | 44.602750 | 7.084 | 15913 | 24270 | 0.505 |

Kernel-backed solvers have equal outcomes/models, logical counters and trace digests within each matched policy group. The search-only control may explore a different tree. SAT models are checked against original equations; UNSAT verification requires a completed independent search reference. UNKNOWN is censored, not UNSAT.

This driver solves bounded generated Boolean systems. It does not benchmark the repository inherited-F4 implementation, accept curve targets, recover scalars or establish index-calculus performance. Production and cryptanalytic costs stay null.

Selective one-word degree-2 full-solve gate: **REJECTED**. The small_flat and word_tail methods share a policy and must match exact traces/counters. They are compared separately with the always-degree-3 methods and search-only.

Ordered-specialization complete-solve gate: **REJECTED**. The search and merge_search arms have identical outcomes/models, logical work and trace digests. Only representation work in specialization changes. The gate compares merge_search with the fastest of all six prior full-solve methods.

Fixed-quadratic complete-solve gate: **REJECTED**. Compilation is charged to each cold solve. Search, merge_search and quadratic_state must match outcomes/models, logical counters and trace digests. The reference is the pointwise fastest of all seven prior methods.

| Variables | Family | Fastest prior / quadratic | 95% paired interval | Gate |
|---:|---|---:|---|---|
| 16 | planted | 2.455 | 2.420–2.521 | PASS |
| 16 | cross_planted | 2.117 | 1.945–2.304 | REJECTED |
| 16 | unplanted | 2.411 | 2.315–2.475 | PASS |
| 20 | planted | 2.390 | 2.375–2.440 | PASS |
| 20 | cross_planted | 2.388 | 2.353–2.417 | PASS |
| 20 | unplanted | 2.449 | 2.438–2.469 | PASS |
| 24 | planted | 2.397 | 2.382–2.410 | PASS |
| 24 | cross_planted | 2.386 | 2.375–2.417 | PASS |
| 24 | unplanted | 2.427 | 2.398–2.441 | PASS |

These paired intervals describe repeated timings on two fixed holdout fixtures per size/family, not a confidence interval over a population of Boolean systems. Fresh-worker RSS includes all arms and is not candidate-specific memory.

Packed-state current-frontier 2x gate: **REJECTED**. Incremental improvement gate (>1x lower bound): **PASS**. Cumulative historical-seven 2x gate: **PASS**.

Both fresh holdouts and the entire prior confirmation regression grid must pass. The current frontier includes the already faster quadratic_state method; the historical-seven comparison cannot substitute for it.

| Split | Variables | Family | Reference | Paired ratio | 95% interval | >2x gate |
|---|---:|---|---|---:|---|---|
| regression | 16 | planted | current_frontier | 1.648 | 1.629–1.670 | REJECTED |
| regression | 16 | planted | historical_seven | 3.744 | 3.634–3.823 | PASS |
| regression | 16 | cross_planted | current_frontier | 1.660 | 1.645–1.680 | REJECTED |
| regression | 16 | cross_planted | historical_seven | 3.564 | 3.443–3.786 | PASS |
| regression | 16 | unplanted | current_frontier | 1.623 | 1.596–1.647 | REJECTED |
| regression | 16 | unplanted | historical_seven | 4.058 | 3.998–4.111 | PASS |
| regression | 20 | planted | current_frontier | 1.656 | 1.647–1.671 | REJECTED |
| regression | 20 | planted | historical_seven | 4.020 | 3.993–4.059 | PASS |
| regression | 20 | cross_planted | current_frontier | 1.666 | 1.646–1.686 | REJECTED |
| regression | 20 | cross_planted | historical_seven | 3.954 | 3.840–4.028 | PASS |
| regression | 20 | unplanted | current_frontier | 1.621 | 1.607–1.641 | REJECTED |
| regression | 20 | unplanted | historical_seven | 3.980 | 3.910–4.048 | PASS |
| regression | 24 | planted | current_frontier | 1.655 | 1.646–1.676 | REJECTED |
| regression | 24 | planted | historical_seven | 3.965 | 3.947–3.991 | PASS |
| regression | 24 | cross_planted | current_frontier | 1.655 | 1.636–1.695 | REJECTED |
| regression | 24 | cross_planted | historical_seven | 3.959 | 3.932–3.980 | PASS |
| regression | 24 | unplanted | current_frontier | 1.653 | 1.648–1.659 | REJECTED |
| regression | 24 | unplanted | historical_seven | 3.989 | 3.972–4.000 | PASS |
| holdout | 16 | planted | current_frontier | 1.651 | 1.605–1.678 | REJECTED |
| holdout | 16 | planted | historical_seven | 4.081 | 4.015–4.123 | PASS |
| holdout | 16 | cross_planted | current_frontier | 1.600 | 1.570–1.657 | REJECTED |
| holdout | 16 | cross_planted | historical_seven | 3.567 | 3.006–3.822 | PASS |
| holdout | 16 | unplanted | current_frontier | 1.650 | 1.635–1.685 | REJECTED |
| holdout | 16 | unplanted | historical_seven | 3.944 | 3.877–4.077 | PASS |
| holdout | 20 | planted | current_frontier | 1.660 | 1.601–1.687 | REJECTED |
| holdout | 20 | planted | historical_seven | 3.958 | 3.906–4.016 | PASS |
| holdout | 20 | cross_planted | current_frontier | 1.664 | 1.632–1.700 | REJECTED |
| holdout | 20 | cross_planted | historical_seven | 3.984 | 3.908–4.025 | PASS |
| holdout | 20 | unplanted | current_frontier | 1.606 | 1.582–1.623 | REJECTED |
| holdout | 20 | unplanted | historical_seven | 3.956 | 3.911–4.000 | PASS |
| holdout | 24 | planted | current_frontier | 1.648 | 1.642–1.679 | REJECTED |
| holdout | 24 | planted | historical_seven | 3.978 | 3.940–4.018 | PASS |
| holdout | 24 | cross_planted | current_frontier | 1.670 | 1.654–1.682 | REJECTED |
| holdout | 24 | cross_planted | historical_seven | 3.981 | 3.936–4.009 | PASS |
| holdout | 24 | unplanted | current_frontier | 1.666 | 1.649–1.678 | REJECTED |
| holdout | 24 | unplanted | historical_seven | 4.016 | 3.977–4.060 | PASS |

## Transported coefficient-row-space comparison

basis_list/basis_wide share a full-RREF policy and branch on that basis. tail_list/tail_wide keep quadratic rows in echelon form, reduce only the affine tail, and branch using carried original-equation residuals. The backends in each pair must match every logical counter, model and trace. New inference policies need not match the older search tree.

Source-row and distinct-column counts price each reduction input, including zero rows. Specialization counts cover every carried representation. Basis histograms use nominal unassigned-variable width. All costs, including the extra original-equation state in the tail policy, are inside cold solve timing.

| Candidate | Dramatic >2x strongest-reference gate | Any improvement strongest-reference gate | Matched backend >1.05x gate |
|---|---|---|---|
| basis_wide | REJECTED | REJECTED | PASS |
| tail_wide | REJECTED | REJECTED | PASS |

| Candidate | Split | Variables | Family | Reference | Paired ratio | 95% interval | Gate |
|---|---|---:|---|---|---:|---|---|
| basis_wide | regression | 16 | planted | strongest | 0.278 | 0.063–0.535 | REJECTED |
| basis_wide | regression | 16 | planted | same_policy | 4.987 | 4.914–5.071 | PASS |
| basis_wide | regression | 16 | cross_planted | strongest | 0.678 | 0.641–0.694 | REJECTED |
| basis_wide | regression | 16 | cross_planted | same_policy | 4.687 | 4.580–4.821 | PASS |
| basis_wide | regression | 16 | unplanted | strongest | 0.534 | 0.498–0.562 | REJECTED |
| basis_wide | regression | 16 | unplanted | same_policy | 5.054 | 4.979–5.144 | PASS |
| basis_wide | regression | 20 | planted | strongest | 0.773 | 0.355–1.181 | REJECTED |
| basis_wide | regression | 20 | planted | same_policy | 4.693 | 4.618–4.761 | PASS |
| basis_wide | regression | 20 | cross_planted | strongest | 1.246 | 1.113–1.436 | REJECTED |
| basis_wide | regression | 20 | cross_planted | same_policy | 4.685 | 4.601–4.712 | PASS |
| basis_wide | regression | 20 | unplanted | strongest | 0.430 | 0.392–0.473 | REJECTED |
| basis_wide | regression | 20 | unplanted | same_policy | 4.808 | 4.776–4.886 | PASS |
| basis_wide | regression | 24 | planted | strongest | 0.264 | 0.205–0.327 | REJECTED |
| basis_wide | regression | 24 | planted | same_policy | 4.332 | 4.202–4.471 | PASS |
| basis_wide | regression | 24 | cross_planted | strongest | 0.233 | 0.170–0.302 | REJECTED |
| basis_wide | regression | 24 | cross_planted | same_policy | 4.430 | 4.328–4.518 | PASS |
| basis_wide | regression | 24 | unplanted | strongest | 0.267 | 0.254–0.285 | REJECTED |
| basis_wide | regression | 24 | unplanted | same_policy | 4.431 | 4.375–4.478 | PASS |
| basis_wide | holdout | 16 | planted | strongest | 0.628 | 0.494–0.734 | REJECTED |
| basis_wide | holdout | 16 | planted | same_policy | 5.202 | 5.094–5.300 | PASS |
| basis_wide | holdout | 16 | cross_planted | strongest | 1.883 | 0.502–3.397 | REJECTED |
| basis_wide | holdout | 16 | cross_planted | same_policy | 4.915 | 4.779–4.966 | PASS |
| basis_wide | holdout | 16 | unplanted | strongest | 0.274 | 0.069–0.504 | REJECTED |
| basis_wide | holdout | 16 | unplanted | same_policy | 5.072 | 4.910–5.165 | PASS |
| basis_wide | holdout | 20 | planted | strongest | 2.667 | 1.177–4.551 | REJECTED |
| basis_wide | holdout | 20 | planted | same_policy | 4.612 | 4.541–4.729 | PASS |
| basis_wide | holdout | 20 | cross_planted | strongest | 2.654 | 1.355–4.572 | REJECTED |
| basis_wide | holdout | 20 | cross_planted | same_policy | 4.616 | 4.546–4.650 | PASS |
| basis_wide | holdout | 20 | unplanted | strongest | 0.391 | 0.305–0.505 | REJECTED |
| basis_wide | holdout | 20 | unplanted | same_policy | 4.563 | 4.534–4.589 | PASS |
| basis_wide | holdout | 24 | planted | strongest | 0.463 | 0.404–1.047 | REJECTED |
| basis_wide | holdout | 24 | planted | same_policy | 4.358 | 4.335–4.615 | PASS |
| basis_wide | holdout | 24 | cross_planted | strongest | 1.572 | 0.331–2.936 | REJECTED |
| basis_wide | holdout | 24 | cross_planted | same_policy | 4.278 | 4.204–4.363 | PASS |
| basis_wide | holdout | 24 | unplanted | strongest | 0.199 | 0.070–0.332 | REJECTED |
| basis_wide | holdout | 24 | unplanted | same_policy | 4.268 | 4.217–4.405 | PASS |
| tail_wide | regression | 16 | planted | strongest | 0.516 | 0.504–0.524 | REJECTED |
| tail_wide | regression | 16 | planted | same_policy | 4.143 | 4.063–4.275 | PASS |
| tail_wide | regression | 16 | cross_planted | strongest | 0.952 | 0.834–1.103 | REJECTED |
| tail_wide | regression | 16 | cross_planted | same_policy | 4.214 | 4.155–4.262 | PASS |
| tail_wide | regression | 16 | unplanted | strongest | 0.536 | 0.522–0.557 | REJECTED |
| tail_wide | regression | 16 | unplanted | same_policy | 3.982 | 3.903–4.082 | PASS |
| tail_wide | regression | 20 | planted | strongest | 0.364 | 0.343–0.408 | REJECTED |
| tail_wide | regression | 20 | planted | same_policy | 3.550 | 3.499–3.596 | PASS |
| tail_wide | regression | 20 | cross_planted | strongest | 0.630 | 0.552–0.701 | REJECTED |
| tail_wide | regression | 20 | cross_planted | same_policy | 3.390 | 3.305–3.431 | PASS |
| tail_wide | regression | 20 | unplanted | strongest | 0.467 | 0.445–0.473 | REJECTED |
| tail_wide | regression | 20 | unplanted | same_policy | 3.475 | 3.351–3.552 | PASS |
| tail_wide | regression | 24 | planted | strongest | 0.391 | 0.383–0.403 | REJECTED |
| tail_wide | regression | 24 | planted | same_policy | 3.181 | 3.085–3.218 | PASS |
| tail_wide | regression | 24 | cross_planted | strongest | 0.534 | 0.410–0.634 | REJECTED |
| tail_wide | regression | 24 | cross_planted | same_policy | 3.066 | 3.022–3.191 | PASS |
| tail_wide | regression | 24 | unplanted | strongest | 0.394 | 0.393–0.397 | REJECTED |
| tail_wide | regression | 24 | unplanted | same_policy | 3.250 | 3.208–3.283 | PASS |
| tail_wide | holdout | 16 | planted | strongest | 0.516 | 0.505–0.534 | REJECTED |
| tail_wide | holdout | 16 | planted | same_policy | 4.150 | 4.083–4.284 | PASS |
| tail_wide | holdout | 16 | cross_planted | strongest | 0.365 | 0.208–0.571 | REJECTED |
| tail_wide | holdout | 16 | cross_planted | same_policy | 4.027 | 3.861–4.126 | PASS |
| tail_wide | holdout | 16 | unplanted | strongest | 0.530 | 0.519–0.538 | REJECTED |
| tail_wide | holdout | 16 | unplanted | same_policy | 4.591 | 4.478–4.689 | PASS |
| tail_wide | holdout | 20 | planted | strongest | 0.235 | 0.081–0.389 | REJECTED |
| tail_wide | holdout | 20 | planted | same_policy | 3.333 | 3.220–3.386 | PASS |
| tail_wide | holdout | 20 | cross_planted | strongest | 0.255 | 0.064–0.454 | REJECTED |
| tail_wide | holdout | 20 | cross_planted | same_policy | 3.339 | 3.291–3.363 | PASS |
| tail_wide | holdout | 20 | unplanted | strongest | 0.465 | 0.441–0.482 | REJECTED |
| tail_wide | holdout | 20 | unplanted | same_policy | 3.295 | 3.138–3.346 | PASS |
| tail_wide | holdout | 24 | planted | strongest | 0.384 | 0.379–0.388 | REJECTED |
| tail_wide | holdout | 24 | planted | same_policy | 3.221 | 3.114–3.244 | PASS |
| tail_wide | holdout | 24 | cross_planted | strongest | 0.238 | 0.144–0.337 | REJECTED |
| tail_wide | holdout | 24 | cross_planted | same_policy | 3.233 | 3.206–3.256 | PASS |
| tail_wide | holdout | 24 | unplanted | strongest | 0.412 | 0.374–0.472 | REJECTED |
| tail_wide | holdout | 24 | unplanted | same_policy | 3.032 | 2.959–3.061 | PASS |
