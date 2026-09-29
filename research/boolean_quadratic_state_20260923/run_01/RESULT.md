# Complete generic Boolean solve comparison

Evidence integrity: **PASS** across 48 cells and 3072 full-solve observations.

All results completed and independently verified: **True**. Full-solve gates against flat: **{'search': 'REJECTED', 'bucket': 'REJECTED', 'hybrid': 'REJECTED', 'small_flat': 'REJECTED', 'word_tail': 'REJECTED', 'merge_search': 'REJECTED', 'quadratic_state': 'PASS'}**. Hybrid versus fastest control: **REJECTED**.

The table reports cold solve plus result-validation milliseconds over two holdout seeds and the declared rotated repetitions. A completion cost is null if any sample in the group is censored or lacks the required verification. Observed capped work is retained separately. Ratios use pooled medians; acceptance uses paired intervals.

| Variables | Family | Arm | Complete cost (ms) | Observed cost (ms) | Flat / arm | Median nodes | Median kernel calls | Median kernel fraction |
|---:|---|---|---:|---:|---:|---:|---:|---:|
| 12 | planted | search | 0.166646 | 0.166646 | 2.235 | 94 | 0 | 0.000 |
| 12 | planted | flat | 0.372479 | 0.372479 | 1.000 | 8 | 8 | 0.956 |
| 12 | planted | bucket | 3.266062 | 3.266062 | 0.114 | 8 | 8 | 0.994 |
| 12 | planted | hybrid | 0.376875 | 0.376875 | 0.988 | 8 | 8 | 0.952 |
| 12 | planted | small_flat | 0.281500 | 0.281500 | 1.323 | 48 | 42 | 0.574 |
| 12 | planted | word_tail | 0.224605 | 0.224605 | 1.658 | 48 | 42 | 0.579 |
| 12 | planted | merge_search | 0.189979 | 0.189979 | 1.961 | 94 | 0 | 0.000 |
| 12 | planted | quadratic_state | 0.068458 | 0.068458 | 5.441 | 94 | 0 | 0.000 |
| 12 | cross_planted | search | 0.138416 | 0.138416 | 1.806 | 70 | 0 | 0.000 |
| 12 | cross_planted | flat | 0.250000 | 0.250000 | 1.000 | 4 | 4 | 0.952 |
| 12 | cross_planted | bucket | 2.185625 | 2.185625 | 0.114 | 4 | 4 | 0.993 |
| 12 | cross_planted | hybrid | 0.265000 | 0.265000 | 0.943 | 4 | 4 | 0.951 |
| 12 | cross_planted | small_flat | 0.163166 | 0.163166 | 1.532 | 28 | 22 | 0.583 |
| 12 | cross_planted | word_tail | 0.120105 | 0.120105 | 2.082 | 28 | 22 | 0.554 |
| 12 | cross_planted | merge_search | 0.126875 | 0.126875 | 1.970 | 70 | 0 | 0.000 |
| 12 | cross_planted | quadratic_state | 0.050750 | 0.050750 | 4.926 | 70 | 0 | 0.000 |
| 12 | unplanted | search | 0.286459 | 0.286459 | 1.764 | 153 | 0 | 0.000 |
| 12 | unplanted | flat | 0.505271 | 0.505271 | 1.000 | 9 | 9 | 0.951 |
| 12 | unplanted | bucket | 4.558624 | 4.558624 | 0.111 | 9 | 9 | 0.994 |
| 12 | unplanted | hybrid | 0.465542 | 0.465542 | 1.085 | 9 | 9 | 0.944 |
| 12 | unplanted | small_flat | 0.420188 | 0.420188 | 1.202 | 77 | 66 | 0.544 |
| 12 | unplanted | word_tail | 0.356083 | 0.356083 | 1.419 | 77 | 66 | 0.549 |
| 12 | unplanted | merge_search | 0.273083 | 0.273083 | 1.850 | 153 | 0 | 0.000 |
| 12 | unplanted | quadratic_state | 0.111146 | 0.111146 | 4.546 | 153 | 0 | 0.000 |
| 16 | planted | search | 2.262605 | 2.262605 | 1.838 | 812 | 0 | 0.000 |
| 16 | planted | flat | 4.158834 | 4.158834 | 1.000 | 68 | 68 | 0.943 |
| 16 | planted | bucket | 37.694250 | 37.694250 | 0.110 | 68 | 68 | 0.993 |
| 16 | planted | hybrid | 11.494771 | 11.494771 | 0.362 | 68 | 68 | 0.979 |
| 16 | planted | small_flat | 3.290688 | 3.290688 | 1.264 | 444 | 348 | 0.555 |
| 16 | planted | word_tail | 2.861854 | 2.861854 | 1.453 | 444 | 348 | 0.515 |
| 16 | planted | merge_search | 2.016291 | 2.016291 | 2.063 | 812 | 0 | 0.000 |
| 16 | planted | quadratic_state | 0.790624 | 0.790624 | 5.260 | 812 | 0 | 0.000 |
| 16 | cross_planted | search | 1.671792 | 1.671792 | 1.309 | 618 | 0 | 0.000 |
| 16 | cross_planted | flat | 2.188791 | 2.188791 | 1.000 | 32 | 32 | 0.946 |
| 16 | cross_planted | bucket | 21.602646 | 21.602646 | 0.101 | 32 | 32 | 0.994 |
| 16 | cross_planted | hybrid | 6.945146 | 6.945146 | 0.315 | 32 | 32 | 0.983 |
| 16 | cross_planted | small_flat | 1.908458 | 1.908458 | 1.147 | 270 | 225 | 0.537 |
| 16 | cross_planted | word_tail | 1.684625 | 1.684625 | 1.299 | 270 | 225 | 0.504 |
| 16 | cross_planted | merge_search | 1.513354 | 1.513354 | 1.446 | 618 | 0 | 0.000 |
| 16 | cross_planted | quadratic_state | 0.617687 | 0.617687 | 3.544 | 618 | 0 | 0.000 |
| 16 | unplanted | search | 4.117438 | 4.117438 | 1.766 | 1493 | 0 | 0.000 |
| 16 | unplanted | flat | 7.272667 | 7.272667 | 1.000 | 123 | 123 | 0.940 |
| 16 | unplanted | bucket | 67.710875 | 67.710875 | 0.107 | 123 | 123 | 0.993 |
| 16 | unplanted | hybrid | 16.115105 | 16.115105 | 0.451 | 123 | 123 | 0.972 |
| 16 | unplanted | small_flat | 6.120166 | 6.120166 | 1.188 | 805 | 622 | 0.568 |
| 16 | unplanted | word_tail | 5.205063 | 5.205063 | 1.397 | 805 | 622 | 0.509 |
| 16 | unplanted | merge_search | 3.685500 | 3.685500 | 1.973 | 1493 | 0 | 0.000 |
| 16 | unplanted | quadratic_state | 1.486438 | 1.486438 | 4.893 | 1493 | 0 | 0.000 |
| 20 | planted | search | 7.372333 | 7.372333 | 3.268 | 2108 | 0 | 0.000 |
| 20 | planted | flat | 24.094666 | 24.094666 | 1.000 | 208 | 222 | 0.961 |
| 20 | planted | bucket | 110.668521 | 110.668521 | 0.218 | 208 | 222 | 0.991 |
| 20 | planted | hybrid | 53.651854 | 53.651854 | 0.449 | 208 | 222 | 0.982 |
| 20 | planted | small_flat | 9.918209 | 9.918209 | 2.429 | 1136 | 579 | 0.527 |
| 20 | planted | word_tail | 8.263292 | 8.263292 | 2.916 | 1136 | 579 | 0.441 |
| 20 | planted | merge_search | 6.447105 | 6.447105 | 3.737 | 2108 | 0 | 0.000 |
| 20 | planted | quadratic_state | 2.696479 | 2.696479 | 8.936 | 2108 | 0 | 0.000 |
| 20 | cross_planted | search | 7.970708 | 7.970708 | 2.300 | 2286 | 0 | 0.000 |
| 20 | cross_planted | flat | 18.331542 | 18.331542 | 1.000 | 205 | 212 | 0.957 |
| 20 | cross_planted | bucket | 81.897688 | 81.897688 | 0.224 | 205 | 212 | 0.989 |
| 20 | cross_planted | hybrid | 37.303667 | 37.303667 | 0.491 | 205 | 212 | 0.979 |
| 20 | cross_planted | small_flat | 8.747417 | 8.747417 | 2.096 | 910 | 539 | 0.528 |
| 20 | cross_planted | word_tail | 7.245312 | 7.245312 | 2.530 | 910 | 539 | 0.442 |
| 20 | cross_planted | merge_search | 6.946209 | 6.946209 | 2.639 | 2286 | 0 | 0.000 |
| 20 | cross_planted | quadratic_state | 2.779396 | 2.779396 | 6.596 | 2286 | 0 | 0.000 |
| 20 | unplanted | search | 15.314312 | 15.314312 | 2.928 | 4224 | 0 | 0.000 |
| 20 | unplanted | flat | 44.843103 | 44.843103 | 1.000 | 340 | 342 | 0.967 |
| 20 | unplanted | bucket | 214.800146 | 214.800146 | 0.209 | 340 | 342 | 0.992 |
| 20 | unplanted | hybrid | 78.920438 | 78.920438 | 0.568 | 340 | 342 | 0.981 |
| 20 | unplanted | small_flat | 17.821771 | 17.821771 | 2.516 | 1980 | 1026 | 0.526 |
| 20 | unplanted | word_tail | 14.991937 | 14.991937 | 2.991 | 1980 | 1026 | 0.448 |
| 20 | unplanted | merge_search | 13.200396 | 13.200396 | 3.397 | 4224 | 0 | 0.000 |
| 20 | unplanted | quadratic_state | 5.384208 | 5.384208 | 8.329 | 4224 | 0 | 0.000 |
| 24 | planted | search | 61.967084 | 61.967084 | 4.315 | 14754 | 0 | 0.000 |
| 24 | planted | flat | 267.374687 | 267.374687 | 1.000 | 1572 | 1608 | 0.971 |
| 24 | planted | bucket | 792.604438 | 792.604438 | 0.337 | 1572 | 1608 | 0.990 |
| 24 | planted | hybrid | 345.599458 | 345.599458 | 0.774 | 1572 | 1608 | 0.978 |
| 24 | planted | small_flat | 72.297813 | 72.297813 | 3.698 | 11326 | 1590 | 0.304 |
| 24 | planted | word_tail | 64.188333 | 64.188333 | 4.165 | 11326 | 1590 | 0.215 |
| 24 | planted | merge_search | 54.717126 | 54.717126 | 4.886 | 14754 | 0 | 0.000 |
| 24 | planted | quadratic_state | 22.657812 | 22.657812 | 11.801 | 14754 | 0 | 0.000 |
| 24 | cross_planted | search | 61.672167 | 61.672167 | 2.488 | 14806 | 0 | 0.000 |
| 24 | cross_planted | flat | 153.435749 | 153.435749 | 1.000 | 888 | 905 | 0.973 |
| 24 | cross_planted | bucket | 503.017229 | 503.017229 | 0.305 | 888 | 905 | 0.991 |
| 24 | cross_planted | hybrid | 198.694854 | 198.694854 | 0.772 | 888 | 905 | 0.979 |
| 24 | cross_planted | small_flat | 72.154250 | 72.154250 | 2.126 | 11140 | 1699 | 0.330 |
| 24 | cross_planted | word_tail | 63.750271 | 63.750271 | 2.407 | 11140 | 1699 | 0.233 |
| 24 | cross_planted | merge_search | 54.606147 | 54.606147 | 2.810 | 14806 | 0 | 0.000 |
| 24 | cross_planted | quadratic_state | 22.433333 | 22.433333 | 6.840 | 14806 | 0 | 0.000 |
| 24 | unplanted | search | 73.110479 | 73.110479 | 4.550 | 17286 | 0 | 0.000 |
| 24 | unplanted | flat | 332.675771 | 332.675771 | 1.000 | 1776 | 1794 | 0.974 |
| 24 | unplanted | bucket | 1075.077958 | 1075.077958 | 0.309 | 1776 | 1794 | 0.991 |
| 24 | unplanted | hybrid | 503.827354 | 503.827354 | 0.660 | 1776 | 1794 | 0.983 |
| 24 | unplanted | small_flat | 83.772417 | 83.772417 | 3.971 | 13998 | 1544 | 0.201 |
| 24 | unplanted | word_tail | 77.538979 | 77.538979 | 4.290 | 13998 | 1544 | 0.145 |
| 24 | unplanted | merge_search | 64.138729 | 64.138729 | 5.187 | 17286 | 0 | 0.000 |
| 24 | unplanted | quadratic_state | 25.485771 | 25.485771 | 13.053 | 17286 | 0 | 0.000 |

Kernel-backed solvers have equal outcomes/models, logical counters and trace digests within each matched policy group. The search-only control may explore a different tree. SAT models are checked against original equations; UNSAT verification requires a completed independent search reference. UNKNOWN is censored, not UNSAT.

This driver solves bounded generated Boolean systems. It does not benchmark the repository inherited-F4 implementation, accept curve targets, recover scalars or establish index-calculus performance. Production and cryptanalytic costs stay null.

Selective one-word degree-2 full-solve gate: **REJECTED**. The small_flat and word_tail methods share a policy and must match exact traces/counters. They are compared separately with the always-degree-3 methods and search-only.

Ordered-specialization complete-solve gate: **REJECTED**. The search and merge_search arms have identical outcomes/models, logical work and trace digests. Only representation work in specialization changes. The gate compares merge_search with the fastest of all six prior full-solve methods.

Fixed-quadratic complete-solve gate: **PASS**. Compilation is charged to each cold solve. Search, merge_search and quadratic_state must match outcomes/models, logical counters and trace digests. The reference is the pointwise fastest of all seven prior methods.

| Variables | Family | Fastest prior / quadratic | 95% paired interval | Gate |
|---:|---|---:|---|---|
| 16 | planted | 2.458 | 2.428–2.519 | PASS |
| 16 | cross_planted | 2.457 | 2.368–2.539 | PASS |
| 16 | unplanted | 2.507 | 2.472–2.540 | PASS |
| 20 | planted | 2.423 | 2.369–2.459 | PASS |
| 20 | cross_planted | 2.472 | 2.426–2.501 | PASS |
| 20 | unplanted | 2.458 | 2.437–2.479 | PASS |
| 24 | planted | 2.412 | 2.405–2.431 | PASS |
| 24 | cross_planted | 2.422 | 2.411–2.447 | PASS |
| 24 | unplanted | 2.401 | 2.367–2.556 | PASS |

These paired intervals describe repeated timings on two fixed holdout fixtures per size/family, not a confidence interval over a population of Boolean systems. Fresh-worker RSS includes all arms and is not candidate-specific memory.
