# Complete generic Boolean solve comparison

Evidence integrity: **PASS** across 120 cells and 58080 full-solve observations.

All results completed and independently verified: **True**. Full-solve gates against flat: **{'search': 'REJECTED', 'bucket': 'REJECTED', 'hybrid': 'REJECTED', 'small_flat': 'REJECTED', 'word_tail': 'REJECTED', 'merge_search': 'REJECTED', 'quadratic_state': 'REJECTED', 'packed_state': 'PASS', 'basis_wide': 'REJECTED', 'tail_wide': 'PASS', 'affine_sl_basis_fast': 'PASS', 'gray_simd': 'PASS', 'packed_gray12_simd': 'PASS', 'packed_gray16_simd': 'PASS'}**. Hybrid versus fastest control: **REJECTED**.

The table reports cold solve plus result-validation milliseconds over two holdout seeds and the declared rotated repetitions. A completion cost is null if any sample in the group is censored or lacks the required verification. Observed capped work is retained separately. Ratios use pooled medians; acceptance uses paired intervals.

| Variables | Family | Arm | Complete cost (ms) | Observed cost (ms) | Flat / arm | Median nodes | Median kernel calls | Median kernel fraction |
|---:|---|---|---:|---:|---:|---:|---:|---:|
| 12 | planted | search | 0.321250 | 0.321250 | 1.671 | 154 | 0 | 0.000 |
| 12 | planted | flat | 0.536937 | 0.536937 | 1.000 | 11 | 11 | 0.945 |
| 12 | planted | bucket | 4.983437 | 4.983437 | 0.108 | 11 | 11 | 0.994 |
| 12 | planted | hybrid | 0.510645 | 0.510645 | 1.051 | 11 | 11 | 0.942 |
| 12 | planted | small_flat | 0.494250 | 0.494250 | 1.086 | 84 | 70 | 0.564 |
| 12 | planted | word_tail | 0.416021 | 0.416021 | 1.291 | 84 | 70 | 0.533 |
| 12 | planted | merge_search | 0.292750 | 0.292750 | 1.834 | 154 | 0 | 0.000 |
| 12 | planted | quadratic_state | 0.116834 | 0.116834 | 4.596 | 154 | 0 | 0.000 |
| 12 | planted | packed_state | 0.075979 | 0.075979 | 7.067 | 154 | 0 | 0.000 |
| 12 | planted | basis_list | 0.400708 | 0.400708 | 1.340 | 90 | 121 | 0.649 |
| 12 | planted | basis_wide | 0.069917 | 0.069917 | 7.680 | 90 | 121 | 0.482 |
| 12 | planted | tail_list | 0.583479 | 0.583479 | 0.920 | 110 | 148 | 0.590 |
| 12 | planted | tail_wide | 0.118813 | 0.118813 | 4.519 | 110 | 148 | 0.392 |
| 12 | planted | packed_untraced | 0.051312 | 0.051312 | 10.464 | 154 | 0 | 0.000 |
| 12 | planted | affine_sl_basis_list | 0.673604 | 0.673604 | 0.797 | 27 | 69 | 0.677 |
| 12 | planted | affine_sl_basis_fast | 0.077458 | 0.077458 | 6.932 | 27 | 69 | 0.223 |
| 12 | planted | gray_scalar | 0.003187 | 0.003187 | 168.451 | 0 | 0 | 0.000 |
| 12 | planted | gray_simd | 0.001375 | 0.001375 | 390.500 | 0 | 0 | 0.000 |
| 12 | planted | packed_gray12_scalar | 0.005188 | 0.005188 | 103.506 | 1 | 0 | 0.000 |
| 12 | planted | packed_gray12_simd | 0.003438 | 0.003438 | 156.200 | 1 | 0 | 0.000 |
| 12 | planted | packed_gray16_scalar | 0.004834 | 0.004834 | 111.087 | 1 | 0 | 0.000 |
| 12 | planted | packed_gray16_simd | 0.003125 | 0.003125 | 171.820 | 1 | 0 | 0.000 |
| 12 | cross_planted | search | 0.191291 | 0.191291 | 1.598 | 106 | 0 | 0.000 |
| 12 | cross_planted | flat | 0.305688 | 0.305688 | 1.000 | 5 | 6 | 0.953 |
| 12 | cross_planted | bucket | 2.384313 | 2.384313 | 0.128 | 5 | 6 | 0.994 |
| 12 | cross_planted | hybrid | 0.279375 | 0.279375 | 1.094 | 5 | 6 | 0.953 |
| 12 | cross_planted | small_flat | 0.210729 | 0.210729 | 1.451 | 38 | 32 | 0.619 |
| 12 | cross_planted | word_tail | 0.160604 | 0.160604 | 1.903 | 38 | 32 | 0.590 |
| 12 | cross_planted | merge_search | 0.178896 | 0.178896 | 1.709 | 106 | 0 | 0.000 |
| 12 | cross_planted | quadratic_state | 0.069771 | 0.069771 | 4.381 | 106 | 0 | 0.000 |
| 12 | cross_planted | packed_state | 0.046583 | 0.046583 | 6.562 | 106 | 0 | 0.000 |
| 12 | cross_planted | basis_list | 0.120188 | 0.120188 | 2.543 | 22 | 31 | 0.661 |
| 12 | cross_planted | basis_wide | 0.022687 | 0.022687 | 13.474 | 22 | 31 | 0.485 |
| 12 | cross_planted | tail_list | 0.213667 | 0.213667 | 1.431 | 44 | 60 | 0.610 |
| 12 | cross_planted | tail_wide | 0.045562 | 0.045562 | 6.709 | 44 | 60 | 0.393 |
| 12 | cross_planted | packed_untraced | 0.029792 | 0.029792 | 10.261 | 106 | 0 | 0.000 |
| 12 | cross_planted | affine_sl_basis_list | 0.434541 | 0.434541 | 0.703 | 16 | 42 | 0.673 |
| 12 | cross_planted | affine_sl_basis_fast | 0.057313 | 0.057313 | 5.334 | 16 | 42 | 0.196 |
| 12 | cross_planted | gray_scalar | 0.003125 | 0.003125 | 97.820 | 0 | 0 | 0.000 |
| 12 | cross_planted | gray_simd | 0.001291 | 0.001291 | 236.692 | 0 | 0 | 0.000 |
| 12 | cross_planted | packed_gray12_scalar | 0.004729 | 0.004729 | 64.641 | 1 | 0 | 0.000 |
| 12 | cross_planted | packed_gray12_simd | 0.003396 | 0.003396 | 90.014 | 1 | 0 | 0.000 |
| 12 | cross_planted | packed_gray16_scalar | 0.004604 | 0.004604 | 66.396 | 1 | 0 | 0.000 |
| 12 | cross_planted | packed_gray16_simd | 0.003167 | 0.003167 | 96.538 | 1 | 0 | 0.000 |
| 12 | unplanted | search | 0.522645 | 0.522645 | 1.205 | 240 | 0 | 0.000 |
| 12 | unplanted | flat | 0.629812 | 0.629812 | 1.000 | 14 | 14 | 0.937 |
| 12 | unplanted | bucket | 6.265854 | 6.265854 | 0.101 | 14 | 14 | 0.993 |
| 12 | unplanted | hybrid | 0.602395 | 0.602395 | 1.046 | 14 | 14 | 0.935 |
| 12 | unplanted | small_flat | 0.705937 | 0.705937 | 0.892 | 125 | 106 | 0.551 |
| 12 | unplanted | word_tail | 0.653645 | 0.653645 | 0.964 | 125 | 106 | 0.531 |
| 12 | unplanted | merge_search | 0.457625 | 0.457625 | 1.376 | 240 | 0 | 0.000 |
| 12 | unplanted | quadratic_state | 0.187688 | 0.187688 | 3.356 | 240 | 0 | 0.000 |
| 12 | unplanted | packed_state | 0.120625 | 0.120625 | 5.221 | 240 | 0 | 0.000 |
| 12 | unplanted | basis_list | 0.764229 | 0.764229 | 0.824 | 164 | 225 | 0.659 |
| 12 | unplanted | basis_wide | 0.133500 | 0.133500 | 4.718 | 164 | 225 | 0.516 |
| 12 | unplanted | tail_list | 0.968458 | 0.968458 | 0.650 | 186 | 248 | 0.584 |
| 12 | unplanted | tail_wide | 0.198187 | 0.198187 | 3.178 | 186 | 248 | 0.401 |
| 12 | unplanted | packed_untraced | 0.078021 | 0.078021 | 8.072 | 240 | 0 | 0.000 |
| 12 | unplanted | affine_sl_basis_list | 1.800146 | 1.800146 | 0.350 | 62 | 161 | 0.703 |
| 12 | unplanted | affine_sl_basis_fast | 0.156812 | 0.156812 | 4.016 | 62 | 161 | 0.262 |
| 12 | unplanted | gray_scalar | 0.006167 | 0.006167 | 102.135 | 0 | 0 | 0.000 |
| 12 | unplanted | gray_simd | 0.001792 | 0.001792 | 351.458 | 0 | 0 | 0.000 |
| 12 | unplanted | packed_gray12_scalar | 0.007938 | 0.007938 | 79.341 | 1 | 0 | 0.000 |
| 12 | unplanted | packed_gray12_simd | 0.003917 | 0.003917 | 160.790 | 1 | 0 | 0.000 |
| 12 | unplanted | packed_gray16_scalar | 0.007521 | 0.007521 | 83.741 | 1 | 0 | 0.000 |
| 12 | unplanted | packed_gray16_simd | 0.003667 | 0.003667 | 171.751 | 1 | 0 | 0.000 |
| 16 | planted | search | 2.098583 | 2.098583 | 2.040 | 758 | 0 | 0.000 |
| 16 | planted | flat | 4.280188 | 4.280188 | 1.000 | 70 | 70 | 0.944 |
| 16 | planted | bucket | 38.209916 | 38.209916 | 0.112 | 70 | 70 | 0.993 |
| 16 | planted | hybrid | 10.091292 | 10.091292 | 0.424 | 70 | 70 | 0.977 |
| 16 | planted | small_flat | 3.098562 | 3.098562 | 1.381 | 416 | 316 | 0.563 |
| 16 | planted | word_tail | 2.608958 | 2.608958 | 1.641 | 416 | 316 | 0.506 |
| 16 | planted | merge_search | 1.865458 | 1.865458 | 2.294 | 758 | 0 | 0.000 |
| 16 | planted | quadratic_state | 0.738958 | 0.738958 | 5.792 | 758 | 0 | 0.000 |
| 16 | planted | packed_state | 0.465438 | 0.465438 | 9.196 | 758 | 0 | 0.000 |
| 16 | planted | basis_list | 2.527479 | 2.527479 | 1.693 | 391 | 546 | 0.684 |
| 16 | planted | basis_wide | 0.480770 | 0.480770 | 8.903 | 391 | 546 | 0.523 |
| 16 | planted | tail_list | 3.693271 | 3.693271 | 1.159 | 577 | 838 | 0.588 |
| 16 | planted | tail_wide | 0.898479 | 0.898479 | 4.764 | 577 | 838 | 0.417 |
| 16 | planted | packed_untraced | 0.317812 | 0.317812 | 13.468 | 758 | 0 | 0.000 |
| 16 | planted | affine_sl_basis_list | 7.202021 | 7.202021 | 0.594 | 194 | 536 | 0.653 |
| 16 | planted | affine_sl_basis_fast | 0.551562 | 0.551562 | 7.760 | 194 | 536 | 0.244 |
| 16 | planted | gray_scalar | 0.018021 | 0.018021 | 237.511 | 0 | 0 | 0.000 |
| 16 | planted | gray_simd | 0.003708 | 0.003708 | 1154.156 | 0 | 0 | 0.000 |
| 16 | planted | packed_gray12_scalar | 0.067917 | 0.067917 | 63.021 | 18 | 0 | 0.000 |
| 16 | planted | packed_gray12_simd | 0.031313 | 0.031313 | 136.693 | 18 | 0 | 0.000 |
| 16 | planted | packed_gray16_scalar | 0.018333 | 0.018333 | 233.463 | 1 | 0 | 0.000 |
| 16 | planted | packed_gray16_simd | 0.006104 | 0.006104 | 701.210 | 1 | 0 | 0.000 |
| 16 | cross_planted | search | 2.751959 | 2.751959 | 0.676 | 1036 | 0 | 0.000 |
| 16 | cross_planted | flat | 1.861479 | 1.861479 | 1.000 | 30 | 31 | 0.947 |
| 16 | cross_planted | bucket | 15.556938 | 15.556938 | 0.120 | 30 | 31 | 0.993 |
| 16 | cross_planted | hybrid | 5.703271 | 5.703271 | 0.326 | 30 | 31 | 0.982 |
| 16 | cross_planted | small_flat | 2.120001 | 2.120001 | 0.878 | 308 | 225 | 0.528 |
| 16 | cross_planted | word_tail | 1.838208 | 1.838208 | 1.013 | 308 | 225 | 0.483 |
| 16 | cross_planted | merge_search | 2.450458 | 2.450458 | 0.760 | 1036 | 0 | 0.000 |
| 16 | cross_planted | quadratic_state | 0.976812 | 0.976812 | 1.906 | 1036 | 0 | 0.000 |
| 16 | cross_planted | packed_state | 0.617563 | 0.617563 | 3.014 | 1036 | 0 | 0.000 |
| 16 | cross_planted | basis_list | 1.059521 | 1.059521 | 1.757 | 179 | 235 | 0.663 |
| 16 | cross_planted | basis_wide | 0.205750 | 0.205750 | 9.047 | 179 | 235 | 0.517 |
| 16 | cross_planted | tail_list | 1.529479 | 1.529479 | 1.217 | 254 | 362 | 0.587 |
| 16 | cross_planted | tail_wide | 0.377937 | 0.377937 | 4.925 | 254 | 362 | 0.417 |
| 16 | cross_planted | packed_untraced | 0.417480 | 0.417480 | 4.459 | 1036 | 0 | 0.000 |
| 16 | cross_planted | affine_sl_basis_list | 2.314792 | 2.314792 | 0.804 | 54 | 164 | 0.614 |
| 16 | cross_planted | affine_sl_basis_fast | 0.212438 | 0.212438 | 8.762 | 54 | 164 | 0.218 |
| 16 | cross_planted | gray_scalar | 0.017250 | 0.017250 | 107.912 | 0 | 0 | 0.000 |
| 16 | cross_planted | gray_simd | 0.003749 | 0.003749 | 496.461 | 0 | 0 | 0.000 |
| 16 | cross_planted | packed_gray12_scalar | 0.089291 | 0.089291 | 20.847 | 24 | 0 | 0.000 |
| 16 | cross_planted | packed_gray12_simd | 0.039334 | 0.039334 | 47.326 | 24 | 0 | 0.000 |
| 16 | cross_planted | packed_gray16_scalar | 0.017687 | 0.017687 | 105.243 | 1 | 0 | 0.000 |
| 16 | cross_planted | packed_gray16_simd | 0.006083 | 0.006083 | 306.013 | 1 | 0 | 0.000 |
| 16 | unplanted | search | 3.250938 | 3.250938 | 1.989 | 1162 | 0 | 0.000 |
| 16 | unplanted | flat | 6.465541 | 6.465541 | 1.000 | 110 | 110 | 0.941 |
| 16 | unplanted | bucket | 55.929708 | 55.929708 | 0.116 | 110 | 110 | 0.993 |
| 16 | unplanted | hybrid | 12.510958 | 12.510958 | 0.517 | 110 | 110 | 0.968 |
| 16 | unplanted | small_flat | 4.663750 | 4.663750 | 1.386 | 647 | 486 | 0.551 |
| 16 | unplanted | word_tail | 4.103854 | 4.103854 | 1.575 | 647 | 486 | 0.502 |
| 16 | unplanted | merge_search | 2.864000 | 2.864000 | 2.258 | 1162 | 0 | 0.000 |
| 16 | unplanted | quadratic_state | 1.146604 | 1.146604 | 5.639 | 1162 | 0 | 0.000 |
| 16 | unplanted | packed_state | 0.725479 | 0.725479 | 8.912 | 1162 | 0 | 0.000 |
| 16 | unplanted | basis_list | 8.989062 | 8.989062 | 0.719 | 1447 | 2034 | 0.687 |
| 16 | unplanted | basis_wide | 1.850646 | 1.850646 | 3.494 | 1447 | 2034 | 0.526 |
| 16 | unplanted | tail_list | 5.830312 | 5.830312 | 1.109 | 921 | 1270 | 0.595 |
| 16 | unplanted | tail_wide | 1.412813 | 1.412813 | 4.576 | 921 | 1270 | 0.428 |
| 16 | unplanted | packed_untraced | 0.481438 | 0.481438 | 13.430 | 1162 | 0 | 0.000 |
| 16 | unplanted | affine_sl_basis_list | 13.915688 | 13.915688 | 0.465 | 338 | 991 | 0.641 |
| 16 | unplanted | affine_sl_basis_fast | 1.072437 | 1.072437 | 6.029 | 338 | 991 | 0.270 |
| 16 | unplanted | gray_scalar | 0.070562 | 0.070562 | 91.629 | 0 | 0 | 0.000 |
| 16 | unplanted | gray_simd | 0.012854 | 0.012854 | 502.998 | 0 | 0 | 0.000 |
| 16 | unplanted | packed_gray12_scalar | 0.105105 | 0.105105 | 61.515 | 28 | 0 | 0.000 |
| 16 | unplanted | packed_gray12_simd | 0.045271 | 0.045271 | 142.819 | 28 | 0 | 0.000 |
| 16 | unplanted | packed_gray16_scalar | 0.075541 | 0.075541 | 85.589 | 1 | 0 | 0.000 |
| 16 | unplanted | packed_gray16_simd | 0.014958 | 0.014958 | 432.232 | 1 | 0 | 0.000 |
| 20 | planted | search | 8.976646 | 8.976646 | 2.876 | 2662 | 0 | 0.000 |
| 20 | planted | flat | 25.820167 | 25.820167 | 1.000 | 263 | 272 | 0.957 |
| 20 | planted | bucket | 136.731917 | 136.731917 | 0.189 | 263 | 272 | 0.991 |
| 20 | planted | hybrid | 73.580625 | 73.580625 | 0.351 | 263 | 272 | 0.983 |
| 20 | planted | small_flat | 12.380979 | 12.380979 | 2.085 | 1400 | 794 | 0.550 |
| 20 | planted | word_tail | 10.084812 | 10.084812 | 2.560 | 1400 | 794 | 0.470 |
| 20 | planted | merge_search | 7.885834 | 7.885834 | 3.274 | 2662 | 0 | 0.000 |
| 20 | planted | quadratic_state | 3.236041 | 3.236041 | 7.979 | 2662 | 0 | 0.000 |
| 20 | planted | packed_state | 2.027750 | 2.027750 | 12.733 | 2662 | 0 | 0.000 |
| 20 | planted | basis_list | 32.994000 | 32.994000 | 0.783 | 4177 | 5958 | 0.715 |
| 20 | planted | basis_wide | 7.051355 | 7.051355 | 3.662 | 4177 | 5958 | 0.567 |
| 20 | planted | tail_list | 15.626292 | 15.626292 | 1.652 | 2090 | 3034 | 0.602 |
| 20 | planted | tail_wide | 4.387521 | 4.387521 | 5.885 | 2090 | 3034 | 0.476 |
| 20 | planted | packed_untraced | 1.414770 | 1.414770 | 18.250 | 2662 | 0 | 0.000 |
| 20 | planted | affine_sl_basis_list | 37.392105 | 37.392105 | 0.691 | 840 | 2302 | 0.580 |
| 20 | planted | affine_sl_basis_fast | 3.444875 | 3.444875 | 7.495 | 840 | 2302 | 0.229 |
| 20 | planted | gray_scalar | 0.565895 | 0.565895 | 45.627 | 0 | 0 | 0.000 |
| 20 | planted | gray_simd | 0.093458 | 0.093458 | 276.276 | 0 | 0 | 0.000 |
| 20 | planted | packed_gray12_scalar | 0.790791 | 0.790791 | 32.651 | 218 | 0 | 0.000 |
| 20 | planted | packed_gray12_simd | 0.332271 | 0.332271 | 77.708 | 218 | 0 | 0.000 |
| 20 | planted | packed_gray16_scalar | 0.552959 | 0.552959 | 46.695 | 15 | 0 | 0.000 |
| 20 | planted | packed_gray16_simd | 0.106897 | 0.106897 | 241.544 | 15 | 0 | 0.000 |
| 20 | cross_planted | search | 6.938834 | 6.938834 | 2.354 | 2110 | 0 | 0.000 |
| 20 | cross_planted | flat | 16.336709 | 16.336709 | 1.000 | 189 | 192 | 0.954 |
| 20 | cross_planted | bucket | 77.658292 | 77.658292 | 0.210 | 189 | 192 | 0.990 |
| 20 | cross_planted | hybrid | 38.875125 | 38.875125 | 0.420 | 189 | 192 | 0.981 |
| 20 | cross_planted | small_flat | 7.928500 | 7.928500 | 2.061 | 902 | 512 | 0.533 |
| 20 | cross_planted | word_tail | 6.546208 | 6.546208 | 2.496 | 902 | 512 | 0.457 |
| 20 | cross_planted | merge_search | 6.085146 | 6.085146 | 2.685 | 2110 | 0 | 0.000 |
| 20 | cross_planted | quadratic_state | 2.527896 | 2.527896 | 6.463 | 2110 | 0 | 0.000 |
| 20 | cross_planted | packed_state | 1.565792 | 1.565792 | 10.434 | 2110 | 0 | 0.000 |
| 20 | cross_planted | basis_list | 24.434354 | 24.434354 | 0.669 | 3119 | 4468 | 0.712 |
| 20 | cross_planted | basis_wide | 5.061396 | 5.061396 | 3.228 | 3119 | 4468 | 0.560 |
| 20 | cross_planted | tail_list | 9.693041 | 9.693041 | 1.685 | 1307 | 1924 | 0.607 |
| 20 | cross_planted | tail_wide | 2.640125 | 2.640125 | 6.188 | 1307 | 1924 | 0.468 |
| 20 | cross_planted | packed_untraced | 1.082125 | 1.082125 | 15.097 | 2110 | 0 | 0.000 |
| 20 | cross_planted | affine_sl_basis_list | 31.889396 | 31.889396 | 0.512 | 708 | 1980 | 0.571 |
| 20 | cross_planted | affine_sl_basis_fast | 2.970458 | 2.970458 | 5.500 | 708 | 1980 | 0.231 |
| 20 | cross_planted | gray_scalar | 0.574875 | 0.574875 | 28.418 | 0 | 0 | 0.000 |
| 20 | cross_planted | gray_simd | 0.096895 | 0.096895 | 168.601 | 0 | 0 | 0.000 |
| 20 | cross_planted | packed_gray12_scalar | 0.664959 | 0.664959 | 24.568 | 188 | 0 | 0.000 |
| 20 | cross_planted | packed_gray12_simd | 0.296854 | 0.296854 | 55.033 | 188 | 0 | 0.000 |
| 20 | cross_planted | packed_gray16_scalar | 0.477500 | 0.477500 | 34.213 | 12 | 0 | 0.000 |
| 20 | cross_planted | packed_gray16_simd | 0.092646 | 0.092646 | 176.336 | 12 | 0 | 0.000 |
| 20 | unplanted | search | 25.016479 | 25.016479 | 1.999 | 7152 | 0 | 0.000 |
| 20 | unplanted | flat | 50.010354 | 50.010354 | 1.000 | 521 | 528 | 0.954 |
| 20 | unplanted | bucket | 362.859146 | 362.859146 | 0.138 | 521 | 528 | 0.993 |
| 20 | unplanted | hybrid | 124.848333 | 124.848333 | 0.401 | 521 | 528 | 0.980 |
| 20 | unplanted | small_flat | 31.518292 | 31.518292 | 1.587 | 3619 | 2116 | 0.534 |
| 20 | unplanted | word_tail | 27.804562 | 27.804562 | 1.799 | 3619 | 2116 | 0.479 |
| 20 | unplanted | merge_search | 22.216355 | 22.216355 | 2.251 | 7152 | 0 | 0.000 |
| 20 | unplanted | quadratic_state | 8.579854 | 8.579854 | 5.829 | 7152 | 0 | 0.000 |
| 20 | unplanted | packed_state | 5.473645 | 5.473645 | 9.137 | 7152 | 0 | 0.000 |
| 20 | unplanted | basis_list | 75.101646 | 75.101646 | 0.666 | 9594 | 13587 | 0.710 |
| 20 | unplanted | basis_wide | 15.313792 | 15.313792 | 3.266 | 9594 | 13587 | 0.565 |
| 20 | unplanted | tail_list | 42.707167 | 42.707167 | 1.171 | 5766 | 8328 | 0.601 |
| 20 | unplanted | tail_wide | 11.922021 | 11.922021 | 4.195 | 5766 | 8328 | 0.476 |
| 20 | unplanted | packed_untraced | 3.773125 | 3.773125 | 13.254 | 7152 | 0 | 0.000 |
| 20 | unplanted | affine_sl_basis_list | 93.609125 | 93.609125 | 0.534 | 1979 | 5419 | 0.556 |
| 20 | unplanted | affine_sl_basis_fast | 8.176250 | 8.176250 | 6.117 | 1979 | 5419 | 0.230 |
| 20 | unplanted | gray_scalar | 1.330000 | 1.330000 | 37.602 | 0 | 0 | 0.000 |
| 20 | unplanted | gray_simd | 0.200583 | 0.200583 | 249.325 | 0 | 0 | 0.000 |
| 20 | unplanted | packed_gray12_scalar | 1.920187 | 1.920187 | 26.045 | 511 | 0 | 0.000 |
| 20 | unplanted | packed_gray12_simd | 0.793125 | 0.793125 | 63.055 | 511 | 0 | 0.000 |
| 20 | unplanted | packed_gray16_scalar | 1.379000 | 1.379000 | 36.266 | 31 | 0 | 0.000 |
| 20 | unplanted | packed_gray16_simd | 0.245063 | 0.245063 | 204.071 | 31 | 0 | 0.000 |
| 24 | planted | search | 79.291375 | 79.291375 | 3.295 | 19200 | 0 | 0.000 |
| 24 | planted | flat | 261.296500 | 261.296500 | 1.000 | 1784 | 1799 | 0.970 |
| 24 | planted | bucket | 932.561834 | 932.561834 | 0.280 | 1784 | 1799 | 0.990 |
| 24 | planted | hybrid | 323.089250 | 323.089250 | 0.809 | 1784 | 1799 | 0.977 |
| 24 | planted | small_flat | 90.098479 | 90.098479 | 2.900 | 12256 | 3042 | 0.344 |
| 24 | planted | word_tail | 78.273562 | 78.273562 | 3.338 | 12256 | 3042 | 0.265 |
| 24 | planted | merge_search | 69.558688 | 69.558688 | 3.756 | 19200 | 0 | 0.000 |
| 24 | planted | quadratic_state | 28.146959 | 28.146959 | 9.283 | 19200 | 0 | 0.000 |
| 24 | planted | packed_state | 17.568792 | 17.568792 | 14.873 | 19200 | 0 | 0.000 |
| 24 | planted | basis_list | 237.259479 | 237.259479 | 1.101 | 25281 | 38189 | 0.725 |
| 24 | planted | basis_wide | 53.390625 | 53.390625 | 4.894 | 25281 | 38189 | 0.577 |
| 24 | planted | tail_list | 125.170687 | 125.170687 | 2.088 | 15496 | 23398 | 0.614 |
| 24 | planted | tail_wide | 40.682292 | 40.682292 | 6.423 | 15496 | 23398 | 0.508 |
| 24 | planted | packed_untraced | 12.401042 | 12.401042 | 21.071 | 19200 | 0 | 0.000 |
| 24 | planted | affine_sl_basis_list | 266.847854 | 266.847854 | 0.979 | 6319 | 15922 | 0.497 |
| 24 | planted | affine_sl_basis_fast | 30.022771 | 30.022771 | 8.703 | 6319 | 15922 | 0.193 |
| 24 | planted | gray_scalar | 5.292791 | 5.292791 | 49.368 | 0 | 0 | 0.000 |
| 24 | planted | gray_simd | 0.802625 | 0.802625 | 325.552 | 0 | 0 | 0.000 |
| 24 | planted | packed_gray12_scalar | 12.976146 | 12.976146 | 20.137 | 3618 | 0 | 0.000 |
| 24 | planted | packed_gray12_simd | 5.622458 | 5.622458 | 46.474 | 3618 | 0 | 0.000 |
| 24 | planted | packed_gray16_scalar | 9.533771 | 9.533771 | 27.407 | 232 | 0 | 0.000 |
| 24 | planted | packed_gray16_simd | 1.691271 | 1.691271 | 154.497 | 232 | 0 | 0.000 |
| 24 | cross_planted | search | 80.110521 | 80.110521 | 1.955 | 20046 | 0 | 0.000 |
| 24 | cross_planted | flat | 156.635917 | 156.635917 | 1.000 | 1319 | 1338 | 0.966 |
| 24 | cross_planted | bucket | 539.046000 | 539.046000 | 0.291 | 1319 | 1338 | 0.989 |
| 24 | cross_planted | hybrid | 197.977687 | 197.977687 | 0.791 | 1319 | 1338 | 0.976 |
| 24 | cross_planted | small_flat | 90.386396 | 90.386396 | 1.733 | 11577 | 3542 | 0.386 |
| 24 | cross_planted | word_tail | 77.630188 | 77.630188 | 2.018 | 11577 | 3542 | 0.298 |
| 24 | cross_planted | merge_search | 70.610833 | 70.610833 | 2.218 | 20046 | 0 | 0.000 |
| 24 | cross_planted | quadratic_state | 28.916208 | 28.916208 | 5.417 | 20046 | 0 | 0.000 |
| 24 | cross_planted | packed_state | 17.828833 | 17.828833 | 8.786 | 20046 | 0 | 0.000 |
| 24 | cross_planted | basis_list | 48.409396 | 48.409396 | 3.236 | 5400 | 8142 | 0.727 |
| 24 | cross_planted | basis_wide | 11.743813 | 11.743813 | 13.338 | 5400 | 8142 | 0.579 |
| 24 | cross_planted | tail_list | 85.458624 | 85.458624 | 1.833 | 10764 | 16078 | 0.612 |
| 24 | cross_planted | tail_wide | 26.818541 | 26.818541 | 5.841 | 10764 | 16078 | 0.505 |
| 24 | cross_planted | packed_untraced | 12.482542 | 12.482542 | 12.548 | 20046 | 0 | 0.000 |
| 24 | cross_planted | affine_sl_basis_list | 60.621437 | 60.621437 | 2.584 | 1460 | 3773 | 0.489 |
| 24 | cross_planted | affine_sl_basis_fast | 7.068459 | 7.068459 | 22.160 | 1460 | 3773 | 0.185 |
| 24 | cross_planted | gray_scalar | 13.534583 | 13.534583 | 11.573 | 0 | 0 | 0.000 |
| 24 | cross_planted | gray_simd | 2.098875 | 2.098875 | 74.629 | 0 | 0 | 0.000 |
| 24 | cross_planted | packed_gray12_scalar | 12.147833 | 12.147833 | 12.894 | 3444 | 0 | 0.000 |
| 24 | cross_planted | packed_gray12_simd | 5.275979 | 5.275979 | 29.689 | 3444 | 0 | 0.000 |
| 24 | cross_planted | packed_gray16_scalar | 9.073272 | 9.073272 | 17.263 | 222 | 0 | 0.000 |
| 24 | cross_planted | packed_gray16_simd | 1.596583 | 1.596583 | 98.107 | 222 | 0 | 0.000 |
| 24 | unplanted | search | 143.931354 | 143.931354 | 3.545 | 35700 | 0 | 0.000 |
| 24 | unplanted | flat | 510.276729 | 510.276729 | 1.000 | 3363 | 3434 | 0.970 |
| 24 | unplanted | bucket | 1861.392605 | 1861.392605 | 0.274 | 3363 | 3434 | 0.991 |
| 24 | unplanted | hybrid | 826.881042 | 826.881042 | 0.617 | 3363 | 3434 | 0.981 |
| 24 | unplanted | small_flat | 165.173479 | 165.173479 | 3.089 | 24350 | 4986 | 0.315 |
| 24 | unplanted | word_tail | 145.692750 | 145.692750 | 3.502 | 24350 | 4986 | 0.241 |
| 24 | unplanted | merge_search | 126.581771 | 126.581771 | 4.031 | 35700 | 0 | 0.000 |
| 24 | unplanted | quadratic_state | 51.876855 | 51.876855 | 9.836 | 35700 | 0 | 0.000 |
| 24 | unplanted | packed_state | 32.550312 | 32.550312 | 15.677 | 35700 | 0 | 0.000 |
| 24 | unplanted | basis_list | 357.864375 | 357.864375 | 1.426 | 37628 | 55785 | 0.733 |
| 24 | unplanted | basis_wide | 82.236834 | 82.236834 | 6.205 | 37628 | 55785 | 0.587 |
| 24 | unplanted | tail_list | 229.339000 | 229.339000 | 2.225 | 26174 | 38643 | 0.624 |
| 24 | unplanted | tail_wide | 70.524333 | 70.524333 | 7.235 | 26174 | 38643 | 0.513 |
| 24 | unplanted | packed_untraced | 22.925854 | 22.925854 | 22.258 | 35700 | 0 | 0.000 |
| 24 | unplanted | affine_sl_basis_list | 471.918396 | 471.918396 | 1.081 | 8526 | 23884 | 0.538 |
| 24 | unplanted | affine_sl_basis_fast | 42.978333 | 42.978333 | 11.873 | 8526 | 23884 | 0.212 |
| 24 | unplanted | gray_scalar | 20.737104 | 20.737104 | 24.607 | 0 | 0 | 0.000 |
| 24 | unplanted | gray_simd | 3.210771 | 3.210771 | 158.927 | 0 | 0 | 0.000 |
| 24 | unplanted | packed_gray12_scalar | 28.111251 | 28.111251 | 18.152 | 7978 | 0 | 0.000 |
| 24 | unplanted | packed_gray12_simd | 12.244438 | 12.244438 | 41.674 | 7978 | 0 | 0.000 |
| 24 | unplanted | packed_gray16_scalar | 21.378792 | 21.378792 | 23.868 | 511 | 0 | 0.000 |
| 24 | unplanted | packed_gray16_simd | 3.761563 | 3.761563 | 135.656 | 511 | 0 | 0.000 |

Kernel-backed solvers have equal outcomes/models, logical counters and trace digests within each matched policy group. The search-only control may explore a different tree. SAT models are checked against original equations; UNSAT verification requires a completed independent search reference. UNKNOWN is censored, not UNSAT.

This driver solves bounded generated Boolean systems. It does not benchmark the repository inherited-F4 implementation, accept curve targets, recover scalars or establish index-calculus performance. Production and cryptanalytic costs stay null.

Selective one-word degree-2 full-solve gate: **REJECTED**. The small_flat and word_tail methods share a policy and must match exact traces/counters. They are compared separately with the always-degree-3 methods and search-only.

Ordered-specialization complete-solve gate: **REJECTED**. The search and merge_search arms have identical outcomes/models, logical work and trace digests. Only representation work in specialization changes. The gate compares merge_search with the fastest of all six prior full-solve methods.

Fixed-quadratic complete-solve gate: **REJECTED**. Compilation is charged to each cold solve. Search, merge_search and quadratic_state must match outcomes/models, logical counters and trace digests. The reference is the pointwise fastest of all seven prior methods.

| Variables | Family | Fastest prior / quadratic | 95% paired interval | Gate |
|---:|---|---:|---|---|
| 16 | planted | 2.465 | 2.446–2.475 | PASS |
| 16 | cross_planted | 1.722 | 1.715–1.734 | REJECTED |
| 16 | unplanted | 2.487 | 2.479–2.501 | PASS |
| 20 | planted | 2.419 | 2.406–2.427 | PASS |
| 20 | cross_planted | 2.435 | 2.412–2.451 | PASS |
| 20 | unplanted | 2.537 | 2.521–2.565 | PASS |
| 24 | planted | 2.443 | 2.417–2.462 | PASS |
| 24 | cross_planted | 2.427 | 2.375–2.440 | PASS |
| 24 | unplanted | 2.459 | 2.434–2.485 | PASS |

These paired intervals describe repeated timings on two fixed holdout fixtures per size/family, not a confidence interval over a population of Boolean systems. Fresh-worker RSS includes all arms and is not candidate-specific memory.

Packed-state current-frontier 2x gate: **REJECTED**. Incremental improvement gate (>1x lower bound): **PASS**. Cumulative historical-seven 2x gate: **PASS**.

Both fresh holdouts and the entire prior confirmation regression grid must pass. The current frontier includes the already faster quadratic_state method; the historical-seven comparison cannot substitute for it.

| Split | Variables | Family | Reference | Paired ratio | 95% interval | >2x gate |
|---|---:|---|---|---:|---|---|
| regression | 16 | planted | current_frontier | 1.592 | 1.584–1.601 | REJECTED |
| regression | 16 | planted | historical_seven | 3.823 | 3.776–3.876 | PASS |
| regression | 16 | cross_planted | current_frontier | 1.602 | 1.591–1.609 | REJECTED |
| regression | 16 | cross_planted | historical_seven | 3.467 | 3.329–3.551 | PASS |
| regression | 16 | unplanted | current_frontier | 1.582 | 1.577–1.588 | REJECTED |
| regression | 16 | unplanted | historical_seven | 3.923 | 3.908–3.939 | PASS |
| regression | 20 | planted | current_frontier | 1.604 | 1.598–1.613 | REJECTED |
| regression | 20 | planted | historical_seven | 3.911 | 3.900–3.924 | PASS |
| regression | 20 | cross_planted | current_frontier | 1.621 | 1.615–1.631 | REJECTED |
| regression | 20 | cross_planted | historical_seven | 3.875 | 3.836–3.919 | PASS |
| regression | 20 | unplanted | current_frontier | 1.582 | 1.569–1.593 | REJECTED |
| regression | 20 | unplanted | historical_seven | 3.935 | 3.915–3.957 | PASS |
| regression | 24 | planted | current_frontier | 1.624 | 1.613–1.636 | REJECTED |
| regression | 24 | planted | historical_seven | 3.911 | 3.897–3.923 | PASS |
| regression | 24 | cross_planted | current_frontier | 1.631 | 1.625–1.638 | REJECTED |
| regression | 24 | cross_planted | historical_seven | 3.915 | 3.895–3.926 | PASS |
| regression | 24 | unplanted | current_frontier | 1.630 | 1.628–1.633 | REJECTED |
| regression | 24 | unplanted | historical_seven | 3.937 | 3.922–3.947 | PASS |
| holdout | 16 | planted | current_frontier | 1.578 | 1.562–1.588 | REJECTED |
| holdout | 16 | planted | historical_seven | 3.896 | 3.857–3.909 | PASS |
| holdout | 16 | cross_planted | current_frontier | 1.574 | 1.564–1.588 | REJECTED |
| holdout | 16 | cross_planted | historical_seven | 2.725 | 2.712–2.746 | PASS |
| holdout | 16 | unplanted | current_frontier | 1.580 | 1.568–1.599 | REJECTED |
| holdout | 16 | unplanted | historical_seven | 3.932 | 3.893–3.964 | PASS |
| holdout | 20 | planted | current_frontier | 1.592 | 1.576–1.598 | REJECTED |
| holdout | 20 | planted | historical_seven | 3.849 | 3.823–3.861 | PASS |
| holdout | 20 | cross_planted | current_frontier | 1.578 | 1.566–1.605 | REJECTED |
| holdout | 20 | cross_planted | historical_seven | 3.844 | 3.828–3.858 | PASS |
| holdout | 20 | unplanted | current_frontier | 1.573 | 1.556–1.584 | REJECTED |
| holdout | 20 | unplanted | historical_seven | 3.971 | 3.940–4.008 | PASS |
| holdout | 24 | planted | current_frontier | 1.616 | 1.605–1.634 | REJECTED |
| holdout | 24 | planted | historical_seven | 3.959 | 3.944–3.969 | PASS |
| holdout | 24 | cross_planted | current_frontier | 1.637 | 1.628–1.688 | REJECTED |
| holdout | 24 | cross_planted | historical_seven | 3.977 | 3.958–3.986 | PASS |
| holdout | 24 | unplanted | current_frontier | 1.592 | 1.573–1.603 | REJECTED |
| holdout | 24 | unplanted | historical_seven | 3.906 | 3.902–3.909 | PASS |

## Transported coefficient-row-space comparison

basis_list/basis_wide share a full-RREF policy and branch on that basis. tail_list/tail_wide keep quadratic rows in echelon form, reduce only the affine tail, and branch using carried original-equation residuals. The backends in each pair must match every logical counter, model and trace. New inference policies need not match the older search tree.

Source-row and distinct-column counts price each reduction input, including zero rows. Specialization counts cover every carried representation. Basis histograms use nominal unassigned-variable width. All costs, including the extra original-equation state in the tail policy, are inside cold solve timing.

| Candidate | Dramatic >2x strongest-reference gate | Any improvement strongest-reference gate | Matched backend >1.05x gate |
|---|---|---|---|
| basis_wide | REJECTED | REJECTED | PASS |
| tail_wide | REJECTED | REJECTED | PASS |

| Candidate | Split | Variables | Family | Reference | Paired ratio | 95% interval | Gate |
|---|---|---:|---|---|---:|---|---|
| basis_wide | regression | 16 | planted | strongest | 0.007 | 0.006–0.007 | REJECTED |
| basis_wide | regression | 16 | planted | same_policy | 5.081 | 5.057–5.120 | PASS |
| basis_wide | regression | 16 | cross_planted | strongest | 0.022 | 0.009–0.038 | REJECTED |
| basis_wide | regression | 16 | cross_planted | same_policy | 4.829 | 4.735–4.914 | PASS |
| basis_wide | regression | 16 | unplanted | strongest | 0.009 | 0.008–0.009 | REJECTED |
| basis_wide | regression | 16 | unplanted | same_policy | 5.119 | 5.081–5.158 | PASS |
| basis_wide | regression | 20 | planted | strongest | 0.018 | 0.016–0.019 | REJECTED |
| basis_wide | regression | 20 | planted | same_policy | 4.706 | 4.683–4.752 | PASS |
| basis_wide | regression | 20 | cross_planted | strongest | 0.053 | 0.034–0.067 | REJECTED |
| basis_wide | regression | 20 | cross_planted | same_policy | 4.660 | 4.635–4.701 | PASS |
| basis_wide | regression | 20 | unplanted | strongest | 0.018 | 0.017–0.018 | REJECTED |
| basis_wide | regression | 20 | unplanted | same_policy | 4.700 | 4.644–4.855 | PASS |
| basis_wide | regression | 24 | planted | strongest | 0.027 | 0.021–0.031 | REJECTED |
| basis_wide | regression | 24 | planted | same_policy | 4.434 | 4.385–4.481 | PASS |
| basis_wide | regression | 24 | cross_planted | strongest | 0.057 | 0.036–0.071 | REJECTED |
| basis_wide | regression | 24 | cross_planted | same_policy | 4.437 | 4.413–4.477 | PASS |
| basis_wide | regression | 24 | unplanted | strongest | 0.029 | 0.028–0.030 | REJECTED |
| basis_wide | regression | 24 | unplanted | same_policy | 4.432 | 4.388–4.456 | PASS |
| basis_wide | holdout | 16 | planted | strongest | 0.008 | 0.007–0.009 | REJECTED |
| basis_wide | holdout | 16 | planted | same_policy | 5.262 | 5.223–5.328 | PASS |
| basis_wide | holdout | 16 | cross_planted | strongest | 0.019 | 0.018–0.021 | REJECTED |
| basis_wide | holdout | 16 | cross_planted | same_policy | 5.273 | 5.201–5.348 | PASS |
| basis_wide | holdout | 16 | unplanted | strongest | 0.008 | 0.007–0.008 | REJECTED |
| basis_wide | holdout | 16 | unplanted | same_policy | 5.276 | 5.233–5.321 | PASS |
| basis_wide | holdout | 20 | planted | strongest | 0.013 | 0.007–0.019 | REJECTED |
| basis_wide | holdout | 20 | planted | same_policy | 4.762 | 4.641–4.854 | PASS |
| basis_wide | holdout | 20 | cross_planted | strongest | 0.015 | 0.014–0.015 | REJECTED |
| basis_wide | holdout | 20 | cross_planted | same_policy | 4.769 | 4.623–4.918 | PASS |
| basis_wide | holdout | 20 | unplanted | strongest | 0.013 | 0.013–0.014 | REJECTED |
| basis_wide | holdout | 20 | unplanted | same_policy | 4.878 | 4.843–4.908 | PASS |
| basis_wide | holdout | 24 | planted | strongest | 0.016 | 0.016–0.016 | REJECTED |
| basis_wide | holdout | 24 | planted | same_policy | 4.450 | 4.439–4.460 | PASS |
| basis_wide | holdout | 24 | cross_planted | strongest | 0.129 | 0.061–0.195 | REJECTED |
| basis_wide | holdout | 24 | cross_planted | same_policy | 4.119 | 4.100–4.128 | PASS |
| basis_wide | holdout | 24 | unplanted | strongest | 0.041 | 0.041–0.042 | REJECTED |
| basis_wide | holdout | 24 | unplanted | same_policy | 4.349 | 4.281–4.411 | PASS |
| tail_wide | regression | 16 | planted | strongest | 0.013 | 0.009–0.017 | REJECTED |
| tail_wide | regression | 16 | planted | same_policy | 4.214 | 4.171–4.250 | PASS |
| tail_wide | regression | 16 | cross_planted | strongest | 0.028 | 0.011–0.051 | REJECTED |
| tail_wide | regression | 16 | cross_planted | same_policy | 4.200 | 4.148–4.234 | PASS |
| tail_wide | regression | 16 | unplanted | strongest | 0.010 | 0.009–0.010 | REJECTED |
| tail_wide | regression | 16 | unplanted | same_policy | 4.209 | 4.118–4.286 | PASS |
| tail_wide | regression | 20 | planted | strongest | 0.011 | 0.005–0.018 | REJECTED |
| tail_wide | regression | 20 | planted | same_policy | 3.557 | 3.463–3.607 | PASS |
| tail_wide | regression | 20 | cross_planted | strongest | 0.018 | 0.006–0.034 | REJECTED |
| tail_wide | regression | 20 | cross_planted | same_policy | 3.421 | 3.403–3.452 | PASS |
| tail_wide | regression | 20 | unplanted | strongest | 0.020 | 0.019–0.021 | REJECTED |
| tail_wide | regression | 20 | unplanted | same_policy | 3.484 | 3.440–3.519 | PASS |
| tail_wide | regression | 24 | planted | strongest | 0.045 | 0.040–0.046 | REJECTED |
| tail_wide | regression | 24 | planted | same_policy | 3.224 | 3.181–3.269 | PASS |
| tail_wide | regression | 24 | cross_planted | strongest | 0.074 | 0.073–0.079 | REJECTED |
| tail_wide | regression | 24 | cross_planted | same_policy | 3.209 | 3.129–3.241 | PASS |
| tail_wide | regression | 24 | unplanted | strongest | 0.042 | 0.040–0.044 | REJECTED |
| tail_wide | regression | 24 | unplanted | same_policy | 3.199 | 3.166–3.257 | PASS |
| tail_wide | holdout | 16 | planted | strongest | 0.007 | 0.002–0.012 | REJECTED |
| tail_wide | holdout | 16 | planted | same_policy | 4.174 | 4.105–4.264 | PASS |
| tail_wide | holdout | 16 | cross_planted | strongest | 0.016 | 0.005–0.030 | REJECTED |
| tail_wide | holdout | 16 | cross_planted | same_policy | 3.956 | 3.892–4.042 | PASS |
| tail_wide | holdout | 16 | unplanted | strongest | 0.010 | 0.007–0.012 | REJECTED |
| tail_wide | holdout | 16 | unplanted | same_policy | 4.023 | 3.906–4.272 | PASS |
| tail_wide | holdout | 20 | planted | strongest | 0.019 | 0.013–0.024 | REJECTED |
| tail_wide | holdout | 20 | planted | same_policy | 3.473 | 3.434–3.519 | PASS |
| tail_wide | holdout | 20 | cross_planted | strongest | 0.033 | 0.023–0.041 | REJECTED |
| tail_wide | holdout | 20 | cross_planted | same_policy | 3.549 | 3.519–3.614 | PASS |
| tail_wide | holdout | 20 | unplanted | strongest | 0.017 | 0.016–0.018 | REJECTED |
| tail_wide | holdout | 20 | unplanted | same_policy | 3.539 | 3.425–3.671 | PASS |
| tail_wide | holdout | 24 | planted | strongest | 0.041 | 0.019–0.064 | REJECTED |
| tail_wide | holdout | 24 | planted | same_policy | 3.090 | 3.078–3.105 | PASS |
| tail_wide | holdout | 24 | cross_planted | strongest | 0.055 | 0.027–0.083 | REJECTED |
| tail_wide | holdout | 24 | cross_planted | same_policy | 3.125 | 3.065–3.173 | PASS |
| tail_wide | holdout | 24 | unplanted | strongest | 0.050 | 0.041–0.061 | REJECTED |
| tail_wide | holdout | 24 | unplanted | same_policy | 3.230 | 3.124–3.372 | PASS |

## Frozen finalist comparison

Reference arms include all thirteen previously retained methods, the packed method without diagnostic hashing, the affine finalist pair, and all scalar enumeration controls. The three SIMD treatments are compared with this fixed reference roster. Comparisons among the SIMD treatments do not redefine that preregistered baseline. The affine finalist is also shown; its matched list implementation is a backend control.

Enumeration uses assignment and block counters, not search-node counts. Its checksum folds each block index and the wrapping sum of all sixteen equation-syndrome words; it is a diagnostic, not a certificate or a canonical search-tree trace. Scalar and SIMD modes must match models, all logical counters and checksums. Tree hybrids also preserve their prefix traces. The untraced packed control must preserve models and logical counters; its absent trace is null.

| Candidate | >2x retained-frontier gate | >1.05x matched-backend gate |
|---|---|---|
| affine_sl_basis_fast | REJECTED | PASS |
| gray_simd | REJECTED | PASS |
| packed_gray12_simd | REJECTED | PASS |
| packed_gray16_simd | REJECTED | PASS |

| Candidate | Split | Variables | Family | Reference | Paired ratio | 95% interval | Gate |
|---|---|---:|---|---|---:|---|---|
| affine_sl_basis_fast | regression | 16 | planted | retained_frontier | 0.032 | 0.025–0.036 | REJECTED |
| affine_sl_basis_fast | regression | 16 | planted | same_policy | 11.740 | 11.374–12.294 | PASS |
| affine_sl_basis_fast | regression | 16 | cross_planted | retained_frontier | 0.075 | 0.054–0.096 | REJECTED |
| affine_sl_basis_fast | regression | 16 | cross_planted | same_policy | 10.383 | 10.288–10.460 | PASS |
| affine_sl_basis_fast | regression | 16 | unplanted | retained_frontier | 0.070 | 0.067–0.075 | REJECTED |
| affine_sl_basis_fast | regression | 16 | unplanted | same_policy | 12.409 | 12.336–12.517 | PASS |
| affine_sl_basis_fast | regression | 20 | planted | retained_frontier | 0.101 | 0.051–0.155 | REJECTED |
| affine_sl_basis_fast | regression | 20 | planted | same_policy | 10.060 | 9.340–10.246 | PASS |
| affine_sl_basis_fast | regression | 20 | cross_planted | retained_frontier | 0.267 | 0.038–0.539 | REJECTED |
| affine_sl_basis_fast | regression | 20 | cross_planted | same_policy | 9.673 | 9.221–10.275 | PASS |
| affine_sl_basis_fast | regression | 20 | unplanted | retained_frontier | 0.161 | 0.160–0.164 | REJECTED |
| affine_sl_basis_fast | regression | 20 | unplanted | same_policy | 10.595 | 10.435–10.739 | PASS |
| affine_sl_basis_fast | regression | 24 | planted | retained_frontier | 0.334 | 0.270–0.424 | REJECTED |
| affine_sl_basis_fast | regression | 24 | planted | same_policy | 9.729 | 9.567–9.863 | PASS |
| affine_sl_basis_fast | regression | 24 | cross_planted | retained_frontier | 0.368 | 0.290–0.460 | REJECTED |
| affine_sl_basis_fast | regression | 24 | cross_planted | same_policy | 8.600 | 8.256–8.963 | PASS |
| affine_sl_basis_fast | regression | 24 | unplanted | retained_frontier | 0.371 | 0.368–0.376 | REJECTED |
| affine_sl_basis_fast | regression | 24 | unplanted | same_policy | 9.560 | 9.268–9.760 | PASS |
| affine_sl_basis_fast | holdout | 16 | planted | retained_frontier | 0.046 | 0.015–0.083 | REJECTED |
| affine_sl_basis_fast | holdout | 16 | planted | same_policy | 11.836 | 10.235–13.589 | PASS |
| affine_sl_basis_fast | holdout | 16 | cross_planted | retained_frontier | 0.120 | 0.039–0.208 | REJECTED |
| affine_sl_basis_fast | holdout | 16 | cross_planted | same_policy | 9.979 | 8.840–11.485 | PASS |
| affine_sl_basis_fast | holdout | 16 | unplanted | retained_frontier | 0.069 | 0.046–0.094 | REJECTED |
| affine_sl_basis_fast | holdout | 16 | unplanted | same_policy | 12.850 | 12.479–13.287 | PASS |
| affine_sl_basis_fast | holdout | 20 | planted | retained_frontier | 0.310 | 0.075–0.673 | REJECTED |
| affine_sl_basis_fast | holdout | 20 | planted | same_policy | 9.758 | 8.611–10.850 | PASS |
| affine_sl_basis_fast | holdout | 20 | cross_planted | retained_frontier | 0.138 | 0.136–0.140 | REJECTED |
| affine_sl_basis_fast | holdout | 20 | cross_planted | same_policy | 10.719 | 10.645–10.767 | PASS |
| affine_sl_basis_fast | holdout | 20 | unplanted | retained_frontier | 0.161 | 0.157–0.164 | REJECTED |
| affine_sl_basis_fast | holdout | 20 | unplanted | same_policy | 11.268 | 10.790–11.514 | PASS |
| affine_sl_basis_fast | holdout | 24 | planted | retained_frontier | 0.252 | 0.144–0.338 | REJECTED |
| affine_sl_basis_fast | holdout | 24 | planted | same_policy | 8.256 | 8.152–8.903 | PASS |
| affine_sl_basis_fast | holdout | 24 | cross_planted | retained_frontier | 0.643 | 0.626–0.680 | REJECTED |
| affine_sl_basis_fast | holdout | 24 | cross_planted | same_policy | 7.674 | 6.627–8.695 | PASS |
| affine_sl_basis_fast | holdout | 24 | unplanted | retained_frontier | 0.429 | 0.428–0.431 | REJECTED |
| affine_sl_basis_fast | holdout | 24 | unplanted | same_policy | 11.009 | 10.187–11.981 | PASS |
| gray_simd | regression | 16 | planted | retained_frontier | 3.265 | 2.767–4.037 | PASS |
| gray_simd | regression | 16 | planted | same_policy | 4.845 | 4.637–5.127 | PASS |
| gray_simd | regression | 16 | cross_planted | retained_frontier | 3.416 | 2.996–3.986 | PASS |
| gray_simd | regression | 16 | cross_planted | same_policy | 4.864 | 4.685–5.113 | PASS |
| gray_simd | regression | 16 | unplanted | retained_frontier | 5.452 | 5.329–5.792 | PASS |
| gray_simd | regression | 16 | unplanted | same_policy | 5.760 | 5.459–6.092 | PASS |
| gray_simd | regression | 20 | planted | retained_frontier | 5.427 | 5.312–5.563 | PASS |
| gray_simd | regression | 20 | planted | same_policy | 5.549 | 5.485–5.663 | PASS |
| gray_simd | regression | 20 | cross_planted | retained_frontier | 5.350 | 5.265–5.458 | PASS |
| gray_simd | regression | 20 | cross_planted | same_policy | 5.486 | 5.445–5.541 | PASS |
| gray_simd | regression | 20 | unplanted | retained_frontier | 5.839 | 5.566–6.435 | PASS |
| gray_simd | regression | 20 | unplanted | same_policy | 6.379 | 5.817–6.534 | PASS |
| gray_simd | regression | 24 | planted | retained_frontier | 5.477 | 5.360–5.527 | PASS |
| gray_simd | regression | 24 | planted | same_policy | 6.541 | 5.949–6.686 | PASS |
| gray_simd | regression | 24 | cross_planted | retained_frontier | 2.364 | 1.978–3.319 | REJECTED |
| gray_simd | regression | 24 | cross_planted | same_policy | 6.593 | 6.509–6.700 | PASS |
| gray_simd | regression | 24 | unplanted | retained_frontier | 6.035 | 5.803–6.524 | PASS |
| gray_simd | regression | 24 | unplanted | same_policy | 6.686 | 6.617–6.803 | PASS |
| gray_simd | holdout | 16 | planted | retained_frontier | 4.141 | 4.026–4.375 | PASS |
| gray_simd | holdout | 16 | planted | same_policy | 4.192 | 4.074–4.408 | PASS |
| gray_simd | holdout | 16 | cross_planted | retained_frontier | 4.286 | 4.210–4.358 | PASS |
| gray_simd | holdout | 16 | cross_planted | same_policy | 4.286 | 4.210–4.358 | PASS |
| gray_simd | holdout | 16 | unplanted | retained_frontier | 5.450 | 5.246–5.704 | PASS |
| gray_simd | holdout | 16 | unplanted | same_policy | 5.492 | 5.255–5.707 | PASS |
| gray_simd | holdout | 20 | planted | retained_frontier | 4.359 | 3.417–5.866 | PASS |
| gray_simd | holdout | 20 | planted | same_policy | 6.425 | 5.866–6.583 | PASS |
| gray_simd | holdout | 20 | cross_planted | retained_frontier | 4.259 | 3.363–5.473 | PASS |
| gray_simd | holdout | 20 | cross_planted | same_policy | 5.613 | 5.506–6.326 | PASS |
| gray_simd | holdout | 20 | unplanted | retained_frontier | 6.599 | 6.249–6.744 | PASS |
| gray_simd | holdout | 20 | unplanted | same_policy | 6.599 | 6.331–6.752 | PASS |
| gray_simd | holdout | 24 | planted | retained_frontier | 4.158 | 2.723–5.550 | PASS |
| gray_simd | holdout | 24 | planted | same_policy | 6.576 | 5.625–6.806 | PASS |
| gray_simd | holdout | 24 | cross_planted | retained_frontier | 2.873 | 0.236–5.913 | REJECTED |
| gray_simd | holdout | 24 | cross_planted | same_policy | 6.783 | 6.431–6.854 | PASS |
| gray_simd | holdout | 24 | unplanted | retained_frontier | 5.418 | 5.158–5.638 | PASS |
| gray_simd | holdout | 24 | unplanted | same_policy | 6.435 | 5.936–6.766 | PASS |
| packed_gray12_simd | regression | 16 | planted | retained_frontier | 1.233 | 1.126–1.320 | REJECTED |
| packed_gray12_simd | regression | 16 | planted | same_policy | 2.132 | 2.068–2.172 | PASS |
| packed_gray12_simd | regression | 16 | cross_planted | retained_frontier | 1.218 | 1.106–1.353 | REJECTED |
| packed_gray12_simd | regression | 16 | cross_planted | same_policy | 2.123 | 2.074–2.163 | PASS |
| packed_gray12_simd | regression | 16 | unplanted | retained_frontier | 1.610 | 1.599–1.623 | REJECTED |
| packed_gray12_simd | regression | 16 | unplanted | same_policy | 2.319 | 2.295–2.359 | PASS |
| packed_gray12_simd | regression | 20 | planted | retained_frontier | 0.950 | 0.313–1.664 | REJECTED |
| packed_gray12_simd | regression | 20 | planted | same_policy | 2.377 | 2.357–2.391 | PASS |
| packed_gray12_simd | regression | 20 | cross_planted | retained_frontier | 0.873 | 0.301–1.535 | REJECTED |
| packed_gray12_simd | regression | 20 | cross_planted | same_policy | 2.404 | 2.393–2.417 | PASS |
| packed_gray12_simd | regression | 20 | unplanted | retained_frontier | 1.656 | 1.650–1.669 | REJECTED |
| packed_gray12_simd | regression | 20 | unplanted | same_policy | 2.395 | 2.381–2.404 | PASS |
| packed_gray12_simd | regression | 24 | planted | retained_frontier | 1.516 | 1.422–1.568 | REJECTED |
| packed_gray12_simd | regression | 24 | planted | same_policy | 2.245 | 2.215–2.272 | PASS |
| packed_gray12_simd | regression | 24 | cross_planted | retained_frontier | 1.363 | 1.324–1.477 | REJECTED |
| packed_gray12_simd | regression | 24 | cross_planted | same_policy | 2.216 | 2.203–2.241 | PASS |
| packed_gray12_simd | regression | 24 | unplanted | retained_frontier | 1.644 | 1.630–1.660 | REJECTED |
| packed_gray12_simd | regression | 24 | unplanted | same_policy | 2.241 | 2.233–2.257 | PASS |
| packed_gray12_simd | holdout | 16 | planted | retained_frontier | 0.646 | 0.303–1.264 | REJECTED |
| packed_gray12_simd | holdout | 16 | planted | same_policy | 2.187 | 2.090–2.243 | PASS |
| packed_gray12_simd | holdout | 16 | cross_planted | retained_frontier | 0.425 | 0.298–0.551 | REJECTED |
| packed_gray12_simd | holdout | 16 | cross_planted | same_policy | 2.234 | 2.205–2.267 | PASS |
| packed_gray12_simd | holdout | 16 | unplanted | retained_frontier | 1.542 | 1.513–1.586 | REJECTED |
| packed_gray12_simd | holdout | 16 | unplanted | same_policy | 2.335 | 2.272–2.374 | PASS |
| packed_gray12_simd | holdout | 20 | planted | retained_frontier | 1.296 | 1.043–1.778 | REJECTED |
| packed_gray12_simd | holdout | 20 | planted | same_policy | 2.388 | 2.340–2.406 | PASS |
| packed_gray12_simd | holdout | 20 | cross_planted | retained_frontier | 1.430 | 1.267–1.741 | REJECTED |
| packed_gray12_simd | holdout | 20 | cross_planted | same_policy | 2.364 | 2.343–2.380 | PASS |
| packed_gray12_simd | holdout | 20 | unplanted | retained_frontier | 1.675 | 1.663–1.705 | REJECTED |
| packed_gray12_simd | holdout | 20 | unplanted | same_policy | 2.431 | 2.409–2.446 | PASS |
| packed_gray12_simd | holdout | 24 | planted | retained_frontier | 1.276 | 0.777–1.460 | REJECTED |
| packed_gray12_simd | holdout | 24 | planted | same_policy | 2.258 | 2.218–2.295 | PASS |
| packed_gray12_simd | holdout | 24 | cross_planted | retained_frontier | 1.063 | 0.816–1.426 | REJECTED |
| packed_gray12_simd | holdout | 24 | cross_planted | same_policy | 2.168 | 2.086–2.295 | PASS |
| packed_gray12_simd | holdout | 24 | unplanted | retained_frontier | 1.493 | 1.336–1.661 | REJECTED |
| packed_gray12_simd | holdout | 24 | unplanted | same_policy | 2.277 | 2.255–2.327 | PASS |
| packed_gray16_simd | regression | 16 | planted | retained_frontier | 1.968 | 1.796–2.287 | REJECTED |
| packed_gray16_simd | regression | 16 | planted | same_policy | 3.687 | 3.473–3.853 | PASS |
| packed_gray16_simd | regression | 16 | cross_planted | retained_frontier | 2.182 | 2.052–2.340 | PASS |
| packed_gray16_simd | regression | 16 | cross_planted | same_policy | 3.573 | 3.360–3.749 | PASS |
| packed_gray16_simd | regression | 16 | unplanted | retained_frontier | 4.610 | 4.545–4.659 | PASS |
| packed_gray16_simd | regression | 16 | unplanted | same_policy | 4.747 | 4.670–4.811 | PASS |
| packed_gray16_simd | regression | 20 | planted | retained_frontier | 2.652 | 0.966–4.795 | REJECTED |
| packed_gray16_simd | regression | 20 | planted | same_policy | 5.435 | 5.266–5.487 | PASS |
| packed_gray16_simd | regression | 20 | cross_planted | retained_frontier | 2.610 | 0.904–4.724 | REJECTED |
| packed_gray16_simd | regression | 20 | cross_planted | same_policy | 5.422 | 5.361–5.493 | PASS |
| packed_gray16_simd | regression | 20 | unplanted | retained_frontier | 5.342 | 5.304–5.407 | PASS |
| packed_gray16_simd | regression | 20 | unplanted | same_policy | 5.540 | 5.492–5.590 | PASS |
| packed_gray16_simd | regression | 24 | planted | retained_frontier | 4.653 | 4.536–4.940 | PASS |
| packed_gray16_simd | regression | 24 | planted | same_policy | 5.587 | 5.552–5.636 | PASS |
| packed_gray16_simd | regression | 24 | cross_planted | retained_frontier | 3.818 | 3.468–4.556 | PASS |
| packed_gray16_simd | regression | 24 | cross_planted | same_policy | 5.584 | 5.558–5.599 | PASS |
| packed_gray16_simd | regression | 24 | unplanted | retained_frontier | 5.220 | 5.068–5.389 | PASS |
| packed_gray16_simd | regression | 24 | unplanted | same_policy | 5.593 | 5.555–5.639 | PASS |
| packed_gray16_simd | holdout | 16 | planted | retained_frontier | 2.484 | 2.326–2.750 | PASS |
| packed_gray16_simd | holdout | 16 | planted | same_policy | 2.789 | 2.601–3.039 | PASS |
| packed_gray16_simd | holdout | 16 | cross_planted | retained_frontier | 2.613 | 2.409–2.824 | PASS |
| packed_gray16_simd | holdout | 16 | cross_planted | same_policy | 2.819 | 2.690–3.056 | PASS |
| packed_gray16_simd | holdout | 16 | unplanted | retained_frontier | 4.602 | 4.387–4.944 | PASS |
| packed_gray16_simd | holdout | 16 | unplanted | same_policy | 4.706 | 4.525–5.055 | PASS |
| packed_gray16_simd | holdout | 20 | planted | retained_frontier | 4.028 | 3.411–5.146 | PASS |
| packed_gray16_simd | holdout | 20 | planted | same_policy | 5.394 | 5.253–5.463 | PASS |
| packed_gray16_simd | holdout | 20 | cross_planted | retained_frontier | 4.456 | 4.351–4.658 | PASS |
| packed_gray16_simd | holdout | 20 | cross_planted | same_policy | 5.380 | 5.285–5.481 | PASS |
| packed_gray16_simd | holdout | 20 | unplanted | retained_frontier | 5.394 | 5.333–5.477 | PASS |
| packed_gray16_simd | holdout | 20 | unplanted | same_policy | 5.604 | 5.541–5.667 | PASS |
| packed_gray16_simd | holdout | 24 | planted | retained_frontier | 3.671 | 2.634–4.125 | PASS |
| packed_gray16_simd | holdout | 24 | planted | same_policy | 5.549 | 5.500–5.606 | PASS |
| packed_gray16_simd | holdout | 24 | cross_planted | retained_frontier | 2.853 | 2.714–3.135 | PASS |
| packed_gray16_simd | holdout | 24 | cross_planted | same_policy | 5.524 | 5.427–5.585 | PASS |
| packed_gray16_simd | holdout | 24 | unplanted | retained_frontier | 4.798 | 4.234–5.415 | PASS |
| packed_gray16_simd | holdout | 24 | unplanted | same_policy | 5.667 | 5.598–5.722 | PASS |
