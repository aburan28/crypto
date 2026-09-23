# Boolean construction plus canonical linear reduction

Correctness **PASS**: 144 cells, 22680 batch-arm samples, 98280 oracle-verified RREF outputs.

Dramatic combined-workload gate: **REJECTED**, 0 / 9 comparisons passed.

Cold batch8 milliseconds include setup, construction, complete forward/backward elimination, compaction, exact validation and destruction. Medians pool two holdout seeds and 30 balanced repetitions. Ratios below use pooled medians; acceptance uses paired ratios against the fastest prior arm, including streaming.

| Variables | Family | Variant | Cold batch (ms) | Sorted / arm | Construction (ms) | Reduction (ms) | Fused (ms) | Reduction word XORs |
|---:|---|---|---:|---:|---:|---:|---:|---:|
| 12 | quadratic | sorted_reduce | 0.399333 | 1.000 | 0.278106 | 0.106521 | null | 58214 |
| 12 | quadratic | ranked_reduce | 0.211458 | 1.888 | 0.093271 | 0.103541 | null | 58214 |
| 12 | quadratic | sparse_reduce | 0.214938 | 1.858 | 0.097770 | 0.102938 | null | 58214 |
| 12 | quadratic | stream_reduce | 0.273626 | 1.459 | null | null | 0.259916 | 57533 |
| 12 | quadratic | incidence_reduce | 0.369417 | 1.081 | 0.095938 | 0.258354 | null | 55667 |
| 12 | quadratic | dense_reduce | 0.187229 | 2.133 | 0.064770 | 0.109229 | null | 58214 |
| 12 | linear_drop | sorted_reduce | 0.416458 | 1.000 | 0.282938 | 0.122542 | null | 61075 |
| 12 | linear_drop | ranked_reduce | 0.226813 | 1.836 | 0.090708 | 0.119750 | null | 61075 |
| 12 | linear_drop | sparse_reduce | 0.229333 | 1.816 | 0.094751 | 0.119790 | null | 61075 |
| 12 | linear_drop | stream_reduce | 0.277167 | 1.503 | null | null | 0.262001 | 60556 |
| 12 | linear_drop | incidence_reduce | 0.392396 | 1.061 | 0.097458 | 0.280855 | null | 55254 |
| 12 | linear_drop | dense_reduce | 0.200459 | 2.078 | 0.062375 | 0.121332 | null | 61075 |
| 12 | restricted_cycle | sorted_reduce | 0.285354 | 1.000 | 0.198125 | 0.075459 | null | 26346 |
| 12 | restricted_cycle | ranked_reduce | 0.162292 | 1.758 | 0.074355 | 0.074437 | null | 26346 |
| 12 | restricted_cycle | sparse_reduce | 0.155958 | 1.830 | 0.068436 | 0.075166 | null | 26346 |
| 12 | restricted_cycle | stream_reduce | 0.203750 | 1.401 | null | null | 0.191563 | 27766 |
| 12 | restricted_cycle | incidence_reduce | 0.238063 | 1.199 | 0.067395 | 0.156790 | null | 21200 |
| 12 | restricted_cycle | dense_reduce | 0.143063 | 1.995 | 0.054791 | 0.075022 | null | 26346 |
| 20 | quadratic | sorted_reduce | 0.729521 | 1.000 | 0.478106 | 0.227334 | null | 183338 |
| 20 | quadratic | ranked_reduce | 0.442458 | 1.649 | 0.192313 | 0.220063 | null | 183338 |
| 20 | quadratic | sparse_reduce | 0.434959 | 1.677 | 0.189813 | 0.215688 | null | 183338 |
| 20 | quadratic | stream_reduce | 0.699937 | 1.042 | null | null | 0.670228 | 267149 |
| 20 | quadratic | incidence_reduce | 0.701167 | 1.040 | 0.182708 | 0.490125 | null | 146092 |
| 20 | linear_drop | sorted_reduce | 0.893084 | 1.000 | 0.536625 | 0.330041 | null | 219928 |
| 20 | linear_drop | ranked_reduce | 0.573480 | 1.557 | 0.217500 | 0.321709 | null | 219928 |
| 20 | linear_drop | sparse_reduce | 0.569750 | 1.568 | 0.219603 | 0.318853 | null | 219928 |
| 20 | linear_drop | stream_reduce | 0.821646 | 1.087 | null | null | 0.788103 | 297056 |
| 20 | linear_drop | incidence_reduce | 0.855916 | 1.043 | 0.212478 | 0.603604 | null | 150952 |
| 20 | restricted_cycle | sorted_reduce | 0.280229 | 1.000 | 0.204583 | 0.065188 | null | 22954 |
| 20 | restricted_cycle | ranked_reduce | 0.172229 | 1.627 | 0.091166 | 0.064354 | null | 22954 |
| 20 | restricted_cycle | sparse_reduce | 0.149562 | 1.874 | 0.070042 | 0.064021 | null | 22954 |
| 20 | restricted_cycle | stream_reduce | 0.219771 | 1.275 | null | null | 0.204250 | 63142 |
| 20 | restricted_cycle | incidence_reduce | 0.223875 | 1.252 | 0.070292 | 0.138625 | null | 16578 |
| 28 | quadratic | sorted_reduce | 1.082187 | 1.000 | 0.676125 | 0.374250 | null | 339362 |
| 28 | quadratic | ranked_reduce | 0.736146 | 1.470 | 0.311603 | 0.362187 | null | 339362 |
| 28 | quadratic | sparse_reduce | 0.669729 | 1.616 | 0.266186 | 0.346751 | null | 339362 |
| 28 | quadratic | stream_reduce | 1.183876 | 0.914 | null | null | 1.118478 | 752766 |
| 28 | quadratic | incidence_reduce | 1.140334 | 0.949 | 0.265979 | 0.808500 | null | 233794 |
| 28 | linear_drop | sorted_reduce | 1.535604 | 1.000 | 0.800416 | 0.684541 | null | 455489 |
| 28 | linear_drop | ranked_reduce | 1.118396 | 1.373 | 0.375980 | 0.667208 | null | 455489 |
| 28 | linear_drop | sparse_reduce | 1.044854 | 1.470 | 0.306833 | 0.661355 | null | 455489 |
| 28 | linear_drop | stream_reduce | 1.566229 | 0.980 | null | null | 1.490498 | 865342 |
| 28 | linear_drop | incidence_reduce | 1.397417 | 1.099 | 0.303916 | 1.015521 | null | 244502 |
| 28 | restricted_cycle | sorted_reduce | 0.277187 | 1.000 | 0.202875 | 0.063042 | null | 20202 |
| 28 | restricted_cycle | ranked_reduce | 0.196500 | 1.411 | 0.109208 | 0.062793 | null | 20202 |
| 28 | restricted_cycle | sparse_reduce | 0.164042 | 1.690 | 0.074728 | 0.062521 | null | 20202 |
| 28 | restricted_cycle | stream_reduce | 0.270521 | 1.025 | null | null | 0.245022 | 123408 |
| 28 | restricted_cycle | incidence_reduce | 0.240375 | 1.153 | 0.074374 | 0.136916 | null | 14756 |
| 36 | quadratic | sorted_reduce | 1.489562 | 1.000 | 0.881478 | 0.559102 | null | 546914 |
| 36 | quadratic | ranked_reduce | 1.168688 | 1.275 | 0.497563 | 0.547188 | null | 546914 |
| 36 | quadratic | sparse_reduce | 0.996125 | 1.495 | 0.336396 | 0.533188 | null | 546914 |
| 36 | quadratic | stream_reduce | 1.912479 | 0.779 | null | null | 1.791750 | 1823976 |
| 36 | quadratic | incidence_reduce | 1.634313 | 0.911 | 0.335938 | 1.172268 | null | 325718 |
| 36 | linear_drop | sorted_reduce | 2.441646 | 1.000 | 1.102812 | 1.240041 | null | 781299 |
| 36 | linear_drop | ranked_reduce | 2.074896 | 1.177 | 0.676458 | 1.222333 | null | 781299 |
| 36 | linear_drop | sparse_reduce | 1.771021 | 1.379 | 0.408042 | 1.185915 | null | 781299 |
| 36 | linear_drop | stream_reduce | 2.882667 | 0.847 | null | null | 2.719938 | 2044122 |
| 36 | linear_drop | incidence_reduce | 2.035062 | 1.200 | 0.414083 | 1.456624 | null | 332952 |
| 36 | restricted_cycle | sorted_reduce | 0.272854 | 1.000 | 0.199210 | 0.062895 | null | 20018 |
| 36 | restricted_cycle | ranked_reduce | 0.257354 | 1.060 | 0.149459 | 0.062751 | null | 20018 |
| 36 | restricted_cycle | sparse_reduce | 0.183812 | 1.484 | 0.077272 | 0.062061 | null | 20018 |
| 36 | restricted_cycle | stream_reduce | 0.343083 | 0.795 | null | null | 0.297332 | 243771 |
| 36 | restricted_cycle | incidence_reduce | 0.256021 | 1.066 | 0.078624 | 0.131749 | null | 13862 |

Fused construction/reduction phases are not separately observable and remain null, not zero. Source rows and logical row-XOR counts agree across all arms; physical word-XOR counts may differ with ambient coordinate width. These counters cover reduction only, not all calibrated operations.

Retained context bytes and end-of-insertion basis storage are recorded separately. Neither is whole-process peak memory. Worker RSS includes all arms and the common reference corpus. Fixture and independent oracle construction are outside arm timing and inside process receipts.

This completes the bounded construction-plus-RREF task. It does not compute complete Groebner closure, enumerate roots, solve the original polynomial system, or measure index-calculus performance. Those costs and rho ratios remain null.

This additive run evaluates incidence_reduce: sparse construction, cached nonzero pivot words and precomputed backward-elimination incidence. All preparation and auxiliary storage are charged; exact canonical output and logical row-XOR counts still match. The primary gate compares it against the pointwise fastest prior pipeline, including the streaming control.
