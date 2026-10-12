# Boolean construction plus canonical linear reduction

Correctness **PASS**: 144 cells, 12240 batch-arm samples, 53040 oracle-verified RREF outputs.

Dramatic combined-workload gate: **REJECTED**, 0 / 9 comparisons passed.

Cold batch8 milliseconds include setup, construction, complete forward/backward elimination, compaction, exact validation and destruction. Medians pool two holdout seeds and 20 balanced repetitions. Ratios below use pooled medians; acceptance uses paired ratios against the fastest staged arm.

| Variables | Family | Variant | Cold batch (ms) | Sorted / arm | Construction (ms) | Reduction (ms) | Fused (ms) | Reduction word XORs |
|---:|---|---|---:|---:|---:|---:|---:|---:|
| 12 | quadratic | sorted_reduce | 0.377979 | 1.000 | 0.269521 | 0.097939 | null | 54087 |
| 12 | quadratic | ranked_reduce | 0.207604 | 1.821 | 0.097459 | 0.095835 | null | 54087 |
| 12 | quadratic | sparse_reduce | 0.216729 | 1.744 | 0.105687 | 0.095770 | null | 54087 |
| 12 | quadratic | stream_reduce | 0.272062 | 1.389 | null | null | 0.256730 | 53621 |
| 12 | quadratic | dense_reduce | 0.177230 | 2.133 | 0.065354 | 0.095063 | null | 54087 |
| 12 | linear_drop | sorted_reduce | 0.420396 | 1.000 | 0.289060 | 0.115917 | null | 57166 |
| 12 | linear_drop | ranked_reduce | 0.230167 | 1.826 | 0.101230 | 0.114854 | null | 57166 |
| 12 | linear_drop | sparse_reduce | 0.248938 | 1.689 | 0.119000 | 0.113168 | null | 57166 |
| 12 | linear_drop | stream_reduce | 0.285562 | 1.472 | null | null | 0.269064 | 56664 |
| 12 | linear_drop | dense_reduce | 0.198354 | 2.119 | 0.069207 | 0.113187 | null | 57166 |
| 12 | restricted_cycle | sorted_reduce | 0.273041 | 1.000 | 0.189897 | 0.070771 | null | 24168 |
| 12 | restricted_cycle | ranked_reduce | 0.156542 | 1.744 | 0.073292 | 0.069269 | null | 24168 |
| 12 | restricted_cycle | sparse_reduce | 0.158625 | 1.721 | 0.075209 | 0.071083 | null | 24168 |
| 12 | restricted_cycle | stream_reduce | 0.192396 | 1.419 | null | null | 0.179125 | 26026 |
| 12 | restricted_cycle | dense_reduce | 0.139916 | 1.951 | 0.055543 | 0.070563 | null | 24168 |
| 20 | quadratic | sorted_reduce | 0.749188 | 1.000 | 0.494999 | 0.230397 | null | 177666 |
| 20 | quadratic | ranked_reduce | 0.463146 | 1.618 | 0.199207 | 0.230500 | null | 177666 |
| 20 | quadratic | sparse_reduce | 0.460938 | 1.625 | 0.213833 | 0.215207 | null | 177666 |
| 20 | quadratic | stream_reduce | 0.737708 | 1.016 | null | null | 0.704812 | 263189 |
| 20 | linear_drop | sorted_reduce | 0.863062 | 1.000 | 0.535021 | 0.306522 | null | 206175 |
| 20 | linear_drop | ranked_reduce | 0.559063 | 1.544 | 0.218645 | 0.307333 | null | 206175 |
| 20 | linear_drop | sparse_reduce | 0.558688 | 1.545 | 0.223291 | 0.297793 | null | 206175 |
| 20 | linear_drop | stream_reduce | 0.834271 | 1.035 | null | null | 0.797438 | 284412 |
| 20 | restricted_cycle | sorted_reduce | 0.277292 | 1.000 | 0.200353 | 0.064458 | null | 20122 |
| 20 | restricted_cycle | ranked_reduce | 0.167854 | 1.652 | 0.086772 | 0.064146 | null | 20122 |
| 20 | restricted_cycle | sparse_reduce | 0.163395 | 1.697 | 0.082978 | 0.063356 | null | 20122 |
| 20 | restricted_cycle | stream_reduce | 0.213875 | 1.297 | null | null | 0.197333 | 56631 |
| 28 | quadratic | sorted_reduce | 1.121354 | 1.000 | 0.698874 | 0.386209 | null | 335368 |
| 28 | quadratic | ranked_reduce | 0.758000 | 1.479 | 0.321208 | 0.363417 | null | 335368 |
| 28 | quadratic | sparse_reduce | 0.721750 | 1.554 | 0.288251 | 0.365187 | null | 335368 |
| 28 | quadratic | stream_reduce | 1.211105 | 0.926 | null | null | 1.140812 | 757311 |
| 28 | linear_drop | sorted_reduce | 1.616271 | 1.000 | 0.838333 | 0.716458 | null | 463156 |
| 28 | linear_drop | ranked_reduce | 1.192209 | 1.356 | 0.399104 | 0.711437 | null | 463156 |
| 28 | linear_drop | sparse_reduce | 1.125333 | 1.436 | 0.331875 | 0.703229 | null | 463156 |
| 28 | linear_drop | stream_reduce | 1.621604 | 0.997 | null | null | 1.532124 | 865478 |
| 28 | restricted_cycle | sorted_reduce | 0.292583 | 1.000 | 0.214688 | 0.066521 | null | 21048 |
| 28 | restricted_cycle | ranked_reduce | 0.211417 | 1.384 | 0.116791 | 0.066792 | null | 21048 |
| 28 | restricted_cycle | sparse_reduce | 0.181188 | 1.615 | 0.087333 | 0.064918 | null | 21048 |
| 28 | restricted_cycle | stream_reduce | 0.270771 | 1.081 | null | null | 0.243896 | 128374 |
| 36 | quadratic | sorted_reduce | 1.574334 | 1.000 | 0.930167 | 0.574896 | null | 513502 |
| 36 | quadratic | ranked_reduce | 1.253812 | 1.256 | 0.536522 | 0.570562 | null | 513502 |
| 36 | quadratic | sparse_reduce | 1.069562 | 1.472 | 0.362918 | 0.561523 | null | 513502 |
| 36 | quadratic | stream_reduce | 2.023604 | 0.778 | null | null | 1.866792 | 1709488 |
| 36 | linear_drop | sorted_reduce | 2.632833 | 1.000 | 1.195729 | 1.320896 | null | 850376 |
| 36 | linear_drop | ranked_reduce | 2.228771 | 1.181 | 0.717291 | 1.325832 | null | 850376 |
| 36 | linear_drop | sparse_reduce | 1.918313 | 1.372 | 0.435125 | 1.309230 | null | 850376 |
| 36 | linear_drop | stream_reduce | 3.052583 | 0.862 | null | null | 2.868521 | 2047441 |
| 36 | restricted_cycle | sorted_reduce | 0.286104 | 1.000 | 0.207937 | 0.065127 | null | 19586 |
| 36 | restricted_cycle | ranked_reduce | 0.268604 | 1.065 | 0.151917 | 0.066084 | null | 19586 |
| 36 | restricted_cycle | sparse_reduce | 0.208437 | 1.373 | 0.093669 | 0.064519 | null | 19586 |
| 36 | restricted_cycle | stream_reduce | 0.347917 | 0.822 | null | null | 0.297000 | 235450 |

Fused construction/reduction phases are not separately observable and remain null, not zero. Source rows and logical row-XOR counts agree across all arms; physical word-XOR counts may differ with ambient coordinate width. These counters cover reduction only, not all calibrated operations.

Retained context bytes and end-of-insertion basis storage are recorded separately. Neither is whole-process peak memory. Worker RSS includes all arms and the common reference corpus. Fixture and independent oracle construction are outside arm timing and inside process receipts.

This completes the bounded construction-plus-RREF task. It does not compute complete Groebner closure, enumerate roots, solve the original polynomial system, or measure index-calculus performance. Those costs and rho ratios remain null.
