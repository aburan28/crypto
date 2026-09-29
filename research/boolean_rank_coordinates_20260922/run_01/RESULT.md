# Combinatorial monomial coordinates: scaling experiment

Correctness **PASS**: 144 cells, 5616 batch-arm samples, 76752 oracle-verified outputs.

Dramatic scalable-construction gate: **REJECTED**, 11 / 18 comparisons passed.

Cold batch32 milliseconds below include setup, applications, allocations, exact equality validation and destruction. Values are medians over two holdout seeds and twelve balanced repetitions. Ratios here compare pooled medians; acceptance uses the per-size paired intervals in results.json.

| Variables | Family | Variant | Cold batch (ms) | Sorted / arm | Retained bytes | Lookup entries | Correctness |
|---:|---|---|---:|---:|---:|---:|---|
| 12 | quadratic | sorted | 1.136562 | 1.000 | 0 | 0 | PASS |
| 12 | quadratic | binary | 0.499063 | 2.277 | 9544 | 0 | PASS |
| 12 | quadratic | ranked | 0.447813 | 2.538 | 9960 | 52 | PASS |
| 12 | quadratic | dense | 0.312167 | 3.641 | 25928 | 4096 | PASS |
| 12 | linear_drop | sorted | 1.178416 | 1.000 | 0 | 0 | PASS |
| 12 | linear_drop | binary | 0.499937 | 2.357 | 9544 | 0 | PASS |
| 12 | linear_drop | ranked | 0.460417 | 2.559 | 9960 | 52 | PASS |
| 12 | linear_drop | dense | 0.316083 | 3.728 | 25928 | 4096 | PASS |
| 12 | restricted_cycle | sorted | 0.840437 | 1.000 | 0 | 0 | PASS |
| 12 | restricted_cycle | binary | 0.390834 | 2.150 | 5960 | 0 | PASS |
| 12 | restricted_cycle | ranked | 0.359625 | 2.337 | 6376 | 52 | PASS |
| 12 | restricted_cycle | dense | 0.261145 | 3.218 | 22344 | 4096 | PASS |
| 20 | quadratic | sorted | 1.983229 | 1.000 | 0 | 0 | PASS |
| 20 | quadratic | binary | 1.054604 | 1.881 | 35272 | 0 | PASS |
| 20 | quadratic | ranked | 0.872854 | 2.272 | 35944 | 84 | PASS |
| 20 | linear_drop | sorted | 2.222666 | 1.000 | 0 | 0 | PASS |
| 20 | linear_drop | binary | 1.196584 | 1.858 | 35272 | 0 | PASS |
| 20 | linear_drop | ranked | 0.988521 | 2.248 | 35944 | 84 | PASS |
| 20 | restricted_cycle | sorted | 0.844584 | 1.000 | 0 | 0 | PASS |
| 20 | restricted_cycle | binary | 0.512208 | 1.649 | 18248 | 0 | PASS |
| 20 | restricted_cycle | ranked | 0.425542 | 1.985 | 18920 | 84 | PASS |
| 28 | quadratic | sorted | 2.765855 | 1.000 | 0 | 0 | PASS |
| 28 | quadratic | binary | 1.747500 | 1.583 | 70088 | 0 | PASS |
| 28 | quadratic | ranked | 1.406230 | 1.967 | 71016 | 116 | PASS |
| 28 | linear_drop | sorted | 3.323583 | 1.000 | 0 | 0 | PASS |
| 28 | linear_drop | binary | 2.204188 | 1.508 | 70088 | 0 | PASS |
| 28 | linear_drop | ranked | 1.708625 | 1.945 | 71016 | 116 | PASS |
| 28 | restricted_cycle | sorted | 0.864958 | 1.000 | 0 | 0 | PASS |
| 28 | restricted_cycle | binary | 0.629938 | 1.373 | 34632 | 0 | PASS |
| 28 | restricted_cycle | ranked | 0.491854 | 1.759 | 35560 | 116 | PASS |
| 36 | quadratic | sorted | 3.619687 | 1.000 | 0 | 0 | PASS |
| 36 | quadratic | binary | 2.831313 | 1.278 | 139976 | 0 | PASS |
| 36 | quadratic | ranked | 2.255750 | 1.605 | 141160 | 148 | PASS |
| 36 | linear_drop | sorted | 4.681813 | 1.000 | 0 | 0 | PASS |
| 36 | linear_drop | binary | 3.987875 | 1.174 | 139976 | 0 | PASS |
| 36 | linear_drop | ranked | 3.132104 | 1.495 | 141160 | 148 | PASS |
| 36 | restricted_cycle | sorted | 0.856437 | 1.000 | 0 | 0 | PASS |
| 36 | restricted_cycle | binary | 0.867583 | 0.987 | 67400 | 0 | PASS |
| 36 | restricted_cycle | ranked | 0.702729 | 1.219 | 68584 | 148 | PASS |

The dense lookup executes only at n=12. Larger dense cells are NOT_EXECUTED by design and retain null costs; no failure or hypothetical runtime is inferred.

Ranked coordinates use O(nD) prefix-count entries, but the ambient basis, multiplier lists and dense row widths remain charged. Whole-worker RSS includes common references and all arms. Fixtures and independent oracle construction are outside arm timing and inside process receipts.

These are generic Boolean matrix-construction diagnostics. No curve inputs, scalar recovery, relation collection, production solver integration, full-solver timing or rho comparison occurs.
