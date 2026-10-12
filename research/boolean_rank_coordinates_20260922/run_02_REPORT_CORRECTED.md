# Combinatorial monomial coordinates: scaling experiment

Correctness **PASS**: 144 cells, 12240 batch-arm samples, 167280 oracle-verified outputs.

Dramatic scalable-construction gate: **REJECTED**, 22 / 27 comparisons passed.

Cold batch32 milliseconds below include setup, applications, allocations, exact equality validation and destruction. Values are medians over two holdout seeds and 20 balanced repetitions. Ratios here compare pooled medians; acceptance uses the per-size paired intervals in results.json.

| Variables | Family | Variant | Cold batch (ms) | Sorted / arm | Retained bytes | Lookup entries | Correctness |
|---:|---|---|---:|---:|---:|---:|---|
| 12 | quadratic | sorted | 1.141812 | 1.000 | 0 | 0 | PASS |
| 12 | quadratic | binary | 0.495125 | 2.306 | 9544 | 0 | PASS |
| 12 | quadratic | ranked | 0.457354 | 2.497 | 9960 | 52 | PASS |
| 12 | quadratic | sparse_rank | 0.560312 | 2.038 | 9960 | 52 | PASS |
| 12 | quadratic | dense | 0.314021 | 3.636 | 25928 | 4096 | PASS |
| 12 | linear_drop | sorted | 1.224479 | 1.000 | 0 | 0 | PASS |
| 12 | linear_drop | binary | 0.515875 | 2.374 | 9544 | 0 | PASS |
| 12 | linear_drop | ranked | 0.480687 | 2.547 | 9960 | 52 | PASS |
| 12 | linear_drop | sparse_rank | 0.590521 | 2.074 | 9960 | 52 | PASS |
| 12 | linear_drop | dense | 0.333354 | 3.673 | 25928 | 4096 | PASS |
| 12 | restricted_cycle | sorted | 0.847771 | 1.000 | 0 | 0 | PASS |
| 12 | restricted_cycle | binary | 0.396208 | 2.140 | 5960 | 0 | PASS |
| 12 | restricted_cycle | ranked | 0.365001 | 2.323 | 6376 | 52 | PASS |
| 12 | restricted_cycle | sparse_rank | 0.431667 | 1.964 | 6376 | 52 | PASS |
| 12 | restricted_cycle | dense | 0.268375 | 3.159 | 22344 | 4096 | PASS |
| 20 | quadratic | sorted | 1.983083 | 1.000 | 0 | 0 | PASS |
| 20 | quadratic | binary | 1.060437 | 1.870 | 35272 | 0 | PASS |
| 20 | quadratic | ranked | 0.858896 | 2.309 | 35944 | 84 | PASS |
| 20 | quadratic | sparse_rank | 0.964625 | 2.056 | 35944 | 84 | PASS |
| 20 | linear_drop | sorted | 2.267291 | 1.000 | 0 | 0 | PASS |
| 20 | linear_drop | binary | 1.218166 | 1.861 | 35272 | 0 | PASS |
| 20 | linear_drop | ranked | 0.987604 | 2.296 | 35944 | 84 | PASS |
| 20 | linear_drop | sparse_rank | 1.075104 | 2.109 | 35944 | 84 | PASS |
| 20 | restricted_cycle | sorted | 0.851583 | 1.000 | 0 | 0 | PASS |
| 20 | restricted_cycle | binary | 0.496896 | 1.714 | 18248 | 0 | PASS |
| 20 | restricted_cycle | ranked | 0.400562 | 2.126 | 18920 | 84 | PASS |
| 20 | restricted_cycle | sparse_rank | 0.416042 | 2.047 | 18920 | 84 | PASS |
| 28 | quadratic | sorted | 2.719313 | 1.000 | 0 | 0 | PASS |
| 28 | quadratic | binary | 1.762417 | 1.543 | 70088 | 0 | PASS |
| 28 | quadratic | ranked | 1.367083 | 1.989 | 71016 | 116 | PASS |
| 28 | quadratic | sparse_rank | 1.266188 | 2.148 | 71016 | 116 | PASS |
| 28 | linear_drop | sorted | 3.301563 | 1.000 | 0 | 0 | PASS |
| 28 | linear_drop | binary | 2.150104 | 1.536 | 70088 | 0 | PASS |
| 28 | linear_drop | ranked | 1.691250 | 1.952 | 71016 | 116 | PASS |
| 28 | linear_drop | sparse_rank | 1.469375 | 2.247 | 71016 | 116 | PASS |
| 28 | restricted_cycle | sorted | 0.835229 | 1.000 | 0 | 0 | PASS |
| 28 | restricted_cycle | binary | 0.625604 | 1.335 | 34632 | 0 | PASS |
| 28 | restricted_cycle | ranked | 0.483375 | 1.728 | 35560 | 116 | PASS |
| 28 | restricted_cycle | sparse_rank | 0.414729 | 2.014 | 35560 | 116 | PASS |
| 36 | quadratic | sorted | 3.581979 | 1.000 | 0 | 0 | PASS |
| 36 | quadratic | binary | 2.868584 | 1.249 | 139976 | 0 | PASS |
| 36 | quadratic | ranked | 2.232979 | 1.604 | 141160 | 148 | PASS |
| 36 | quadratic | sparse_rank | 1.666313 | 2.150 | 141160 | 148 | PASS |
| 36 | linear_drop | sorted | 4.658271 | 1.000 | 0 | 0 | PASS |
| 36 | linear_drop | binary | 3.937604 | 1.183 | 139976 | 0 | PASS |
| 36 | linear_drop | ranked | 3.075438 | 1.515 | 141160 | 148 | PASS |
| 36 | linear_drop | sparse_rank | 2.124417 | 2.193 | 141160 | 148 | PASS |
| 36 | restricted_cycle | sorted | 0.849187 | 1.000 | 0 | 0 | PASS |
| 36 | restricted_cycle | binary | 0.856833 | 0.991 | 67400 | 0 | PASS |
| 36 | restricted_cycle | ranked | 0.688917 | 1.233 | 68584 | 148 | PASS |
| 36 | restricted_cycle | sparse_rank | 0.469021 | 1.811 | 68584 | 148 | PASS |

The dense lookup executes only at n=12. Larger dense cells are NOT_EXECUTED by design and retain null costs; no failure or hypothetical runtime is inferred.

Ranked coordinates use O(nD) prefix-count entries, but the ambient basis, multiplier lists and dense row widths remain charged. Whole-worker RSS includes common references and all arms. Fixtures and independent oracle construction are outside arm timing and inside process receipts.

These are generic Boolean matrix-construction diagnostics. No curve inputs, scalar recovery, relation collection, production solver integration, full-solver timing or rho comparison occurs.

This additive run evaluates sparse_rank against the fresh sorted, binary and ranked controls. Only touched nonzero intermediate words are retained; the final dense matrix and all validation costs remain charged. Its fresh protocol and holdouts do not replace the initial failed scaling run.
