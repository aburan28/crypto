# Boolean symbolic-product schedule experiment

Over F2, complete normalized support specifies every coefficient. An exact full-system support hit therefore repeats an identical system and matrix, not a distinct coefficient instance.

Correctness: **PASS** across 256 fixed cells, 8192 batch-arm samples, and 174080 verified matrix outputs.

Cross-instance schedule hits: **0**. Performance promotion against both retained controls: **REJECTED**.

The table below reports cold batch construction plus output validation in milliseconds, at batch 64 on holdout fixtures. Values are medians over four variable sizes, two seeds and eight repetitions; the per-size acceptance intervals remain in results.json. These are standalone construction diagnostics, not solver or cryptanalytic performance.

| Family | Variant | Cold batch median (ms) | Retained structure median (bytes) | Median hits / 64 |
|---|---|---:|---:|---:|
| repeat | direct | 1.526501 | 0 | 0 |
| repeat | layout | 1.340375 | 3352 | 64 |
| repeat | schedule | 0.195416 | 8302 | 64 |
| repeat | matrix_cache | 0.171000 | 4875 | 64 |
| constant_toggle | direct | 1.555333 | 0 | 0 |
| constant_toggle | layout | 1.457251 | 3352 | 32 |
| constant_toggle | schedule | 0.897500 | 8302 | 32 |
| constant_toggle | matrix_cache | 0.930312 | 4875 | 32 |
| term_toggle | direct | 1.602645 | 0 | 0 |
| term_toggle | layout | 1.637500 | 3352 | 64 |
| term_toggle | schedule | 1.077271 | 8302 | 32 |
| term_toggle | matrix_cache | 1.056000 | 4875 | 32 |
| degree_drop | direct | 1.578729 | 0 | 0 |
| degree_drop | layout | 1.450813 | 3352 | 64 |
| degree_drop | schedule | 0.910521 | 8302 | 32 |
| degree_drop | matrix_cache | 0.898791 | 4875 | 32 |

Cache setup, exact-key checks, fallbacks, output allocation, output validation and destruction are charged. Common fixture and independent reference generation are outside timing; fresh-process receipts cover the whole worker, including those costs. Retained bytes exclude allocator metadata; RSS includes all variants and the common reference corpus.

Changing a coefficient from 1 to 0 over F2 changes support. Reusing a schedule across such changes would require a different, parameterized support-envelope contract with explicit cancellation and degree-drop rules. This experiment does not implement or validate that different contract.

No production solver path is changed. The full-support cache key is evaluated exactly as proposed, including the stronger packed-matrix-cache control. No curve targets, key work, relation collection, or rho comparison are part of this study.
