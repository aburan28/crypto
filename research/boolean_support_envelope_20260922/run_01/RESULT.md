# Parameterized Boolean support-envelope experiment

A fixed support envelope E_j represents p_j(c)=sum_{m in E_j} c_{j,m} x^m over F2. Coefficients may change. Compile a parity mask per product monomial and multiplier. At application recompute actual generator degrees, select all degree-eligible multipliers, evaluate parity, remove zero rows and compact actual output columns. Exact-support caching cannot reuse changed inputs; this different contract can.

Correctness: **PASS** across 256 fixed cells, 12800 batch-arm samples, and 272000 verified matrix outputs.

Changed-input envelope hits: **35520**; support-escape fallbacks: **3360**. Exact-support schedule changed hits: **0**.

Cold performance promotion: **REJECTED**, 0 / 48 gates passed. Newly eligible multipliers before cancellation, over the fixture grid without arm/repetition multiplication: **60816**.

The table reports cold batch construction plus output validation in milliseconds, at batch 64 on holdouts. Medians pool four sizes, two seeds and ten repetitions; paired intervals remain separate by size and family in results.json. Ratios below are descriptive ratios of pooled medians against direct construction and packed-matrix caching. These are construction-stage diagnostics.

| Family | Variant | Cold batch (ms) | Direct / arm | Matrix cache / arm | Retained bytes | Hits / 64 | Correctness |
|---|---|---:|---:|---:|---:|---:|---|
| repeat | direct | 1.657479 | 1.000 | 0.100 | 0 | 0 | PASS |
| repeat | layout | 1.399625 | 1.184 | 0.118 | 3352 | 64 | PASS |
| repeat | schedule | 0.215458 | 7.693 | 0.770 | 5234 | 64 | PASS |
| repeat | matrix_cache | 0.165813 | 9.996 | 1.000 | 4878 | 64 | PASS |
| repeat | envelope | 1.610500 | 1.029 | 0.103 | 100715 | 64 | PASS |
| coefficients | direct | 1.594291 | 1.000 | 0.993 | 0 | 0 | PASS |
| coefficients | layout | 1.659000 | 0.961 | 0.954 | 3352 | 1 | PASS |
| coefficients | schedule | 1.631145 | 0.977 | 0.970 | 5234 | 1 | PASS |
| coefficients | matrix_cache | 1.582396 | 1.008 | 1.000 | 4878 | 1 | PASS |
| coefficients | envelope | 2.081271 | 0.766 | 0.760 | 100715 | 64 | PASS |
| degree_cycle | direct | 1.825187 | 1.000 | 0.953 | 0 | 0 | PASS |
| degree_cycle | layout | 2.159729 | 0.845 | 0.805 | 3352 | 1 | PASS |
| degree_cycle | schedule | 1.816083 | 1.005 | 0.958 | 5234 | 1 | PASS |
| degree_cycle | matrix_cache | 1.739000 | 1.050 | 1.000 | 4878 | 1 | PASS |
| degree_cycle | envelope | 2.184021 | 0.836 | 0.796 | 100715 | 64 | PASS |
| escape | direct | 1.567375 | 1.000 | 1.031 | 0 | 0 | PASS |
| escape | layout | 1.562709 | 1.003 | 1.034 | 3352 | 1 | PASS |
| escape | schedule | 1.550729 | 1.011 | 1.042 | 5234 | 1 | PASS |
| escape | matrix_cache | 1.615521 | 0.970 | 1.000 | 4878 | 1 | PASS |
| escape | envelope | 2.195896 | 0.714 | 0.736 | 100715 | 48 | PASS |

Cold totals charge compilation, guards, applications, fallbacks, fresh output, validation and destruction. Common fixture/envelope generation and independent reference generation are outside timing; process receipts include them. Retained bytes exclude allocator metadata. Worker RSS covers all variants plus the reference corpus, not candidate-specific peak usage.

The retained envelope explicitly represents varying coefficients. It is not a cached completed matrix. Product parity, actual generator degree, newly eligible multipliers and actual output columns are evaluated at each application. Support escapes leave the retained plan unchanged and rebuild directly.

No production solver path changes. This is finite generic Boolean matrix construction, with no curve inputs, solving, relation collection, target import, scalar recovery or rho comparison. Full-pipeline costs and normalized cryptanalytic ratios remain null.
