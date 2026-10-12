# Packed Boolean matrix construction experiment

Same finite Boolean multiplication contract as the parameterized envelope study. A packed coordinate map is shared by both new constructors. The actual output support may shrink with coefficients and cancellations; an occupancy bitmap determines the only columns returned. Precomputed degree-eligible views preserve numeric multiplier order. Compare envelope-specific reuse against packed direct multiplication, not only the older sort-based construction.

Correctness: **PASS** across 256 fixed cells, 25088 batch-arm samples, and 533120 verified matrix outputs.

Changed-input envelope hits: **49728**; support-escape fallbacks: **4704**. Exact-support schedule changed hits: **0**.

Cold performance promotion: **REJECTED**, 45 / 48 gates passed. Newly eligible multipliers before cancellation, over the fixture grid without arm/repetition multiplication: **60816**.

Dramatic 2x construction gate: **{'packed_envelope': 'REJECTED', 'packed_direct': 'PASS'}**. Envelope-specific gain over packed direct construction: **REJECTED**.

The table reports cold batch construction plus output validation in milliseconds, at batch 64 on holdouts. Medians pool four sizes, two seeds and fourteen repetitions; paired intervals remain separate by size and family in results.json. Ratios below are descriptive ratios of pooled medians against direct construction and packed-matrix caching. These are construction-stage diagnostics.

| Family | Variant | Cold batch (ms) | Direct / arm | Matrix cache / arm | Retained bytes | Hits / 64 | Correctness |
|---|---|---:|---:|---:|---:|---:|---|
| repeat | direct | 1.527208 | 1.000 | 0.108 | 0 | 0 | PASS |
| repeat | layout | 1.365522 | 1.118 | 0.121 | 3352 | 64 | PASS |
| repeat | schedule | 0.196729 | 7.763 | 0.839 | 5221 | 64 | PASS |
| repeat | matrix_cache | 0.165104 | 9.250 | 1.000 | 4875 | 64 | PASS |
| repeat | envelope | 1.538896 | 0.992 | 0.107 | 100929 | 64 | PASS |
| repeat | packed_envelope | 0.697812 | 2.189 | 0.237 | 96881 | 64 | PASS |
| repeat | packed_direct | 0.423021 | 3.610 | 0.390 | 2386 | 64 | PASS |
| coefficients | direct | 1.448730 | 1.000 | 1.003 | 0 | 0 | PASS |
| coefficients | layout | 1.540104 | 0.941 | 0.943 | 3352 | 1 | PASS |
| coefficients | schedule | 1.438396 | 1.007 | 1.010 | 5221 | 1 | PASS |
| coefficients | matrix_cache | 1.453021 | 0.997 | 1.000 | 4875 | 1 | PASS |
| coefficients | envelope | 1.836625 | 0.789 | 0.791 | 100929 | 64 | PASS |
| coefficients | packed_envelope | 0.586687 | 2.469 | 2.477 | 96881 | 64 | PASS |
| coefficients | packed_direct | 0.319625 | 4.533 | 4.546 | 2386 | 64 | PASS |
| degree_cycle | direct | 1.609520 | 1.000 | 0.995 | 0 | 0 | PASS |
| degree_cycle | layout | 1.661208 | 0.969 | 0.964 | 3352 | 1 | PASS |
| degree_cycle | schedule | 1.623917 | 0.991 | 0.986 | 5221 | 1 | PASS |
| degree_cycle | matrix_cache | 1.601625 | 1.005 | 1.000 | 4875 | 1 | PASS |
| degree_cycle | envelope | 1.982729 | 0.812 | 0.808 | 100929 | 64 | PASS |
| degree_cycle | packed_envelope | 0.647167 | 2.487 | 2.475 | 96881 | 64 | PASS |
| degree_cycle | packed_direct | 0.371125 | 4.337 | 4.316 | 2386 | 64 | PASS |
| escape | direct | 1.453626 | 1.000 | 0.967 | 0 | 0 | PASS |
| escape | layout | 1.472729 | 0.987 | 0.955 | 3352 | 1 | PASS |
| escape | schedule | 1.436917 | 1.012 | 0.979 | 5221 | 1 | PASS |
| escape | matrix_cache | 1.406083 | 1.034 | 1.000 | 4875 | 1 | PASS |
| escape | envelope | 1.768812 | 0.822 | 0.795 | 100929 | 48 | PASS |
| escape | packed_envelope | 0.823937 | 1.764 | 1.707 | 96881 | 48 | PASS |
| escape | packed_direct | 0.313521 | 4.636 | 4.485 | 2386 | 64 | PASS |

Cold totals charge compilation, guards, applications, fallbacks, fresh output, validation and destruction. Common fixture/envelope generation and independent reference generation are outside timing; process receipts include them. Retained bytes exclude allocator metadata. Worker RSS covers all variants plus the reference corpus, not candidate-specific peak usage.

Both new arms use packed ambient coordinates and occupancy-based output compaction. The packed-direct arm performs current polynomial products by XOR and does not retain coefficient routes or generator envelopes. The retained envelope explicitly represents varying coefficients. It is not a cached completed matrix. Product parity, actual generator degree, newly eligible multipliers and actual output columns are evaluated at each application. Support escapes leave the retained plan unchanged and rebuild directly.

No production solver path changes. This is finite generic Boolean matrix construction, with no curve inputs, solving, relation collection, target import, scalar recovery or rho comparison. Full-pipeline costs and normalized cryptanalytic ratios remain null.
