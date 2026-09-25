# Packed Boolean matrix construction experiment

Same finite Boolean multiplication contract as the parameterized envelope study. A packed coordinate map is shared by both new constructors. The actual output support may shrink with coefficients and cancellations; an occupancy bitmap determines the only columns returned. Precomputed degree-eligible views preserve numeric multiplier order. Compare envelope-specific reuse against packed direct multiplication, not only the older sort-based construction.

Correctness: **PASS** across 256 fixed cells, 25088 batch-arm samples, and 533120 verified matrix outputs.

Changed-input envelope hits: **49728**; support-escape fallbacks: **4704**. Exact-support schedule changed hits: **0**.

Cold performance promotion: **REJECTED**, 45 / 48 gates passed. Newly eligible multipliers before cancellation, over the fixture grid without arm/repetition multiplication: **60816**.

Dramatic 2x construction gate: **{'packed_envelope': 'REJECTED', 'packed_direct': 'PASS'}**. Envelope-specific gain over packed direct construction: **REJECTED**.

The table reports cold batch construction plus output validation in milliseconds, at batch 64 on holdouts. Medians pool four sizes, two seeds and fourteen repetitions; paired intervals remain separate by size and family in results.json. Ratios below are descriptive ratios of pooled medians against direct construction and packed-matrix caching. These are construction-stage diagnostics.

| Family | Variant | Cold batch (ms) | Direct / arm | Matrix cache / arm | Retained bytes | Hits / 64 | Correctness |
|---|---|---:|---:|---:|---:|---:|---|
| repeat | direct | 1.513688 | 1.000 | 0.110 | 0 | 0 | PASS |
| repeat | layout | 1.355354 | 1.117 | 0.123 | 3352 | 64 | PASS |
| repeat | schedule | 0.194688 | 7.775 | 0.858 | 5256 | 64 | PASS |
| repeat | matrix_cache | 0.167021 | 9.063 | 1.000 | 4878 | 64 | PASS |
| repeat | envelope | 1.539625 | 0.983 | 0.108 | 100939 | 64 | PASS |
| repeat | packed_envelope | 0.663104 | 2.283 | 0.252 | 96891 | 64 | PASS |
| repeat | packed_direct | 0.401292 | 3.772 | 0.416 | 2386 | 64 | PASS |
| coefficients | direct | 1.571417 | 1.000 | 0.923 | 0 | 0 | PASS |
| coefficients | layout | 1.809438 | 0.868 | 0.802 | 3352 | 1 | PASS |
| coefficients | schedule | 1.759542 | 0.893 | 0.825 | 5256 | 1 | PASS |
| coefficients | matrix_cache | 1.450812 | 1.083 | 1.000 | 4878 | 1 | PASS |
| coefficients | envelope | 1.847271 | 0.851 | 0.785 | 100939 | 64 | PASS |
| coefficients | packed_envelope | 0.775875 | 2.025 | 1.870 | 96891 | 64 | PASS |
| coefficients | packed_direct | 0.332417 | 4.727 | 4.364 | 2386 | 64 | PASS |
| degree_cycle | direct | 1.639750 | 1.000 | 0.984 | 0 | 0 | PASS |
| degree_cycle | layout | 1.637417 | 1.001 | 0.985 | 3352 | 1 | PASS |
| degree_cycle | schedule | 1.614562 | 1.016 | 0.999 | 5256 | 1 | PASS |
| degree_cycle | matrix_cache | 1.613021 | 1.017 | 1.000 | 4878 | 1 | PASS |
| degree_cycle | envelope | 1.989396 | 0.824 | 0.811 | 100939 | 64 | PASS |
| degree_cycle | packed_envelope | 0.613354 | 2.673 | 2.630 | 96891 | 64 | PASS |
| degree_cycle | packed_direct | 0.343000 | 4.781 | 4.703 | 2386 | 64 | PASS |
| escape | direct | 1.445709 | 1.000 | 0.979 | 0 | 0 | PASS |
| escape | layout | 1.787583 | 0.809 | 0.792 | 3352 | 1 | PASS |
| escape | schedule | 1.480979 | 0.976 | 0.956 | 5256 | 1 | PASS |
| escape | matrix_cache | 1.415104 | 1.022 | 1.000 | 4878 | 1 | PASS |
| escape | envelope | 1.785208 | 0.810 | 0.793 | 100939 | 48 | PASS |
| escape | packed_envelope | 0.839458 | 1.722 | 1.686 | 96891 | 48 | PASS |
| escape | packed_direct | 0.338438 | 4.272 | 4.181 | 2386 | 64 | PASS |

Cold totals charge compilation, guards, applications, fallbacks, fresh output, validation and destruction. Common fixture/envelope generation and independent reference generation are outside timing; process receipts include them. Retained bytes exclude allocator metadata. Worker RSS covers all variants plus the reference corpus, not candidate-specific peak usage.

Both new arms use packed ambient coordinates and occupancy-based output compaction. The packed-direct arm performs current polynomial products by XOR and does not retain coefficient routes or generator envelopes. The retained envelope explicitly represents varying coefficients. It is not a cached completed matrix. Product parity, actual generator degree, newly eligible multipliers and actual output columns are evaluated at each application. Support escapes leave the retained plan unchanged and rebuild directly.

No production solver path changes. This is finite generic Boolean matrix construction, with no curve inputs, solving, relation collection, target import, scalar recovery or rho comparison. Full-pipeline costs and normalized cryptanalytic ratios remain null.
