# Exact minors verified; the fixed prepass does not justify a complete candidate

The four selected augmented minors are constructed exactly, including Boolean
product cancellations and degree drops. Every minor truth bit matches independent
numeric elimination. The necessary affine masks also match their separate oracle,
and every original solution survives both filters. All **5,280 observations across
12 discovery systems and 55 arms** complete and verify.

This representation does not justify a complete candidate. In all six measured
size/family groups, even granting all subsequent recovery free leaves the optimistic
reference/preprocessing upper confidence bound below **0.082**. The complete
reference already finishes sooner than this measured prepass. No new complete
solver or holdout campaign was launched, and no speedup is claimed.

## Conditional cost boundary

For a full solver retaining this fixed complete-mask prepass,
`T_complete >= T_prepass`. The table divides the fastest retained complete-reference
total by the symbolic arm's preprocessing time, excluding external mask-validation
time from the denominator to favor the candidate. Intervals use the fixed paired
discovery observations. Instrumentation and cap checks remain part of the measured
compiler; the boundary is conditional on keeping this implementation and architecture.
It does not establish a bound for another compiler, representation or early exit.

| Variables | Family | Median reference / prepass | 95% interval | Decision |
|---:|---|---:|---|---|
| 12 | cross_planted | 0.011738 | [0.009875, 0.013883] | BLOCKS_MEASURED_FIXED_PREPROCESSING |
| 12 | planted | 0.010903 | [0.008474, 0.015090] | BLOCKS_MEASURED_FIXED_PREPROCESSING |
| 12 | unplanted | 0.012081 | [0.009786, 0.013399] | BLOCKS_MEASURED_FIXED_PREPROCESSING |
| 16 | cross_planted | 0.019185 | [0.011461, 0.056777] | BLOCKS_MEASURED_FIXED_PREPROCESSING |
| 16 | planted | 0.022413 | [0.007921, 0.056610] | BLOCKS_MEASURED_FIXED_PREPROCESSING |
| 16 | unplanted | 0.047025 | [0.034008, 0.081790] | BLOCKS_MEASURED_FIXED_PREPROCESSING |

## Complete solver and necessary-filter costs

All values are milliseconds per cold arm plus external validation, pooled over
the two n16 discovery seeds and eight random/reverse repetitions. The kind column
separates complete solves from filter-only work. The final ratio is descriptive:
the retained dispatcher's planted median divided by the arm's planted median.
It is not a complete-solver speedup for a filter row. All 52 retained solvers remain
in the same binary alongside the three new filter arms.

| Arm | Kind | Planted (ms) | Cross-planted (ms) | Unplanted (ms) | Dispatcher / arm, planted | Correctness |
|---|---|---:|---:|---:|---:|---|
| search | complete solve | 3.161792 | 3.622229 | 3.848584 | 0.0015 | PASS |
| flat | complete solve | 6.724417 | 3.217792 | 7.329375 | 0.0007 | PASS |
| bucket | complete solve | 63.313041 | 28.751771 | 60.618333 | 0.0001 | PASS |
| hybrid | complete solve | 14.092563 | 7.393916 | 15.796083 | 0.0003 | PASS |
| small_flat | complete solve | 4.851021 | 2.852708 | 5.503917 | 0.0010 | PASS |
| word_tail | complete solve | 4.267812 | 2.452062 | 4.680438 | 0.0011 | PASS |
| merge_search | complete solve | 2.786104 | 3.276916 | 3.520771 | 0.0017 | PASS |
| quadratic_state | complete solve | 1.140750 | 1.227062 | 1.284312 | 0.0042 | PASS |
| packed_state | complete solve | 0.693875 | 0.695875 | 0.867395 | 0.0069 | PASS |
| basis_list | complete solve | 2.276937 | 1.975521 | 9.229458 | 0.0021 | PASS |
| basis_wide | complete solve | 0.459854 | 0.462166 | 1.790209 | 0.0105 | PASS |
| tail_list | complete solve | 5.717791 | 3.070208 | 6.134750 | 0.0008 | PASS |
| tail_wide | complete solve | 1.346812 | 0.748834 | 1.615541 | 0.0036 | PASS |
| packed_untraced | complete solve | 0.481687 | 0.557000 | 0.594771 | 0.0100 | PASS |
| affine_sl_basis_list | complete solve | 10.589729 | 6.691959 | 17.195958 | 0.0005 | PASS |
| affine_sl_basis_fast | complete solve | 0.802208 | 0.546333 | 1.251105 | 0.0060 | PASS |
| gray_scalar | complete solve | 0.031771 | 0.025229 | 0.083104 | 0.1515 | PASS |
| gray_simd | complete solve | 0.005895 | 0.005834 | 0.015709 | 0.8163 | PASS |
| packed_gray12_scalar | complete solve | 0.109917 | 0.114625 | 0.123458 | 0.0438 | PASS |
| packed_gray12_simd | complete solve | 0.050188 | 0.044916 | 0.054937 | 0.0959 | PASS |
| packed_gray16_scalar | complete solve | 0.034959 | 0.031355 | 0.086146 | 0.1377 | PASS |
| packed_gray16_simd | complete solve | 0.010333 | 0.009458 | 0.019291 | 0.4657 | PASS |
| gray_delta_scalar | complete solve | 0.025916 | 0.031626 | 0.083896 | 0.1857 | PASS |
| gray_delta_simd | complete solve | 0.007104 | 0.005959 | 0.013999 | 0.6774 | PASS |
| initial_list | complete solve | 0.027354 | 0.045437 | 0.085104 | 0.1759 | PASS |
| initial_simd | complete solve | 0.007187 | 0.013750 | 0.015354 | 0.6696 | PASS |
| packed_gray16_delta_scalar | complete solve | 0.034354 | 0.036333 | 0.088041 | 0.1401 | PASS |
| packed_gray16_delta_simd | complete solve | 0.009876 | 0.009874 | 0.017229 | 0.4873 | PASS |
| fiber_rows | complete solve | 0.034792 | 0.022396 | 0.075521 | 0.1383 | PASS |
| fiber_columns | complete solve | 0.032771 | 0.023334 | 0.070708 | 0.1469 | PASS |
| fiber_simd | complete solve | 0.014250 | 0.011938 | 0.037583 | 0.3377 | PASS |
| fiber_zero_simd | complete solve | 0.014209 | 0.013833 | 0.055083 | 0.3387 | PASS |
| gray_quiet | complete solve | 0.006021 | 0.006167 | 0.016063 | 0.7993 | PASS |
| wide64_quiet | complete solve | 0.005916 | 0.004938 | 0.010021 | 0.8135 | PASS |
| byte_scalar | complete solve | 0.015459 | 0.013104 | 0.034042 | 0.3113 | PASS |
| byte_simd | complete solve | 0.006042 | 0.005521 | 0.013416 | 0.7965 | PASS |
| leaf16_quiet | complete solve | 0.009000 | 0.010042 | 0.019166 | 0.5347 | PASS |
| leaf16_byte_scalar | complete solve | 0.018021 | 0.015583 | 0.036166 | 0.2670 | PASS |
| leaf16_byte_simd | complete solve | 0.009021 | 0.009437 | 0.017167 | 0.5335 | PASS |
| byte_single_quiet | complete solve | 0.005396 | 0.006208 | 0.014020 | 0.8919 | PASS |
| leaf16_single_quiet | complete solve | 0.009666 | 0.009187 | 0.017542 | 0.4979 | PASS |
| byte_planes | complete solve | 0.007062 | 0.006625 | 0.012729 | 0.6815 | PASS |
| leaf16_byte_planes | complete solve | 0.009104 | 0.009270 | 0.016459 | 0.5286 | PASS |
| byte_unrolled | complete solve | 0.009166 | 0.007625 | 0.012542 | 0.5250 | PASS |
| wide64_unrolled | complete solve | 0.007374 | 0.007333 | 0.009521 | 0.6526 | PASS |
| leaf16_byte_unrolled | complete solve | 0.011730 | 0.011833 | 0.016958 | 0.4103 | PASS |
| word16_unrolled | complete solve | 0.004063 | 0.003916 | 0.008292 | 1.1846 | PASS |
| word_dispatch | complete solve | 0.004812 | 0.004916 | 0.008833 | 1.0000 | PASS |
| leaf16_word_unrolled | complete solve | 0.008688 | 0.007812 | 0.011500 | 0.5540 | PASS |
| projected4 | complete solve | 0.015041 | 0.013771 | 0.032250 | 0.3199 | PASS |
| projected5 | complete solve | 0.014958 | 0.010500 | 0.033562 | 0.3217 | PASS |
| projected6 | complete solve | 0.022708 | 0.020688 | 0.112792 | 0.2119 | PASS |
| numeric_minor4 | necessary filter | 0.605666 | 0.608313 | 0.674958 | 0.0079 | PASS |
| symbolic_minor4 | necessary filter | 0.210646 | 0.183667 | 0.149229 | 0.0228 | PASS |
| affine_filter | necessary filter | 0.043229 | 0.041521 | 0.042271 | 0.1113 | PASS |

## Exact structural evidence

The fixed minors leave **1,410–1,920** of 4,096 outside assignments at n16, while
exact affine consistency leaves **5–15**. Both predicates remain necessary only;
the final column gives independently enumerated original solution counts. Four
minors drop to degree five in the cross-planted cases; the remaining 44 have
degree six. No chosen minor is zero or duplicated on this discovery grid. Separate
tests retain those degeneracies and the rank-deficient insufficiency counterexample.

| Variables | Seed | Family | Outside assignments | Minor survivors | Affine survivors | Original solutions | Degrees | ANF supports |
|---:|---:|---|---:|---:|---:|---:|---|---|
| 12 | 17 | planted | 256 | 103 | 15 | 1 | [6, 6, 6, 6] | [116, 121, 130, 124] |
| 12 | 17 | cross_planted | 256 | 119 | 25 | 1 | [5, 6, 6, 6] | [76, 130, 130, 124] |
| 12 | 17 | unplanted | 256 | 104 | 25 | 2 | [6, 6, 6, 6] | [122, 114, 120, 136] |
| 12 | 937 | planted | 256 | 101 | 14 | 1 | [6, 6, 6, 6] | [108, 90, 114, 106] |
| 12 | 937 | cross_planted | 256 | 106 | 14 | 1 | [5, 6, 6, 6] | [88, 106, 114, 106] |
| 12 | 937 | unplanted | 256 | 99 | 7 | 1 | [6, 6, 6, 6] | [114, 122, 132, 134] |
| 16 | 17 | planted | 4096 | 1855 | 5 | 1 | [6, 6, 6, 6] | [314, 313, 427, 388] |
| 16 | 17 | cross_planted | 4096 | 1920 | 7 | 1 | [5, 6, 6, 6] | [148, 501, 427, 388] |
| 16 | 17 | unplanted | 4096 | 1410 | 13 | 0 | [6, 6, 6, 6] | [696, 700, 574, 536] |
| 16 | 937 | planted | 4096 | 1667 | 15 | 1 | [6, 6, 6, 6] | [1079, 1025, 1184, 1238] |
| 16 | 937 | cross_planted | 4096 | 1762 | 11 | 1 | [5, 6, 6, 6] | [108, 1096, 1138, 1200] |
| 16 | 937 | unplanted | 4096 | 1496 | 8 | 0 | [6, 6, 6, 6] | [1128, 1124, 978, 918] |

Individual rejection masks, their exact pairwise intersections and their joint
union are retained in the result. Their rates are not added or assumed independent.
The verifier checks **417,792 original assignments** across the complete discovery
grid. A vanishing determinant cannot substitute for an original solution.

The campaign took **6.472 seconds**. Whole-worker peak RSS was
**6,881,280 bytes**, covering all arms and oracle preparation.
The symbolic DP coefficient payload peaks at 1,024 bytes at n12 and 16,384 bytes
at n16; these are that component's word arrays, not total candidate memory or
allocator peaks. Other candidate-specific memory remains unmeasured.

Validation comprises 76 Rust tests and 15 Python evidence checks, including
independent truth-table replay, explicit permutation parity, cancellation and
degree-drop cases, false models, changed sources/receipts, overlap accounting and
censored-construction rejection. These are producer checks, not external review.

`README.md` states the exact scope and timing exclusions. The frozen protocol,
source, raw samples, receipts and result are sealed in `run_01/manifest.json`;
`RUN_LEDGER.json` binds this report. The earlier projected-fiber result remains
unchanged. `NEXT_EXPERIMENT.md` records a distinct, unimplemented 16-bit syndrome
representation hypothesis, with mandatory original-equation checks and fresh
holdouts for any promising complete candidate.

Full-IC, production, calibrated-operation and rho costs remain **null**. This
construction diagnostic does not meet the active dramatic-gain objective.
