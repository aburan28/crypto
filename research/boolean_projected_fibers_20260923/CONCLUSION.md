# Complete recovery is verified; the dramatic-gain gate is rejected

This implements the exact projected-fiber contract on bounded generated Boolean
quadratic systems. All **89,856 observations on 216 systems and
52 methods** completed with verified outcomes. Each candidate computes a necessary
affine system, recovers its complete free-variable space, and checks the original
equations before accepting a solution. Degree changes and coefficient cancellations
are exact XOR operations; a rank drop never discards required extensions.

| Candidate | Groups above 2.0 | Groups above 1.0 | Dramatic decision |
|---|---:|---:|---|
| projected4 | 0/18 | 0/18 | REJECTED |
| projected5 | 0/18 | 0/18 | REJECTED |
| projected6 | 0/18 | 0/18 | REJECTED |

The gate compares against the pointwise fastest of all 49 retained references,
with a paired 95% lower confidence bound in every one of 18 regression/holdout
groups. No source was tuned after holdout timing. These are finite generic-solver
measurements, not calibrated operation ratios, asymptotic results or a full
index-calculus result. A positive group does not establish a universal gain.

## Complete cold comparison on n24 fresh holdouts

Times are milliseconds per complete cold solve plus validation, pooled across two
fresh seeds and eight random/reverse repetitions in each family. The last ratio is
the retained dispatcher's pooled planted median divided by the displayed arm's
median; it is descriptive and is not the acceptance statistic. All 52 methods
remain visible. The original fixtures and all 192 predecessor cases are retained.

| Method | Planted (ms) | Cross-planted (ms) | Unplanted (ms) | Dispatcher / arm, planted | Correctness |
|---|---:|---:|---:|---:|---|
| search | 80.001708 | 47.048313 | 144.082667 | 0.006 | PASS |
| flat | 293.787812 | 149.567542 | 551.156271 | 0.002 | PASS |
| bucket | 1063.146729 | 522.085438 | 2073.274041 | 0.000 | PASS |
| hybrid | 514.655333 | 217.043645 | 929.048666 | 0.001 | PASS |
| small_flat | 95.769917 | 55.280020 | 168.129188 | 0.005 | PASS |
| word_tail | 78.956625 | 48.750083 | 148.793167 | 0.006 | PASS |
| merge_search | 67.071146 | 41.841020 | 125.967833 | 0.007 | PASS |
| quadratic_state | 27.550375 | 16.624250 | 50.572479 | 0.018 | PASS |
| packed_state | 16.797020 | 10.198187 | 31.142854 | 0.029 | PASS |
| basis_list | 282.821000 | 145.705105 | 473.644334 | 0.002 | PASS |
| basis_wide | 62.981354 | 33.804167 | 101.966166 | 0.008 | PASS |
| tail_list | 132.341771 | 63.002187 | 261.885167 | 0.004 | PASS |
| tail_wide | 40.335729 | 19.494583 | 74.902167 | 0.012 | PASS |
| packed_untraced | 12.085312 | 7.211729 | 22.438917 | 0.040 | PASS |
| affine_sl_basis_list | 261.098875 | 109.415937 | 570.792292 | 0.002 | PASS |
| affine_sl_basis_fast | 22.341000 | 11.336166 | 59.270583 | 0.022 | PASS |
| gray_scalar | 12.095083 | 12.198167 | 20.990667 | 0.040 | PASS |
| gray_simd | 1.786687 | 1.735084 | 3.127355 | 0.271 | PASS |
| packed_gray12_scalar | 14.836375 | 10.110437 | 27.358458 | 0.033 | PASS |
| packed_gray12_simd | 6.599166 | 4.629062 | 12.401625 | 0.073 | PASS |
| packed_gray16_scalar | 11.702291 | 9.840834 | 21.972292 | 0.041 | PASS |
| packed_gray16_simd | 2.139625 | 1.859104 | 3.849729 | 0.226 | PASS |
| gray_delta_scalar | 11.774125 | 11.606979 | 20.919500 | 0.041 | PASS |
| gray_delta_simd | 1.632521 | 1.648708 | 2.863208 | 0.297 | PASS |
| initial_list | 11.813999 | 5.895625 | 21.233896 | 0.041 | PASS |
| initial_simd | 1.626167 | 0.865980 | 2.890875 | 0.298 | PASS |
| packed_gray16_delta_scalar | 11.444187 | 9.950750 | 22.081146 | 0.042 | PASS |
| packed_gray16_delta_simd | 2.035855 | 1.855229 | 3.784687 | 0.238 | PASS |
| fiber_rows | 4.458437 | 4.303229 | 8.958771 | 0.109 | PASS |
| fiber_columns | 4.115854 | 4.017104 | 8.536625 | 0.118 | PASS |
| fiber_simd | 1.554646 | 1.470374 | 3.321729 | 0.312 | PASS |
| fiber_zero_simd | 2.315125 | 2.065333 | 4.517771 | 0.209 | PASS |
| gray_quiet | 1.618208 | 1.631375 | 2.724959 | 0.299 | PASS |
| wide64_quiet | 1.013438 | 1.018417 | 1.798479 | 0.478 | PASS |
| byte_scalar | 4.101854 | 4.145062 | 7.269125 | 0.118 | PASS |
| byte_simd | 1.143271 | 1.085000 | 1.937355 | 0.424 | PASS |
| leaf16_quiet | 1.929312 | 1.649083 | 3.457562 | 0.251 | PASS |
| leaf16_byte_scalar | 4.417916 | 3.702125 | 8.282917 | 0.110 | PASS |
| leaf16_byte_simd | 1.428667 | 1.224229 | 2.714583 | 0.339 | PASS |
| byte_single_quiet | 1.304167 | 1.289146 | 2.367167 | 0.372 | PASS |
| leaf16_single_quiet | 1.688063 | 1.388875 | 3.103000 | 0.287 | PASS |
| byte_planes | 1.149563 | 1.148270 | 2.101666 | 0.422 | PASS |
| leaf16_byte_planes | 1.493667 | 1.286229 | 2.807771 | 0.324 | PASS |
| byte_unrolled | 0.882062 | 0.897208 | 1.582542 | 0.549 | PASS |
| wide64_unrolled | 0.481604 | 0.487708 | 0.847584 | 1.006 | PASS |
| leaf16_byte_unrolled | 1.238458 | 1.030750 | 2.281375 | 0.391 | PASS |
| word16_unrolled | 0.897208 | 0.873667 | 1.550771 | 0.540 | PASS |
| word_dispatch | 0.484604 | 0.485542 | 0.848666 | 1.000 | PASS |
| leaf16_word_unrolled | 1.228167 | 1.086000 | 2.281417 | 0.395 | PASS |
| projected4 | 1.575000 | 1.532229 | 3.701479 | 0.308 | PASS |
| projected5 | 1.586854 | 1.449938 | 2.955208 | 0.305 | PASS |
| projected6 | 1.518646 | 1.398583 | 2.735750 | 0.319 | PASS |

## Mathematical work and limits

Projection removes the quadratic terms internal to the chosen low block, making
the projected subsystem affine. It can also erase constraints. The solver therefore
returns a particular affine answer **and a complete kernel basis**, enumerates the
remaining freedom, and verifies all original equations. The `xy+1` counterexample
is a regression test: a zero projected system must not be accepted at x=y=0.
Caps remain UNKNOWN. The measured suite completed without censored observations.

`SUMMARY.json` records n24 holdout screening, affine elimination and original-check
counts once per fixture, without multiplying by benchmark repetitions. Its ranks
and quotient dimensions are exact structural quantities; they are not predictions
of success or counts of independent conditional constraints.

The campaign took **1972.255 seconds**. Whole-worker peak RSS was
**24,625,152 bytes**, including all methods and common references.
Candidate-specific allocator peaks and calibrated operation costs are unmeasured.
Temporary solver-workspace destruction occurs inside solve time. Fixture generation,
reference preparation, serialization and returned-diagnostic destruction are outside
arm clocks and inside worker process receipts. No phase cost is subtracted.

## Paired gates

| Candidate | Split | Variables | Family | Median reference / candidate | 95% interval | Above 2.0 |
|---|---|---:|---|---:|---|---|
| projected4 | regression | 16 | planted | 0.2425 | [0.2343, 0.2554] | REJECTED |
| projected4 | regression | 16 | cross_planted | 0.2788 | [0.2639, 0.2967] | REJECTED |
| projected4 | regression | 16 | unplanted | 0.2296 | [0.2205, 0.2424] | REJECTED |
| projected4 | regression | 20 | planted | 0.2112 | [0.1995, 0.2217] | REJECTED |
| projected4 | regression | 20 | cross_planted | 0.2577 | [0.2447, 0.2740] | REJECTED |
| projected4 | regression | 20 | unplanted | 0.1883 | [0.1800, 0.1969] | REJECTED |
| projected4 | regression | 24 | planted | 0.2290 | [0.2140, 0.2621] | REJECTED |
| projected4 | regression | 24 | cross_planted | 0.2633 | [0.2461, 0.2696] | REJECTED |
| projected4 | regression | 24 | unplanted | 0.2299 | [0.2205, 0.2582] | REJECTED |
| projected4 | holdout | 16 | planted | 0.2595 | [0.2199, 0.2791] | REJECTED |
| projected4 | holdout | 16 | cross_planted | 0.2405 | [0.2182, 0.2808] | REJECTED |
| projected4 | holdout | 16 | unplanted | 0.2405 | [0.2192, 0.2466] | REJECTED |
| projected4 | holdout | 20 | planted | 0.1463 | [0.1196, 0.2167] | REJECTED |
| projected4 | holdout | 20 | cross_planted | 0.2289 | [0.2075, 0.2380] | REJECTED |
| projected4 | holdout | 20 | unplanted | 0.1754 | [0.1623, 0.1835] | REJECTED |
| projected4 | holdout | 24 | planted | 0.2966 | [0.2564, 0.3171] | REJECTED |
| projected4 | holdout | 24 | cross_planted | 0.3138 | [0.3072, 0.3190] | REJECTED |
| projected4 | holdout | 24 | unplanted | 0.2242 | [0.2176, 0.2944] | REJECTED |
| projected5 | regression | 16 | planted | 0.2354 | [0.2124, 0.2419] | REJECTED |
| projected5 | regression | 16 | cross_planted | 0.2471 | [0.2319, 0.2678] | REJECTED |
| projected5 | regression | 16 | unplanted | 0.2156 | [0.1998, 0.2313] | REJECTED |
| projected5 | regression | 20 | planted | 0.2595 | [0.2472, 0.2722] | REJECTED |
| projected5 | regression | 20 | cross_planted | 0.3226 | [0.2971, 0.3380] | REJECTED |
| projected5 | regression | 20 | unplanted | 0.2404 | [0.2347, 0.2472] | REJECTED |
| projected5 | regression | 24 | planted | 0.2713 | [0.2598, 0.2791] | REJECTED |
| projected5 | regression | 24 | cross_planted | 0.3021 | [0.2713, 0.3187] | REJECTED |
| projected5 | regression | 24 | unplanted | 0.2710 | [0.2640, 0.2789] | REJECTED |
| projected5 | holdout | 16 | planted | 0.2035 | [0.1827, 0.2089] | REJECTED |
| projected5 | holdout | 16 | cross_planted | 0.1855 | [0.1588, 0.2115] | REJECTED |
| projected5 | holdout | 16 | unplanted | 0.2002 | [0.1795, 0.2155] | REJECTED |
| projected5 | holdout | 20 | planted | 0.1846 | [0.1259, 0.2594] | REJECTED |
| projected5 | holdout | 20 | cross_planted | 0.2563 | [0.2060, 0.2632] | REJECTED |
| projected5 | holdout | 20 | unplanted | 0.2255 | [0.2078, 0.2401] | REJECTED |
| projected5 | holdout | 24 | planted | 0.2886 | [0.2680, 0.3123] | REJECTED |
| projected5 | holdout | 24 | cross_planted | 0.3461 | [0.3246, 0.3719] | REJECTED |
| projected5 | holdout | 24 | unplanted | 0.2830 | [0.2748, 0.2881] | REJECTED |
| projected6 | regression | 16 | planted | 0.1255 | [0.1079, 0.1377] | REJECTED |
| projected6 | regression | 16 | cross_planted | 0.1456 | [0.1310, 0.1645] | REJECTED |
| projected6 | regression | 16 | unplanted | 0.0840 | [0.0768, 0.0913] | REJECTED |
| projected6 | regression | 20 | planted | 0.2439 | [0.2257, 0.2546] | REJECTED |
| projected6 | regression | 20 | cross_planted | 0.3163 | [0.3013, 0.3320] | REJECTED |
| projected6 | regression | 20 | unplanted | 0.2142 | [0.2022, 0.2271] | REJECTED |
| projected6 | regression | 24 | planted | 0.2534 | [0.2459, 0.2636] | REJECTED |
| projected6 | regression | 24 | cross_planted | 0.3597 | [0.3168, 0.3729] | REJECTED |
| projected6 | regression | 24 | unplanted | 0.2402 | [0.2341, 0.2713] | REJECTED |
| projected6 | holdout | 16 | planted | 0.0793 | [0.0640, 0.1011] | REJECTED |
| projected6 | holdout | 16 | cross_planted | 0.0623 | [0.0517, 0.0767] | REJECTED |
| projected6 | holdout | 16 | unplanted | 0.0536 | [0.0426, 0.0729] | REJECTED |
| projected6 | holdout | 20 | planted | 0.1675 | [0.0462, 0.2688] | REJECTED |
| projected6 | holdout | 20 | cross_planted | 0.2538 | [0.2388, 0.2616] | REJECTED |
| projected6 | holdout | 20 | unplanted | 0.1709 | [0.1435, 0.1922] | REJECTED |
| projected6 | holdout | 24 | planted | 0.2717 | [0.2425, 0.3389] | REJECTED |
| projected6 | holdout | 24 | cross_planted | 0.3653 | [0.3403, 0.3777] | REJECTED |
| projected6 | holdout | 24 | unplanted | 0.3087 | [0.2803, 0.3422] | REJECTED |

The 69 Rust tests include exhaustive complete affine fibers, every pair of
three-variable quadratics, changing coefficients and rank, direct projection
identities, original-equation recovery and caps, together with all retained tests.
The 15 Python evidence checks replay the result, independently regenerate the fixtures,
verify source custody and reject altered models, work, receipts and capped claims.
These are producer checks, not an external review.

The earlier support-envelope construction and construction-cost diagnostic remain
unchanged. This result implements the previously prospective scan; it does not
convert their negative or inconclusive timing cells into gains. See `README.md`
for the exact domain, fixed selection, accounting and replay contract.

`NEXT_EXPERIMENT.md` records a separate unimplemented hypothesis: capped implicit
determinant filters for affine consistency. It requires exact changing-coefficient
products, Boolean cancellations, rank-deficient controls and charged construction
before any complete-solver comparison. No performance conclusion follows from it.

Production solver, full index-calculus, calibrated-operation and rho costs remain
**null**. The original dramatic-gain goal remains open.
