# Bounded IC round 2

Round 2: retained. Selected challenger: stop5_word; retained/promoted winner: incumbent. 3480/3480 native/profile pairs independently verified by the measured round's frozen checker. Fresh Linux transport replay is the evidence-PR merge gate. All results are limited to the registered synthetic toy panel.

[Provenance, interpretation and reproduction commands](EVIDENCE.md).

Primary metric: one supplied target, after reusable preparation through scalar replay; fixture generation is outside both timed algorithms. Cold instruction and native process costs are supplementary promotion gates. Each table uses one cost unit. Values are equal-cell geometric means of three-process per-point medians, with no target amortization.

Cold IC ratios use the qualified `pairinv` incumbent. The online table names each IC denominator: challengers use the separately qualified `ic_online` role; reference/control rows use `incumbent`. Online rho / variant is computed directly from matched online costs, never by dividing ratios with different IC denominators. `rho` and `rho_online` are the separately qualified cold-instruction and online-time references. Ratios to rho are descriptive. The K-instruction floor applies only to this full-rank collector; it is not a generic IC lower bound.

The class column labels engineering experiments and accounting controls. No asymptotic advance is claimed. Variant names are readable aliases; the machine-readable export retains every canonical candidate, workload and run ID. `stop3` uses the legacy adaptive orbit bound, which is not a universal three-column guarantee. Actual admitted bases and columns below are authoritative.

## Confirmation

Verified 1080/1080 pairs.

### Single-target online time

| Variant | Online ms | IC denominator | / IC reference | Descriptive 95% interval | Online rho / variant | Verified / scheduled | Class |
|---|---|---|---|---|---|---|---|
| incumbent | 0.0300712 | incumbent | 1 | reference | 5.53271 | 216/216 | accounting |
| ic_online | 0.0302927 | incumbent | 1.00737 | [0.995238, 1.01899] | 5.49226 | 216/216 | accounting |
| rho | 0.225603 | incumbent | 7.50231 | [6.79684, 8.33391] | 0.737468 | 216/216 | accounting |
| rho_online | 0.166375 | incumbent | 5.53271 | [5.11118, 6.01516] | 1 | 216/216 | accounting |
| stop5_word | 0.0331297 | ic_online | 1.09365 | [1.02861, 1.15887] | 5.02193 | 216/216 | engineering |

### Complete cold instructions

| Variant | Complete cold Ir | S = Ir / sqrt(r) | / IC reference | Descriptive 95% interval | / cold rho | / K floor | Verified / scheduled | Class |
|---|---|---|---|---|---|---|---|---|
| incumbent | 1.78088e+06 | 2788.01 | 1 | reference | 0.706919 | 222610 | 216/216 | accounting |
| ic_online | 1.80886e+06 | 2831.81 | 1.01571 | [1.00708, 1.0248] | 0.718025 | 226108 | 216/216 | accounting |
| rho | 2.51922e+06 | 3943.89 | 1.41459 | [1.28029, 1.56319] | 1 | not applicable | 216/216 | accounting |
| rho_online | 4.58321e+06 | 7175.12 | 2.57356 | [2.07447, 3.108] | 1.8193 | not applicable | 216/216 | accounting |
| stop5_word | 1.64326e+06 | 2572.57 | 0.922725 | [0.891417, 0.954667] | 0.652292 | 328653 | 216/216 | engineering |

### Complete cold native time

| Variant | Complete cold ms | / IC reference | Descriptive 95% interval | Verified / scheduled | Class |
|---|---|---|---|---|---|
| incumbent | 0.776046 | 1 | reference | 216/216 | accounting |
| ic_online | 0.785151 | 1.01173 | [0.978519, 1.04686] | 216/216 | accounting |
| rho | 0.823813 | 1.06155 | [1.02343, 1.10013] | 216/216 | accounting |
| rho_online | 1.60659 | 2.07023 | [1.98018, 2.17453] | 216/216 | accounting |
| stop5_word | 0.765974 | 0.987022 | [0.959326, 1.01544] | 216/216 | engineering |

### Actual base and matrix sizes

| Variant | Curve cell | Usable points B / folded columns K / final rank |
|---|---|---|
| incumbent | n17a1 | 272 / 8 / 8 |
| incumbent | n19a0 | 304 / 8 / 8 |
| incumbent | n23a0 | 368 / 8 / 8 |
| incumbent | n23a1 | 368 / 8 / 8 |
| incumbent | n29a1 | 464 / 8 / 8 |
| incumbent | n31a0 | 496 / 8 / 8 |
| ic_online | n17a1 | 272 / 8 / 8 |
| ic_online | n19a0 | 304 / 8 / 8 |
| ic_online | n23a0 | 368 / 8 / 8 |
| ic_online | n23a1 | 368 / 8 / 8 |
| ic_online | n29a1 | 464 / 8 / 8 |
| ic_online | n31a0 | 496 / 8 / 8 |
| stop5_word | n17a1 | 170 / 5 / 5 |
| stop5_word | n19a0 | 190 / 5 / 5 |
| stop5_word | n23a0 | 230 / 5 / 5 |
| stop5_word | n23a1 | 230 / 5 / 5 |
| stop5_word | n29a1 | 290 / 5 / 5 |
| stop5_word | n31a0 | 310 / 5 / 5 |

## Replay

Verified 1080/1080 pairs.

### Single-target online time

| Variant | Online ms | IC denominator | / IC reference | Descriptive 95% interval | Online rho / variant | Verified / scheduled | Class |
|---|---|---|---|---|---|---|---|
| incumbent | 0.0301787 | incumbent | 1 | reference | 5.60364 | 216/216 | accounting |
| ic_online | 0.0304443 | incumbent | 1.0088 | [0.994649, 1.0208] | 5.55476 | 216/216 | accounting |
| rho | 0.227286 | incumbent | 7.53133 | [6.83563, 8.34729] | 0.744044 | 216/216 | accounting |
| rho_online | 0.169111 | incumbent | 5.60364 | [5.18522, 6.08137] | 1 | 216/216 | accounting |
| stop5_word | 0.0329592 | ic_online | 1.08261 | [1.01826, 1.14675] | 5.13091 | 216/216 | engineering |

### Complete cold instructions

| Variant | Complete cold Ir | S = Ir / sqrt(r) | / IC reference | Descriptive 95% interval | / cold rho | / K floor | Verified / scheduled | Class |
|---|---|---|---|---|---|---|---|---|
| incumbent | 1.78088e+06 | 2788.01 | 1 | reference | 0.70692 | 222610 | 216/216 | accounting |
| ic_online | 1.80886e+06 | 2831.81 | 1.01571 | [1.00708, 1.0248] | 0.718027 | 226107 | 216/216 | accounting |
| rho | 2.51921e+06 | 3943.88 | 1.41459 | [1.28029, 1.5632] | 1 | not applicable | 216/216 | accounting |
| rho_online | 4.58322e+06 | 7175.13 | 2.57357 | [2.07451, 3.10801] | 1.81931 | not applicable | 216/216 | accounting |
| stop5_word | 1.64326e+06 | 2572.57 | 0.922726 | [0.891419, 0.954669] | 0.652294 | 328653 | 216/216 | engineering |

### Complete cold native time

| Variant | Complete cold ms | / IC reference | Descriptive 95% interval | Verified / scheduled | Class |
|---|---|---|---|---|---|
| incumbent | 0.782402 | 1 | reference | 216/216 | accounting |
| ic_online | 0.807039 | 1.03149 | [0.994054, 1.06489] | 216/216 | accounting |
| rho | 0.838329 | 1.07148 | [1.03267, 1.11189] | 216/216 | accounting |
| rho_online | 1.6424 | 2.09918 | [1.99526, 2.18742] | 216/216 | accounting |
| stop5_word | 0.785639 | 1.00414 | [0.979479, 1.0295] | 216/216 | engineering |

### Actual base and matrix sizes

| Variant | Curve cell | Usable points B / folded columns K / final rank |
|---|---|---|
| incumbent | n17a1 | 272 / 8 / 8 |
| incumbent | n19a0 | 304 / 8 / 8 |
| incumbent | n23a0 | 368 / 8 / 8 |
| incumbent | n23a1 | 368 / 8 / 8 |
| incumbent | n29a1 | 464 / 8 / 8 |
| incumbent | n31a0 | 496 / 8 / 8 |
| ic_online | n17a1 | 272 / 8 / 8 |
| ic_online | n19a0 | 304 / 8 / 8 |
| ic_online | n23a0 | 368 / 8 / 8 |
| ic_online | n23a1 | 368 / 8 / 8 |
| ic_online | n29a1 | 464 / 8 / 8 |
| ic_online | n31a0 | 496 / 8 / 8 |
| stop5_word | n17a1 | 170 / 5 / 5 |
| stop5_word | n19a0 | 190 / 5 / 5 |
| stop5_word | n23a0 | 230 / 5 / 5 |
| stop5_word | n23a1 | 230 / 5 / 5 |
| stop5_word | n29a1 | 290 / 5 / 5 |
| stop5_word | n31a0 | 310 / 5 / 5 |

## Development

Verified 630/630 pairs.

### Single-target online time

| Variant | Online ms | IC denominator | / IC reference | Descriptive 95% interval | Online rho / variant | Verified / scheduled | Class |
|---|---|---|---|---|---|---|---|
| incumbent | 0.0324391 | incumbent | 1 | reference | 5.4419 | 45/45 | accounting |
| batch4_stop6_word | 0.031925 | ic_online | 0.988777 | [0.958135, 1.01936] | 5.52952 | 45/45 | engineering |
| compatibility | 0.0322854 | ic_online | 0.999939 | [0.952806, 1.07447] | 5.4678 | 45/45 | accounting |
| cover_table | 0.0308678 | ic_online | 0.956032 | [0.894381, 1.01261] | 5.71891 | 45/45 | engineering |
| half_table | 0.0307428 | ic_online | 0.95216 | [0.910193, 0.979171] | 5.74217 | 45/45 | engineering |
| ic_online | 0.0322874 | incumbent | 0.995324 | [0.9508, 1.03375] | 5.46746 | 45/45 | accounting |
| rho | 0.240285 | incumbent | 7.40729 | [6.20784, 8.86932] | 0.734668 | 45/45 | accounting |
| rho_online | 0.17653 | incumbent | 5.4419 | [4.60466, 6.34739] | 1 | 45/45 | accounting |
| stop5_word | 0.0368943 | ic_online | 1.14268 | [1.04381, 1.25609] | 4.78475 | 45/45 | engineering |
| stop6_half | 0.0343731 | ic_online | 1.0646 | [0.969192, 1.24539] | 5.1357 | 45/45 | engineering |
| stop6_word | 0.0348506 | ic_online | 1.07939 | [0.984671, 1.24676] | 5.06534 | 45/45 | engineering |
| stop6_word_half | 0.0349149 | ic_online | 1.08138 | [0.979282, 1.25134] | 5.05602 | 45/45 | engineering |
| stop7_half | 0.0326492 | ic_online | 1.01121 | [0.951769, 1.08414] | 5.40687 | 45/45 | engineering |
| word_half | 0.0310254 | ic_online | 0.960914 | [0.910474, 1.00331] | 5.68986 | 45/45 | engineering |

### Complete cold instructions

| Variant | Complete cold Ir | S = Ir / sqrt(r) | / IC reference | Descriptive 95% interval | / cold rho | / K floor | Verified / scheduled | Class |
|---|---|---|---|---|---|---|---|---|
| incumbent | 1.84007e+06 | 2297.31 | 1 | reference | 0.695632 | 230009 | 45/45 | accounting |
| batch4_stop6_word | 1.84334e+06 | 2301.39 | 1.00178 | [1.00092, 1.00258] | 0.696868 | 230417 | 45/45 | engineering |
| compatibility | 1.84043e+06 | 2297.76 | 1.00019 | [0.999721, 1.0007] | 0.695767 | 230053 | 45/45 | accounting |
| cover_table | 2.92164e+06 | 3647.65 | 1.58779 | [1.26462, 1.9492] | 1.10452 | 365205 | 45/45 | engineering |
| half_table | 2.32126e+06 | 2898.08 | 1.26151 | [1.07197, 1.45119] | 0.877544 | 290157 | 45/45 | engineering |
| ic_online | 1.87399e+06 | 2339.66 | 1.01843 | [1.00803, 1.02973] | 0.708454 | 234248 | 45/45 | accounting |
| rho | 2.64518e+06 | 3302.48 | 1.43754 | [1.17677, 1.72102] | 1 | not applicable | 45/45 | accounting |
| rho_online | 4.3996e+06 | 5492.87 | 2.391 | [1.77221, 3.0715] | 1.66326 | not applicable | 45/45 | accounting |
| stop5_word | 1.7647e+06 | 2203.21 | 0.959039 | [0.889063, 1.07556] | 0.667138 | 352940 | 45/45 | engineering |
| stop6_half | 1.99384e+06 | 2489.3 | 1.08357 | [0.991588, 1.14722] | 0.753764 | 332307 | 45/45 | engineering |
| stop6_word | 1.76314e+06 | 2201.27 | 0.958194 | [0.896201, 1.03968] | 0.66655 | 293857 | 45/45 | engineering |
| stop6_word_half | 1.99171e+06 | 2486.64 | 1.08241 | [0.99081, 1.14589] | 0.75296 | 331952 | 45/45 | engineering |
| stop7_half | 2.20854e+06 | 2757.35 | 1.20025 | [1.04028, 1.35899] | 0.834931 | 315506 | 45/45 | engineering |
| word_half | 2.31756e+06 | 2893.46 | 1.2595 | [1.07083, 1.44848] | 0.876147 | 289695 | 45/45 | engineering |

### Complete cold native time

| Variant | Complete cold ms | / IC reference | Descriptive 95% interval | Verified / scheduled | Class |
|---|---|---|---|---|---|
| incumbent | 0.78059 | 1 | reference | 45/45 | accounting |
| batch4_stop6_word | 0.786094 | 1.00705 | [0.9669, 1.05491] | 45/45 | engineering |
| compatibility | 0.788452 | 1.01007 | [0.954645, 1.05953] | 45/45 | accounting |
| cover_table | 0.884694 | 1.13336 | [1.05871, 1.23197] | 45/45 | engineering |
| half_table | 0.801723 | 1.02707 | [0.943319, 1.11723] | 45/45 | engineering |
| ic_online | 0.801967 | 1.02739 | [0.957903, 1.09589] | 45/45 | accounting |
| rho | 0.844018 | 1.08126 | [1.02359, 1.14318] | 45/45 | accounting |
| rho_online | 1.56977 | 2.011 | [1.86215, 2.16501] | 45/45 | accounting |
| stop5_word | 0.764284 | 0.97911 | [0.900617, 1.06604] | 45/45 | engineering |
| stop6_half | 0.796603 | 1.02051 | [0.951836, 1.0937] | 45/45 | engineering |
| stop6_word | 0.771569 | 0.988442 | [0.926275, 1.06167] | 45/45 | engineering |
| stop6_word_half | 0.806693 | 1.03344 | [0.97319, 1.10112] | 45/45 | engineering |
| stop7_half | 0.771435 | 0.988271 | [0.917037, 1.07152] | 45/45 | engineering |
| word_half | 0.82124 | 1.05207 | [0.982673, 1.13446] | 45/45 | engineering |

Retained portfolio:

- `half_table`: single-target online-time leader.
- `stop6_word`: complete instruction-cost leader.
- `stop5_word`: complete native-time leader.
- `batch4_stop6_word`: non-dominated implementation family.
- `stop6_half`: non-dominated implementation family.
- `word_half`: predeclared exploration slot.

### Actual base and matrix sizes

| Variant | Curve cell | Usable points B / folded columns K / final rank |
|---|---|---|
| incumbent | n17a1 | 272 / 8 / 8 |
| incumbent | n19a0 | 304 / 8 / 8 |
| incumbent | n23a0 | 368 / 8 / 8 |
| incumbent | n23a1 | 368 / 8 / 8 |
| incumbent | n31a0 | 496 / 8 / 8 |
| batch4_stop6_word | n17a1 | 272 / 8 / 8 |
| batch4_stop6_word | n19a0 | 304 / 8 / 8 |
| batch4_stop6_word | n23a0 | 368 / 8 / 8 |
| batch4_stop6_word | n23a1 | 368 / 8 / 8 |
| batch4_stop6_word | n31a0 | 496 / 8 / 8 |
| compatibility | n17a1 | 272 / 8 / 8 |
| compatibility | n19a0 | 304 / 8 / 8 |
| compatibility | n23a0 | 368 / 8 / 8 |
| compatibility | n23a1 | 368 / 8 / 8 |
| compatibility | n31a0 | 496 / 8 / 8 |
| cover_table | n17a1 | 272 / 8 / 8 |
| cover_table | n19a0 | 304 / 8 / 8 |
| cover_table | n23a0 | 368 / 8 / 8 |
| cover_table | n23a1 | 368 / 8 / 8 |
| cover_table | n31a0 | 496 / 8 / 8 |
| half_table | n17a1 | 272 / 8 / 8 |
| half_table | n19a0 | 304 / 8 / 8 |
| half_table | n23a0 | 368 / 8 / 8 |
| half_table | n23a1 | 368 / 8 / 8 |
| half_table | n31a0 | 496 / 8 / 8 |
| ic_online | n17a1 | 272 / 8 / 8 |
| ic_online | n19a0 | 304 / 8 / 8 |
| ic_online | n23a0 | 368 / 8 / 8 |
| ic_online | n23a1 | 368 / 8 / 8 |
| ic_online | n31a0 | 496 / 8 / 8 |
| stop5_word | n17a1 | 170 / 5 / 5 |
| stop5_word | n19a0 | 190 / 5 / 5 |
| stop5_word | n23a0 | 230 / 5 / 5 |
| stop5_word | n23a1 | 230 / 5 / 5 |
| stop5_word | n31a0 | 310 / 5 / 5 |
| stop6_half | n17a1 | 204 / 6 / 6 |
| stop6_half | n19a0 | 228 / 6 / 6 |
| stop6_half | n23a0 | 276 / 6 / 6 |
| stop6_half | n23a1 | 276 / 6 / 6 |
| stop6_half | n31a0 | 372 / 6 / 6 |
| stop6_word | n17a1 | 204 / 6 / 6 |
| stop6_word | n19a0 | 228 / 6 / 6 |
| stop6_word | n23a0 | 276 / 6 / 6 |
| stop6_word | n23a1 | 276 / 6 / 6 |
| stop6_word | n31a0 | 372 / 6 / 6 |
| stop6_word_half | n17a1 | 204 / 6 / 6 |
| stop6_word_half | n19a0 | 228 / 6 / 6 |
| stop6_word_half | n23a0 | 276 / 6 / 6 |
| stop6_word_half | n23a1 | 276 / 6 / 6 |
| stop6_word_half | n31a0 | 372 / 6 / 6 |
| stop7_half | n17a1 | 238 / 7 / 7 |
| stop7_half | n19a0 | 266 / 7 / 7 |
| stop7_half | n23a0 | 322 / 7 / 7 |
| stop7_half | n23a1 | 322 / 7 / 7 |
| stop7_half | n31a0 | 434 / 7 / 7 |
| word_half | n17a1 | 272 / 8 / 8 |
| word_half | n19a0 | 304 / 8 / 8 |
| word_half | n23a0 | 368 / 8 / 8 |
| word_half | n23a1 | 368 / 8 / 8 |
| word_half | n31a0 | 496 / 8 / 8 |

## Smoke

Verified 210/210 pairs.

### Single-target online time

| Variant | Online ms | IC denominator | / IC reference | Descriptive 95% interval | Online rho / variant | Verified / scheduled | Class |
|---|---|---|---|---|---|---|---|
| incumbent | 0.0343879 | incumbent | 1 | reference | 5.77961 | 15/15 | accounting |
| batch4_stop6_word | 0.0340717 | ic_online | 0.969299 | [0.932079, 0.998657] | 5.83326 | 15/15 | engineering |
| compatibility | 0.0341205 | ic_online | 0.970689 | [0.921854, 1.00575] | 5.8249 | 15/15 | accounting |
| cover_table | 0.0336985 | ic_online | 0.958684 | [0.939401, 0.972397] | 5.89785 | 15/15 | engineering |
| half_table | 0.0334549 | ic_online | 0.951754 | [0.904468, 0.981786] | 5.94079 | 15/15 | engineering |
| ic_online | 0.0351508 | incumbent | 1.02219 | [0.978404, 1.08059] | 5.65417 | 15/15 | accounting |
| rho | 0.264357 | incumbent | 7.68749 | [6.55676, 9.20841] | 0.75182 | 15/15 | accounting |
| rho_online | 0.198749 | incumbent | 5.77961 | [5.148, 6.79648] | 1 | 15/15 | accounting |
| stop5_word | 0.0365239 | ic_online | 1.03906 | [0.950518, 1.20629] | 5.4416 | 15/15 | engineering |
| stop6_half | 0.0331454 | ic_online | 0.942949 | [0.903755, 0.983814] | 5.99626 | 15/15 | engineering |
| stop6_word | 0.0345409 | ic_online | 0.982649 | [0.909099, 1.06425] | 5.75401 | 15/15 | engineering |
| stop6_word_half | 0.033814 | ic_online | 0.961968 | [0.911507, 1.02496] | 5.87771 | 15/15 | engineering |
| stop7_half | 0.0328168 | ic_online | 0.9336 | [0.902179, 0.966443] | 6.05631 | 15/15 | engineering |
| word_half | 0.0345434 | ic_online | 0.98272 | [0.907491, 1.08555] | 5.75359 | 15/15 | engineering |

### Complete cold instructions

| Variant | Complete cold Ir | S = Ir / sqrt(r) | / IC reference | Descriptive 95% interval | / cold rho | / K floor | Verified / scheduled | Class |
|---|---|---|---|---|---|---|---|---|
| incumbent | 1.95591e+06 | 2441.95 | 1 | reference | 0.668756 | 244489 | 15/15 | accounting |
| batch4_stop6_word | 1.95893e+06 | 2445.72 | 1.00154 | [1.00033, 1.00278] | 0.669789 | 244867 | 15/15 | engineering |
| compatibility | 1.95601e+06 | 2442.06 | 1.00005 | [0.999089, 1.00101] | 0.668789 | 244501 | 15/15 | accounting |
| cover_table | 2.8749e+06 | 3589.29 | 1.46985 | [1.12193, 1.85791] | 0.982971 | 359362 | 15/15 | engineering |
| half_table | 2.22654e+06 | 2779.83 | 1.13837 | [0.900282, 1.38731] | 0.761289 | 278318 | 15/15 | engineering |
| ic_online | 1.99761e+06 | 2494.01 | 1.02132 | [1.00994, 1.03359] | 0.683015 | 249702 | 15/15 | accounting |
| rho | 2.9247e+06 | 3651.47 | 1.49531 | [1.41065, 1.57653] | 1 | not applicable | 15/15 | accounting |
| rho_online | 4.65417e+06 | 5810.7 | 2.37954 | [1.93471, 2.79046] | 1.59133 | not applicable | 15/15 | accounting |
| stop5_word | 1.64788e+06 | 2057.37 | 0.842514 | [0.739185, 0.937249] | 0.563436 | 329577 | 15/15 | engineering |
| stop6_half | 1.78125e+06 | 2223.88 | 0.910698 | [0.717783, 1.07938] | 0.609035 | 296875 | 15/15 | engineering |
| stop6_word | 1.60664e+06 | 2005.88 | 0.821426 | [0.695266, 0.925538] | 0.549334 | 267773 | 15/15 | engineering |
| stop6_word_half | 1.77942e+06 | 2221.59 | 0.909764 | [0.717054, 1.07808] | 0.60841 | 296570 | 15/15 | engineering |
| stop7_half | 2.0637e+06 | 2576.51 | 1.05511 | [0.808529, 1.30343] | 0.705609 | 294814 | 15/15 | engineering |
| word_half | 2.22329e+06 | 2775.76 | 1.1367 | [0.899087, 1.38493] | 0.760175 | 277911 | 15/15 | engineering |

### Complete cold native time

| Variant | Complete cold ms | / IC reference | Descriptive 95% interval | Verified / scheduled | Class |
|---|---|---|---|---|---|
| incumbent | 0.883762 | 1 | reference | 15/15 | accounting |
| batch4_stop6_word | 0.806024 | 0.912038 | [0.824216, 0.989386] | 15/15 | engineering |
| compatibility | 0.769606 | 0.87083 | [0.822997, 0.926405] | 15/15 | accounting |
| cover_table | 0.87437 | 0.989373 | [0.914658, 1.10576] | 15/15 | engineering |
| half_table | 0.836014 | 0.945973 | [0.841427, 1.0603] | 15/15 | engineering |
| ic_online | 0.834291 | 0.944023 | [0.867752, 1.027] | 15/15 | accounting |
| rho | 0.834498 | 0.944257 | [0.907044, 0.982996] | 15/15 | accounting |
| rho_online | 1.64785 | 1.86459 | [1.79109, 1.9345] | 15/15 | accounting |
| stop5_word | 0.751737 | 0.850611 | [0.784426, 0.914588] | 15/15 | engineering |
| stop6_half | 0.779652 | 0.882198 | [0.817026, 0.961222] | 15/15 | engineering |
| stop6_word | 0.795277 | 0.899877 | [0.864115, 0.936088] | 15/15 | engineering |
| stop6_word_half | 0.785345 | 0.888639 | [0.819207, 0.978242] | 15/15 | engineering |
| stop7_half | 0.783953 | 0.887064 | [0.814203, 0.956984] | 15/15 | engineering |
| word_half | 0.893904 | 1.01148 | [0.909177, 1.12712] | 15/15 | engineering |

### Actual base and matrix sizes

| Variant | Curve cell | Usable points B / folded columns K / final rank |
|---|---|---|
| incumbent | n17a1 | 272 / 8 / 8 |
| incumbent | n19a0 | 304 / 8 / 8 |
| incumbent | n23a0 | 368 / 8 / 8 |
| incumbent | n23a1 | 368 / 8 / 8 |
| incumbent | n31a0 | 496 / 8 / 8 |
| batch4_stop6_word | n17a1 | 272 / 8 / 8 |
| batch4_stop6_word | n19a0 | 304 / 8 / 8 |
| batch4_stop6_word | n23a0 | 368 / 8 / 8 |
| batch4_stop6_word | n23a1 | 368 / 8 / 8 |
| batch4_stop6_word | n31a0 | 496 / 8 / 8 |
| compatibility | n17a1 | 272 / 8 / 8 |
| compatibility | n19a0 | 304 / 8 / 8 |
| compatibility | n23a0 | 368 / 8 / 8 |
| compatibility | n23a1 | 368 / 8 / 8 |
| compatibility | n31a0 | 496 / 8 / 8 |
| cover_table | n17a1 | 272 / 8 / 8 |
| cover_table | n19a0 | 304 / 8 / 8 |
| cover_table | n23a0 | 368 / 8 / 8 |
| cover_table | n23a1 | 368 / 8 / 8 |
| cover_table | n31a0 | 496 / 8 / 8 |
| half_table | n17a1 | 272 / 8 / 8 |
| half_table | n19a0 | 304 / 8 / 8 |
| half_table | n23a0 | 368 / 8 / 8 |
| half_table | n23a1 | 368 / 8 / 8 |
| half_table | n31a0 | 496 / 8 / 8 |
| ic_online | n17a1 | 272 / 8 / 8 |
| ic_online | n19a0 | 304 / 8 / 8 |
| ic_online | n23a0 | 368 / 8 / 8 |
| ic_online | n23a1 | 368 / 8 / 8 |
| ic_online | n31a0 | 496 / 8 / 8 |
| stop5_word | n17a1 | 170 / 5 / 5 |
| stop5_word | n19a0 | 190 / 5 / 5 |
| stop5_word | n23a0 | 230 / 5 / 5 |
| stop5_word | n23a1 | 230 / 5 / 5 |
| stop5_word | n31a0 | 310 / 5 / 5 |
| stop6_half | n17a1 | 204 / 6 / 6 |
| stop6_half | n19a0 | 228 / 6 / 6 |
| stop6_half | n23a0 | 276 / 6 / 6 |
| stop6_half | n23a1 | 276 / 6 / 6 |
| stop6_half | n31a0 | 372 / 6 / 6 |
| stop6_word | n17a1 | 204 / 6 / 6 |
| stop6_word | n19a0 | 228 / 6 / 6 |
| stop6_word | n23a0 | 276 / 6 / 6 |
| stop6_word | n23a1 | 276 / 6 / 6 |
| stop6_word | n31a0 | 372 / 6 / 6 |
| stop6_word_half | n17a1 | 204 / 6 / 6 |
| stop6_word_half | n19a0 | 228 / 6 / 6 |
| stop6_word_half | n23a0 | 276 / 6 / 6 |
| stop6_word_half | n23a1 | 276 / 6 / 6 |
| stop6_word_half | n31a0 | 372 / 6 / 6 |
| stop7_half | n17a1 | 238 / 7 / 7 |
| stop7_half | n19a0 | 266 / 7 / 7 |
| stop7_half | n23a0 | 322 / 7 / 7 |
| stop7_half | n23a1 | 322 / 7 / 7 |
| stop7_half | n31a0 | 434 / 7 / 7 |
| word_half | n17a1 | 272 / 8 / 8 |
| word_half | n19a0 | 304 / 8 / 8 |
| word_half | n23a0 | 368 / 8 / 8 |
| word_half | n23a1 | 368 / 8 / 8 |
| word_half | n31a0 | 496 / 8 / 8 |

## Selection

Verified 450/450 pairs.

### Single-target online time

| Variant | Online ms | IC denominator | / IC reference | Descriptive 95% interval | Online rho / variant | Verified / scheduled | Class |
|---|---|---|---|---|---|---|---|
| incumbent | 0.032541 | incumbent | 1 | reference | 5.55286 | 45/45 | accounting |
| batch4_stop6_word | 0.0321056 | ic_online | 0.981566 | [0.961966, 1.00298] | 5.62817 | 45/45 | engineering |
| half_table | 0.0314828 | ic_online | 0.962525 | [0.941738, 0.982505] | 5.73951 | 45/45 | engineering |
| ic_online | 0.0327085 | incumbent | 1.00515 | [0.993509, 1.01776] | 5.52442 | 45/45 | accounting |
| rho | 0.245161 | incumbent | 7.5339 | [6.45513, 8.72947] | 0.737049 | 45/45 | accounting |
| rho_online | 0.180695 | incumbent | 5.55286 | [4.85652, 6.28377] | 1 | 45/45 | accounting |
| stop5_word | 0.0334425 | ic_online | 1.02244 | [0.951767, 1.13196] | 5.40316 | 45/45 | engineering |
| stop6_half | 0.0330198 | ic_online | 1.00952 | [0.961085, 1.11163] | 5.47234 | 45/45 | engineering |
| stop6_word | 0.0330347 | ic_online | 1.00997 | [0.955874, 1.10989] | 5.46988 | 45/45 | engineering |
| word_half | 0.0316654 | ic_online | 0.96811 | [0.939361, 1.00195] | 5.7064 | 45/45 | engineering |

### Complete cold instructions

| Variant | Complete cold Ir | S = Ir / sqrt(r) | / IC reference | Descriptive 95% interval | / cold rho | / K floor | Verified / scheduled | Class |
|---|---|---|---|---|---|---|---|---|
| incumbent | 1.75294e+06 | 2188.54 | 1 | reference | 0.651111 | 219118 | 45/45 | accounting |
| batch4_stop6_word | 1.75608e+06 | 2192.45 | 1.00179 | [1.00074, 1.0026] | 0.652276 | 219510 | 45/45 | engineering |
| half_table | 2.24448e+06 | 2802.22 | 1.2804 | [1.12456, 1.43251] | 0.833686 | 280560 | 45/45 | engineering |
| ic_online | 1.78023e+06 | 2222.61 | 1.01557 | [1.00748, 1.02473] | 0.661248 | 222529 | 45/45 | accounting |
| rho | 2.69223e+06 | 3361.24 | 1.53584 | [1.35252, 1.68147] | 1 | not applicable | 45/45 | accounting |
| rho_online | 4.44488e+06 | 5549.4 | 2.53566 | [2.04978, 3.03767] | 1.651 | not applicable | 45/45 | accounting |
| stop5_word | 1.61737e+06 | 2019.28 | 0.922662 | [0.840074, 1.0685] | 0.600756 | 323475 | 45/45 | engineering |
| stop6_half | 1.84908e+06 | 2308.57 | 1.05484 | [0.982786, 1.11766] | 0.686821 | 308180 | 45/45 | engineering |
| stop6_word | 1.72975e+06 | 2159.58 | 0.986767 | [0.902847, 1.10735] | 0.642495 | 288291 | 45/45 | engineering |
| word_half | 2.24079e+06 | 2797.61 | 1.2783 | [1.12293, 1.42972] | 0.832315 | 280098 | 45/45 | engineering |

### Complete cold native time

| Variant | Complete cold ms | / IC reference | Descriptive 95% interval | Verified / scheduled | Class |
|---|---|---|---|---|---|
| incumbent | 0.807803 | 1 | reference | 45/45 | accounting |
| batch4_stop6_word | 0.776076 | 0.960725 | [0.903317, 1.02945] | 45/45 | engineering |
| half_table | 0.837177 | 1.03636 | [0.991041, 1.07899] | 45/45 | engineering |
| ic_online | 0.842292 | 1.04269 | [0.994276, 1.09398] | 45/45 | accounting |
| rho | 0.874824 | 1.08297 | [1.03424, 1.13639] | 45/45 | accounting |
| rho_online | 1.65626 | 2.05032 | [1.9013, 2.16493] | 45/45 | accounting |
| stop5_word | 0.751534 | 0.930343 | [0.852299, 0.999863] | 45/45 | engineering |
| stop6_half | 0.797849 | 0.987678 | [0.937866, 1.032] | 45/45 | engineering |
| stop6_word | 0.774784 | 0.959125 | [0.88991, 1.03793] | 45/45 | engineering |
| word_half | 0.87816 | 1.0871 | [1.04196, 1.13625] | 45/45 | engineering |

### Actual base and matrix sizes

| Variant | Curve cell | Usable points B / folded columns K / final rank |
|---|---|---|
| incumbent | n17a1 | 272 / 8 / 8 |
| incumbent | n19a0 | 304 / 8 / 8 |
| incumbent | n23a0 | 368 / 8 / 8 |
| incumbent | n23a1 | 368 / 8 / 8 |
| incumbent | n31a0 | 496 / 8 / 8 |
| batch4_stop6_word | n17a1 | 272 / 8 / 8 |
| batch4_stop6_word | n19a0 | 304 / 8 / 8 |
| batch4_stop6_word | n23a0 | 368 / 8 / 8 |
| batch4_stop6_word | n23a1 | 368 / 8 / 8 |
| batch4_stop6_word | n31a0 | 496 / 8 / 8 |
| half_table | n17a1 | 272 / 8 / 8 |
| half_table | n19a0 | 304 / 8 / 8 |
| half_table | n23a0 | 368 / 8 / 8 |
| half_table | n23a1 | 368 / 8 / 8 |
| half_table | n31a0 | 496 / 8 / 8 |
| ic_online | n17a1 | 272 / 8 / 8 |
| ic_online | n19a0 | 304 / 8 / 8 |
| ic_online | n23a0 | 368 / 8 / 8 |
| ic_online | n23a1 | 368 / 8 / 8 |
| ic_online | n31a0 | 496 / 8 / 8 |
| stop5_word | n17a1 | 170 / 5 / 5 |
| stop5_word | n19a0 | 190 / 5 / 5 |
| stop5_word | n23a0 | 230 / 5 / 5 |
| stop5_word | n23a1 | 230 / 5 / 5 |
| stop5_word | n31a0 | 310 / 5 / 5 |
| stop6_half | n17a1 | 204 / 6 / 6 |
| stop6_half | n19a0 | 228 / 6 / 6 |
| stop6_half | n23a0 | 276 / 6 / 6 |
| stop6_half | n23a1 | 276 / 6 / 6 |
| stop6_half | n31a0 | 372 / 6 / 6 |
| stop6_word | n17a1 | 204 / 6 / 6 |
| stop6_word | n19a0 | 228 / 6 / 6 |
| stop6_word | n23a0 | 276 / 6 / 6 |
| stop6_word | n23a1 | 276 / 6 / 6 |
| stop6_word | n31a0 | 372 / 6 / 6 |
| word_half | n17a1 | 272 / 8 / 8 |
| word_half | n19a0 | 304 / 8 / 8 |
| word_half | n23a0 | 368 / 8 / 8 |
| word_half | n23a1 | 368 / 8 / 8 |
| word_half | n31a0 | 496 / 8 / 8 |

## Aa

Verified 30/30 pairs.

### Single-target online time

| Variant | Online ms | IC denominator | / IC reference | Descriptive 95% interval | Online rho / variant | Verified / scheduled | Class |
|---|---|---|---|---|---|---|---|
| incumbent | 0.0304482 | incumbent | 1 | reference | unknown | 15/15 | accounting |
| aa_control | 0.0300838 | incumbent | 0.988032 | [0.957747, 1.01113] | unknown | 15/15 | accounting |

### Complete cold instructions

| Variant | Complete cold Ir | S = Ir / sqrt(r) | / IC reference | Descriptive 95% interval | / cold rho | / K floor | Verified / scheduled | Class |
|---|---|---|---|---|---|---|---|---|
| incumbent | 1.68199e+06 | 2099.95 | 1 | reference | unknown | 210248 | 15/15 | accounting |
| aa_control | 1.68199e+06 | 2099.95 | 0.999999 | [0.999997, 1] | unknown | 210248 | 15/15 | accounting |

### Complete cold native time

| Variant | Complete cold ms | / IC reference | Descriptive 95% interval | Verified / scheduled | Class |
|---|---|---|---|---|---|
| incumbent | 0.875529 | 1 | reference | 15/15 | accounting |
| aa_control | 0.820886 | 0.937588 | [0.813701, 1.07041] | 15/15 | accounting |

### Actual base and matrix sizes

| Variant | Curve cell | Usable points B / folded columns K / final rank |
|---|---|---|
| incumbent | n17a1 | 272 / 8 / 8 |
| incumbent | n19a0 | 304 / 8 / 8 |
| incumbent | n23a0 | 368 / 8 / 8 |
| incumbent | n23a1 | 368 / 8 / 8 |
| incumbent | n31a0 | 496 / 8 / 8 |
| aa_control | n17a1 | 272 / 8 / 8 |
| aa_control | n19a0 | 304 / 8 / 8 |
| aa_control | n23a0 | 368 / 8 / 8 |
| aa_control | n23a1 | 368 / 8 / 8 |
| aa_control | n31a0 | 496 / 8 / 8 |

## Confirmation rule and scope

The selected challenger must pass both final stages: cold Ir and native ratios <=0.8, online ratio <=1, all three predeclared one-sided familywise upper bounds <1, and each metric's largest per-cell ratio <=1.1. The nominal familywise rule allocates alpha across three attempts, two stages and three metrics. Bootstrap coverage is approximate, and replay repeats the same targets in new processes. Ordinary 95% table intervals are descriptive; they do not replace this gate.

### Confirmation promotion evidence

| Metric | Candidate / reference | Familywise upper | Largest cell ratio |
|---|---|---|---|
| online_ns | 1.09365 | 1.18502 | 1.50242 |
| instructions | 0.922725 | 0.967882 | 1.1279 |
| cold_ns | 0.987022 | 1.027 | 1.05604 |

### Replay promotion evidence

| Metric | Candidate / reference | Familywise upper | Largest cell ratio |
|---|---|---|---|
| online_ns | 1.08261 | 1.17291 | 1.46751 |
| instructions | 0.922726 | 0.967883 | 1.1279 |
| cold_ns | 1.00414 | 1.04003 | 1.06157 |

Decision reasons:

- No challenger passed every confirmation and replay threshold.

Complete run keys, admitted B/column counts, exclusive instruction ledgers, online phase clocks, peak process RSS and collection statistics are retained in `RUNS.csv` and `RESULTS.json`; raw certificates, sources, binaries and profiles are in the archive. Base-only memory was not instrumented and remains unknown. Yield intervals resample whole target/walk clusters, never individual dependent queries; few-target diagnostic intervals have weak coverage.

F4/F5/SAT and large sparse LA backends were not executed in this round. The tested mechanisms are pair-table PDP policies, Frobenius-orbit base construction and scalar-field Gaussian row kernels. No result here establishes a production-curve crossover or a globally fastest IC method.
