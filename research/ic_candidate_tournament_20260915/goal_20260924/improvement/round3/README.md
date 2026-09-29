# Bounded IC round 3

Round 3: retained. Selected challenger: stop7_word; retained/promoted winner: incumbent. 3480/3480 native/profile pairs independently verified by the measured round's frozen checker. Fresh Linux transport replay is the evidence-PR merge gate. All results are limited to the registered synthetic toy panel.

[Provenance, interpretation and reproduction commands](EVIDENCE.md).

Primary metric: one supplied target, after reusable preparation through scalar replay; fixture generation is outside both timed algorithms. Cold instruction and native process costs are supplementary promotion gates. Each table uses one cost unit. Values are equal-cell geometric means of three-process per-point medians, with no target amortization.

Cold IC ratios use the qualified `pairinv` incumbent. The online table names each IC denominator: challengers use the separately qualified `ic_online` role; reference/control rows use `incumbent`. Online rho / variant is computed directly from matched online costs, never by dividing ratios with different IC denominators. `rho` and `rho_online` are the separately qualified cold-instruction and online-time references. Ratios to rho are descriptive. The K-instruction floor applies only to this full-rank collector; it is not a generic IC lower bound.

The class column labels engineering experiments and accounting controls. No asymptotic advance is claimed. Variant names are readable aliases; the machine-readable export retains every canonical candidate, workload and run ID. `stop3` uses the legacy adaptive orbit bound, which is not a universal three-column guarantee. Actual admitted bases and columns below are authoritative.

## Confirmation

Verified 1080/1080 pairs.

### Single-target online time

| Variant | Online ms | IC denominator | / IC reference | Descriptive 95% interval | Online rho / variant | Verified / scheduled | Class |
|---|---|---|---|---|---|---|---|
| incumbent | 0.0339263 | incumbent | 1 | reference | 7.44684 | 216/216 | accounting |
| ic_online | 0.0348801 | incumbent | 1.02811 | [0.982522, 1.0798] | 7.2432 | 216/216 | accounting |
| rho | 0.286506 | incumbent | 8.44497 | [7.78213, 9.1975] | 0.881808 | 216/216 | accounting |
| rho_online | 0.252643 | incumbent | 7.44684 | [6.94635, 8.02028] | 1 | 216/216 | accounting |
| stop7_word | 0.03395 | ic_online | 0.973335 | [0.932577, 1.0161] | 7.44163 | 216/216 | engineering |

### Complete cold instructions

| Variant | Complete cold Ir | S = Ir / sqrt(r) | / IC reference | Descriptive 95% interval | / cold rho | / K floor | Verified / scheduled | Class |
|---|---|---|---|---|---|---|---|---|
| incumbent | 1.80249e+06 | 2821.83 | 1 | reference | 0.72135 | 225311 | 216/216 | accounting |
| ic_online | 1.83093e+06 | 2866.36 | 1.01578 | [1.00701, 1.02497] | 0.732733 | 228866 | 216/216 | accounting |
| rho | 2.49877e+06 | 3911.88 | 1.38629 | [1.2411, 1.53301] | 1 | not applicable | 216/216 | accounting |
| rho_online | 4.56307e+06 | 7143.59 | 2.53154 | [1.99963, 3.10109] | 1.82613 | not applicable | 216/216 | accounting |
| stop7_word | 1.78479e+06 | 2794.13 | 0.990182 | [0.963649, 1.02062] | 0.714267 | 254970 | 216/216 | engineering |

### Complete cold native time

| Variant | Complete cold ms | / IC reference | Descriptive 95% interval | Verified / scheduled | Class |
|---|---|---|---|---|---|
| incumbent | 0.888366 | 1 | reference | 216/216 | accounting |
| ic_online | 0.885258 | 0.996502 | [0.971982, 1.01839] | 216/216 | accounting |
| rho | 0.859978 | 0.968046 | [0.915412, 1.02749] | 216/216 | accounting |
| rho_online | 1.90207 | 2.14109 | [1.96908, 2.31144] | 216/216 | accounting |
| stop7_word | 0.875659 | 0.985697 | [0.962282, 1.0102] | 216/216 | engineering |

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
| stop7_word | n17a1 | 238 / 7 / 7 |
| stop7_word | n19a0 | 266 / 7 / 7 |
| stop7_word | n23a0 | 322 / 7 / 7 |
| stop7_word | n23a1 | 322 / 7 / 7 |
| stop7_word | n29a1 | 406 / 7 / 7 |
| stop7_word | n31a0 | 434 / 7 / 7 |

## Replay

Verified 1080/1080 pairs.

### Single-target online time

| Variant | Online ms | IC denominator | / IC reference | Descriptive 95% interval | Online rho / variant | Verified / scheduled | Class |
|---|---|---|---|---|---|---|---|
| incumbent | 0.0342533 | incumbent | 1 | reference | 7.49015 | 216/216 | accounting |
| ic_online | 0.0347676 | incumbent | 1.01502 | [0.967392, 1.06221] | 7.37934 | 216/216 | accounting |
| rho | 0.281739 | incumbent | 8.22517 | [7.3824, 9.14448] | 0.910637 | 216/216 | accounting |
| rho_online | 0.256562 | incumbent | 7.49015 | [6.87004, 8.22436] | 1 | 216/216 | accounting |
| stop7_word | 0.0340288 | ic_online | 0.97875 | [0.936818, 1.02233] | 7.53956 | 216/216 | engineering |

### Complete cold instructions

| Variant | Complete cold Ir | S = Ir / sqrt(r) | / IC reference | Descriptive 95% interval | / cold rho | / K floor | Verified / scheduled | Class |
|---|---|---|---|---|---|---|---|---|
| incumbent | 1.80248e+06 | 2821.83 | 1 | reference | 0.721351 | 225311 | 216/216 | accounting |
| ic_online | 1.83093e+06 | 2866.36 | 1.01578 | [1.00702, 1.02497] | 0.732735 | 228866 | 216/216 | accounting |
| rho | 2.49876e+06 | 3911.87 | 1.38629 | [1.2411, 1.53299] | 1 | not applicable | 216/216 | accounting |
| rho_online | 4.56309e+06 | 7143.62 | 2.53156 | [1.99964, 3.10111] | 1.82614 | not applicable | 216/216 | accounting |
| stop7_word | 1.78479e+06 | 2794.13 | 0.990183 | [0.963649, 1.02062] | 0.71427 | 254970 | 216/216 | engineering |

### Complete cold native time

| Variant | Complete cold ms | / IC reference | Descriptive 95% interval | Verified / scheduled | Class |
|---|---|---|---|---|---|
| incumbent | 0.889031 | 1 | reference | 216/216 | accounting |
| ic_online | 0.904588 | 1.0175 | [0.992254, 1.04324] | 216/216 | accounting |
| rho | 0.856718 | 0.963654 | [0.924107, 1.00131] | 216/216 | accounting |
| rho_online | 1.89573 | 2.13236 | [1.98661, 2.2674] | 216/216 | accounting |
| stop7_word | 0.871766 | 0.98058 | [0.957679, 1.00429] | 216/216 | engineering |

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
| stop7_word | n17a1 | 238 / 7 / 7 |
| stop7_word | n19a0 | 266 / 7 / 7 |
| stop7_word | n23a0 | 322 / 7 / 7 |
| stop7_word | n23a1 | 322 / 7 / 7 |
| stop7_word | n29a1 | 406 / 7 / 7 |
| stop7_word | n31a0 | 434 / 7 / 7 |

## Development

Verified 630/630 pairs.

### Single-target online time

| Variant | Online ms | IC denominator | / IC reference | Descriptive 95% interval | Online rho / variant | Verified / scheduled | Class |
|---|---|---|---|---|---|---|---|
| incumbent | 0.0385771 | incumbent | 1 | reference | 6.94774 | 45/45 | accounting |
| batch2_stop4_word_half | 0.0394742 | ic_online | 1.0987 | [0.921522, 1.44079] | 6.78985 | 45/45 | engineering |
| ic_online | 0.035928 | incumbent | 0.931328 | [0.867625, 1.00847] | 7.46004 | 45/45 | accounting |
| rho | 0.307113 | incumbent | 7.96102 | [6.86011, 9.91793] | 0.872721 | 45/45 | accounting |
| rho_online | 0.268024 | incumbent | 6.94774 | [5.92489, 8.67992] | 1 | 45/45 | accounting |
| stop3_word_half | 0.0462711 | ic_online | 1.28789 | [0.86859, 2.4968] | 5.79247 | 45/45 | engineering |
| stop4_word_cover | 0.0377981 | ic_online | 1.05205 | [0.859708, 1.4178] | 7.09094 | 45/45 | engineering |
| stop4_word_half | 0.0392488 | ic_online | 1.09243 | [0.929816, 1.42591] | 6.82884 | 45/45 | engineering |
| stop5_half | 0.0368815 | ic_online | 1.02654 | [0.921422, 1.15678] | 7.26717 | 45/45 | engineering |
| stop5_word_half | 0.0361499 | ic_online | 1.00618 | [0.895144, 1.15229] | 7.41425 | 45/45 | engineering |
| stop6_bounded_half | 0.0348761 | ic_online | 0.970723 | [0.888893, 1.10354] | 7.68503 | 45/45 | engineering |
| stop6_word_cover | 0.0353002 | ic_online | 0.982528 | [0.856516, 1.17785] | 7.5927 | 45/45 | engineering |
| stop7_word | 0.0361961 | ic_online | 1.00746 | [0.919718, 1.14042] | 7.40478 | 45/45 | engineering |
| word_cover | 0.0338846 | ic_online | 0.943125 | [0.880026, 1.02493] | 7.90991 | 45/45 | engineering |

### Complete cold instructions

| Variant | Complete cold Ir | S = Ir / sqrt(r) | / IC reference | Descriptive 95% interval | / cold rho | / K floor | Verified / scheduled | Class |
|---|---|---|---|---|---|---|---|---|
| incumbent | 1.80161e+06 | 2249.3 | 1 | reference | 0.675047 | 225201 | 45/45 | accounting |
| batch2_stop4_word_half | 1.73491e+06 | 2166.02 | 0.962977 | [0.816092, 1.2722] | 0.650055 | 433727 | 45/45 | engineering |
| ic_online | 1.83208e+06 | 2287.34 | 1.01691 | [1.00892, 1.0248] | 0.686465 | 229010 | 45/45 | accounting |
| rho | 2.66886e+06 | 3332.06 | 1.48138 | [1.32453, 1.65186] | 1 | not applicable | 45/45 | accounting |
| rho_online | 4.46939e+06 | 5580.01 | 2.48078 | [2.15523, 2.83222] | 1.67464 | not applicable | 45/45 | accounting |
| stop3_word_half | 1.69919e+06 | 2121.43 | 0.943152 | [0.753007, 1.39372] | 0.636672 | 566397 | 45/45 | engineering |
| stop4_word_cover | 1.68469e+06 | 2103.33 | 0.935106 | [0.835785, 1.09658] | 0.63124 | 421174 | 45/45 | engineering |
| stop4_word_half | 1.71935e+06 | 2146.59 | 0.95434 | [0.807554, 1.26345] | 0.644224 | 429837 | 45/45 | engineering |
| stop5_half | 1.79106e+06 | 2236.13 | 0.994147 | [0.854393, 1.1671] | 0.671096 | 358212 | 45/45 | engineering |
| stop5_word_half | 1.78963e+06 | 2234.34 | 0.993351 | [0.853671, 1.16643] | 0.670558 | 357926 | 45/45 | engineering |
| stop6_bounded_half | 1.87069e+06 | 2335.55 | 1.03835 | [0.947918, 1.12218] | 0.700932 | 311782 | 45/45 | engineering |
| stop6_word_cover | 2.14198e+06 | 2674.24 | 1.18892 | [1.05644, 1.31977] | 0.80258 | 356996 | 45/45 | engineering |
| stop7_word | 1.7553e+06 | 2191.48 | 0.974295 | [0.918202, 1.0335] | 0.657695 | 250757 | 45/45 | engineering |
| word_cover | 2.86971e+06 | 3582.82 | 1.59286 | [1.35931, 1.85432] | 1.07526 | 358714 | 45/45 | engineering |

### Complete cold native time

| Variant | Complete cold ms | / IC reference | Descriptive 95% interval | Verified / scheduled | Class |
|---|---|---|---|---|---|
| incumbent | 0.879348 | 1 | reference | 45/45 | accounting |
| batch2_stop4_word_half | 0.819806 | 0.932288 | [0.837927, 1.07888] | 45/45 | engineering |
| ic_online | 0.864572 | 0.983196 | [0.929643, 1.05316] | 45/45 | accounting |
| rho | 0.874583 | 0.994581 | [0.921921, 1.07544] | 45/45 | accounting |
| rho_online | 1.89854 | 2.15903 | [2.00355, 2.34797] | 45/45 | accounting |
| stop3_word_half | 0.816267 | 0.928264 | [0.805095, 1.17803] | 45/45 | engineering |
| stop4_word_cover | 0.814317 | 0.926046 | [0.859071, 1.02391] | 45/45 | engineering |
| stop4_word_half | 0.832169 | 0.946347 | [0.829622, 1.10855] | 45/45 | engineering |
| stop5_half | 0.84849 | 0.964907 | [0.863484, 1.08911] | 45/45 | engineering |
| stop5_word_half | 0.84102 | 0.956412 | [0.856136, 1.08061] | 45/45 | engineering |
| stop6_bounded_half | 0.839366 | 0.954532 | [0.88299, 1.03332] | 45/45 | engineering |
| stop6_word_cover | 0.88988 | 1.01198 | [0.955936, 1.06995] | 45/45 | engineering |
| stop7_word | 0.858469 | 0.976256 | [0.916536, 1.0548] | 45/45 | engineering |
| word_cover | 1.01002 | 1.1486 | [1.05044, 1.24527] | 45/45 | engineering |

Retained portfolio:

- `word_cover`: single-target online-time leader.
- `stop4_word_cover`: complete instruction-cost leader.
- `stop3_word_half`: cell specialist: cold instructions n17a1.
- `stop4_word_half`: non-dominated implementation family.
- `batch2_stop4_word_half`: non-dominated implementation family.
- `stop7_word`: predeclared exploration slot.

### Actual base and matrix sizes

| Variant | Curve cell | Usable points B / folded columns K / final rank |
|---|---|---|
| incumbent | n17a1 | 272 / 8 / 8 |
| incumbent | n19a0 | 304 / 8 / 8 |
| incumbent | n23a0 | 368 / 8 / 8 |
| incumbent | n23a1 | 368 / 8 / 8 |
| incumbent | n31a0 | 496 / 8 / 8 |
| batch2_stop4_word_half | n17a1 | 136 / 4 / 4 |
| batch2_stop4_word_half | n19a0 | 152 / 4 / 4 |
| batch2_stop4_word_half | n23a0 | 184 / 4 / 4 |
| batch2_stop4_word_half | n23a1 | 184 / 4 / 4 |
| batch2_stop4_word_half | n31a0 | 248 / 4 / 4 |
| ic_online | n17a1 | 272 / 8 / 8 |
| ic_online | n19a0 | 304 / 8 / 8 |
| ic_online | n23a0 | 368 / 8 / 8 |
| ic_online | n23a1 | 368 / 8 / 8 |
| ic_online | n31a0 | 496 / 8 / 8 |
| stop3_word_half | n17a1 | 102 / 3 / 3 |
| stop3_word_half | n19a0 | 114 / 3 / 3 |
| stop3_word_half | n23a0 | 138 / 3 / 3 |
| stop3_word_half | n23a1 | 138 / 3 / 3 |
| stop3_word_half | n31a0 | 186 / 3 / 3 |
| stop4_word_cover | n17a1 | 136 / 4 / 4 |
| stop4_word_cover | n19a0 | 152 / 4 / 4 |
| stop4_word_cover | n23a0 | 184 / 4 / 4 |
| stop4_word_cover | n23a1 | 184 / 4 / 4 |
| stop4_word_cover | n31a0 | 248 / 4 / 4 |
| stop4_word_half | n17a1 | 136 / 4 / 4 |
| stop4_word_half | n19a0 | 152 / 4 / 4 |
| stop4_word_half | n23a0 | 184 / 4 / 4 |
| stop4_word_half | n23a1 | 184 / 4 / 4 |
| stop4_word_half | n31a0 | 248 / 4 / 4 |
| stop5_half | n17a1 | 170 / 5 / 5 |
| stop5_half | n19a0 | 190 / 5 / 5 |
| stop5_half | n23a0 | 230 / 5 / 5 |
| stop5_half | n23a1 | 230 / 5 / 5 |
| stop5_half | n31a0 | 310 / 5 / 5 |
| stop5_word_half | n17a1 | 170 / 5 / 5 |
| stop5_word_half | n19a0 | 190 / 5 / 5 |
| stop5_word_half | n23a0 | 230 / 5 / 5 |
| stop5_word_half | n23a1 | 230 / 5 / 5 |
| stop5_word_half | n31a0 | 310 / 5 / 5 |
| stop6_bounded_half | n17a1 | 204 / 6 / 6 |
| stop6_bounded_half | n19a0 | 228 / 6 / 6 |
| stop6_bounded_half | n23a0 | 276 / 6 / 6 |
| stop6_bounded_half | n23a1 | 276 / 6 / 6 |
| stop6_bounded_half | n31a0 | 372 / 6 / 6 |
| stop6_word_cover | n17a1 | 204 / 6 / 6 |
| stop6_word_cover | n19a0 | 228 / 6 / 6 |
| stop6_word_cover | n23a0 | 276 / 6 / 6 |
| stop6_word_cover | n23a1 | 276 / 6 / 6 |
| stop6_word_cover | n31a0 | 372 / 6 / 6 |
| stop7_word | n17a1 | 238 / 7 / 7 |
| stop7_word | n19a0 | 266 / 7 / 7 |
| stop7_word | n23a0 | 322 / 7 / 7 |
| stop7_word | n23a1 | 322 / 7 / 7 |
| stop7_word | n31a0 | 434 / 7 / 7 |
| word_cover | n17a1 | 272 / 8 / 8 |
| word_cover | n19a0 | 304 / 8 / 8 |
| word_cover | n23a0 | 368 / 8 / 8 |
| word_cover | n23a1 | 368 / 8 / 8 |
| word_cover | n31a0 | 496 / 8 / 8 |

## Smoke

Verified 210/210 pairs.

### Single-target online time

| Variant | Online ms | IC denominator | / IC reference | Descriptive 95% interval | Online rho / variant | Verified / scheduled | Class |
|---|---|---|---|---|---|---|---|
| incumbent | 0.0471979 | incumbent | 1 | reference | 5.68062 | 15/15 | accounting |
| batch2_stop4_word_half | 0.0382659 | ic_online | 0.823242 | [0.733182, 0.944586] | 7.00658 | 15/15 | engineering |
| ic_online | 0.0464819 | incumbent | 0.984831 | [0.861925, 1.14768] | 5.76811 | 15/15 | accounting |
| rho | 0.306729 | incumbent | 6.4988 | [3.70671, 9.2523] | 0.874103 | 15/15 | accounting |
| rho_online | 0.268113 | incumbent | 5.68062 | [3.2452, 7.86101] | 1 | 15/15 | accounting |
| stop3_word_half | 0.0386649 | ic_online | 0.831827 | [0.68968, 0.982107] | 6.93427 | 15/15 | engineering |
| stop4_word_cover | 0.0378888 | ic_online | 0.81513 | [0.714919, 0.989613] | 7.07631 | 15/15 | engineering |
| stop4_word_half | 0.0402608 | ic_online | 0.86616 | [0.78551, 0.954756] | 6.65941 | 15/15 | engineering |
| stop5_half | 0.0405233 | ic_online | 0.871808 | [0.735675, 1.01081] | 6.61626 | 15/15 | engineering |
| stop5_word_half | 0.0388704 | ic_online | 0.836248 | [0.718074, 1.03284] | 6.89761 | 15/15 | engineering |
| stop6_bounded_half | 0.0386008 | ic_online | 0.830447 | [0.718027, 0.947932] | 6.94579 | 15/15 | engineering |
| stop6_word_cover | 0.0378828 | ic_online | 0.815001 | [0.666526, 1.0143] | 7.07743 | 15/15 | engineering |
| stop7_word | 0.0409453 | ic_online | 0.880887 | [0.791374, 1.00072] | 6.54808 | 15/15 | engineering |
| word_cover | 0.0402694 | ic_online | 0.866344 | [0.747914, 0.955202] | 6.65799 | 15/15 | engineering |

### Complete cold instructions

| Variant | Complete cold Ir | S = Ir / sqrt(r) | / IC reference | Descriptive 95% interval | / cold rho | / K floor | Verified / scheduled | Class |
|---|---|---|---|---|---|---|---|---|
| incumbent | 1.85252e+06 | 2312.86 | 1 | reference | 0.687295 | 231565 | 15/15 | accounting |
| batch2_stop4_word_half | 1.73965e+06 | 2171.95 | 0.939076 | [0.796341, 1.10824] | 0.645422 | 434913 | 15/15 | engineering |
| ic_online | 1.88716e+06 | 2356.1 | 1.0187 | [1.00926, 1.02869] | 0.700146 | 235894 | 15/15 | accounting |
| rho | 2.69538e+06 | 3365.16 | 1.45498 | [1.09588, 1.78785] | 1 | not applicable | 15/15 | accounting |
| rho_online | 4.44911e+06 | 5554.69 | 2.40166 | [1.72904, 2.99334] | 1.65065 | not applicable | 15/15 | accounting |
| stop3_word_half | 2.01557e+06 | 2516.42 | 1.08802 | [0.879223, 1.3329] | 0.747787 | 671856 | 15/15 | engineering |
| stop4_word_cover | 1.82197e+06 | 2274.72 | 0.983511 | [0.846609, 1.14931] | 0.675962 | 455493 | 15/15 | engineering |
| stop4_word_half | 1.72865e+06 | 2158.2 | 0.933134 | [0.797793, 1.10131] | 0.641338 | 432162 | 15/15 | engineering |
| stop5_half | 1.68805e+06 | 2107.52 | 0.911218 | [0.799333, 1.02804] | 0.626275 | 337609 | 15/15 | engineering |
| stop5_word_half | 1.68663e+06 | 2105.75 | 0.910454 | [0.79884, 1.02702] | 0.62575 | 337326 | 15/15 | engineering |
| stop6_bounded_half | 1.88797e+06 | 2357.12 | 1.01914 | [0.928693, 1.08933] | 0.700448 | 314662 | 15/15 | engineering |
| stop6_word_cover | 2.17511e+06 | 2715.61 | 1.17414 | [1.04958, 1.31348] | 0.806979 | 362519 | 15/15 | engineering |
| stop7_word | 1.83854e+06 | 2295.41 | 0.992455 | [0.962621, 1.02734] | 0.682109 | 262649 | 15/15 | engineering |
| word_cover | 2.90913e+06 | 3632.03 | 1.57037 | [1.3082, 1.8648] | 1.0793 | 363641 | 15/15 | engineering |

### Complete cold native time

| Variant | Complete cold ms | / IC reference | Descriptive 95% interval | Verified / scheduled | Class |
|---|---|---|---|---|---|
| incumbent | 0.926555 | 1 | reference | 15/15 | accounting |
| batch2_stop4_word_half | 0.797957 | 0.861208 | [0.801228, 0.921624] | 15/15 | engineering |
| ic_online | 0.901561 | 0.973024 | [0.894155, 1.09305] | 15/15 | accounting |
| rho | 0.875588 | 0.944993 | [0.814792, 1.06446] | 15/15 | accounting |
| rho_online | 1.90139 | 2.05211 | [1.74631, 2.35645] | 15/15 | accounting |
| stop3_word_half | 0.871412 | 0.940485 | [0.886429, 0.992239] | 15/15 | engineering |
| stop4_word_cover | 0.813774 | 0.878279 | [0.802844, 0.955157] | 15/15 | engineering |
| stop4_word_half | 0.817891 | 0.882722 | [0.855644, 0.918773] | 15/15 | engineering |
| stop5_half | 0.812514 | 0.876919 | [0.79573, 0.974211] | 15/15 | engineering |
| stop5_word_half | 0.819389 | 0.884338 | [0.781968, 0.996074] | 15/15 | engineering |
| stop6_bounded_half | 0.823936 | 0.889246 | [0.810976, 0.972651] | 15/15 | engineering |
| stop6_word_cover | 0.912822 | 0.985178 | [0.872509, 1.1101] | 15/15 | engineering |
| stop7_word | 0.838382 | 0.904837 | [0.833541, 0.982232] | 15/15 | engineering |
| word_cover | 0.983447 | 1.0614 | [0.963368, 1.15041] | 15/15 | engineering |

### Actual base and matrix sizes

| Variant | Curve cell | Usable points B / folded columns K / final rank |
|---|---|---|
| incumbent | n17a1 | 272 / 8 / 8 |
| incumbent | n19a0 | 304 / 8 / 8 |
| incumbent | n23a0 | 368 / 8 / 8 |
| incumbent | n23a1 | 368 / 8 / 8 |
| incumbent | n31a0 | 496 / 8 / 8 |
| batch2_stop4_word_half | n17a1 | 136 / 4 / 4 |
| batch2_stop4_word_half | n19a0 | 152 / 4 / 4 |
| batch2_stop4_word_half | n23a0 | 184 / 4 / 4 |
| batch2_stop4_word_half | n23a1 | 184 / 4 / 4 |
| batch2_stop4_word_half | n31a0 | 248 / 4 / 4 |
| ic_online | n17a1 | 272 / 8 / 8 |
| ic_online | n19a0 | 304 / 8 / 8 |
| ic_online | n23a0 | 368 / 8 / 8 |
| ic_online | n23a1 | 368 / 8 / 8 |
| ic_online | n31a0 | 496 / 8 / 8 |
| stop3_word_half | n17a1 | 102 / 3 / 3 |
| stop3_word_half | n19a0 | 114 / 3 / 3 |
| stop3_word_half | n23a0 | 138 / 3 / 3 |
| stop3_word_half | n23a1 | 138 / 3 / 3 |
| stop3_word_half | n31a0 | 186 / 3 / 3 |
| stop4_word_cover | n17a1 | 136 / 4 / 4 |
| stop4_word_cover | n19a0 | 152 / 4 / 4 |
| stop4_word_cover | n23a0 | 184 / 4 / 4 |
| stop4_word_cover | n23a1 | 184 / 4 / 4 |
| stop4_word_cover | n31a0 | 248 / 4 / 4 |
| stop4_word_half | n17a1 | 136 / 4 / 4 |
| stop4_word_half | n19a0 | 152 / 4 / 4 |
| stop4_word_half | n23a0 | 184 / 4 / 4 |
| stop4_word_half | n23a1 | 184 / 4 / 4 |
| stop4_word_half | n31a0 | 248 / 4 / 4 |
| stop5_half | n17a1 | 170 / 5 / 5 |
| stop5_half | n19a0 | 190 / 5 / 5 |
| stop5_half | n23a0 | 230 / 5 / 5 |
| stop5_half | n23a1 | 230 / 5 / 5 |
| stop5_half | n31a0 | 310 / 5 / 5 |
| stop5_word_half | n17a1 | 170 / 5 / 5 |
| stop5_word_half | n19a0 | 190 / 5 / 5 |
| stop5_word_half | n23a0 | 230 / 5 / 5 |
| stop5_word_half | n23a1 | 230 / 5 / 5 |
| stop5_word_half | n31a0 | 310 / 5 / 5 |
| stop6_bounded_half | n17a1 | 204 / 6 / 6 |
| stop6_bounded_half | n19a0 | 228 / 6 / 6 |
| stop6_bounded_half | n23a0 | 276 / 6 / 6 |
| stop6_bounded_half | n23a1 | 276 / 6 / 6 |
| stop6_bounded_half | n31a0 | 372 / 6 / 6 |
| stop6_word_cover | n17a1 | 204 / 6 / 6 |
| stop6_word_cover | n19a0 | 228 / 6 / 6 |
| stop6_word_cover | n23a0 | 276 / 6 / 6 |
| stop6_word_cover | n23a1 | 276 / 6 / 6 |
| stop6_word_cover | n31a0 | 372 / 6 / 6 |
| stop7_word | n17a1 | 238 / 7 / 7 |
| stop7_word | n19a0 | 266 / 7 / 7 |
| stop7_word | n23a0 | 322 / 7 / 7 |
| stop7_word | n23a1 | 322 / 7 / 7 |
| stop7_word | n31a0 | 434 / 7 / 7 |
| word_cover | n17a1 | 272 / 8 / 8 |
| word_cover | n19a0 | 304 / 8 / 8 |
| word_cover | n23a0 | 368 / 8 / 8 |
| word_cover | n23a1 | 368 / 8 / 8 |
| word_cover | n31a0 | 496 / 8 / 8 |

## Selection

Verified 450/450 pairs.

### Single-target online time

| Variant | Online ms | IC denominator | / IC reference | Descriptive 95% interval | Online rho / variant | Verified / scheduled | Class |
|---|---|---|---|---|---|---|---|
| incumbent | 0.0403974 | incumbent | 1 | reference | 6.681 | 45/45 | accounting |
| batch2_stop4_word_half | 0.0487036 | ic_online | 1.15915 | [0.806892, 1.75376] | 5.54159 | 45/45 | engineering |
| ic_online | 0.0420165 | incumbent | 1.04008 | [0.946188, 1.1466] | 6.42356 | 45/45 | accounting |
| rho | 0.290987 | incumbent | 7.2031 | [5.80541, 8.70796] | 0.927517 | 45/45 | accounting |
| rho_online | 0.269895 | incumbent | 6.681 | [5.23493, 8.2041] | 1 | 45/45 | accounting |
| stop3_word_half | 0.0487007 | ic_online | 1.15908 | [0.871092, 1.77326] | 5.54192 | 45/45 | engineering |
| stop4_word_cover | 0.0462867 | ic_online | 1.10163 | [0.767603, 1.68454] | 5.83095 | 45/45 | engineering |
| stop4_word_half | 0.0471902 | ic_online | 1.12314 | [0.798865, 1.69673] | 5.71931 | 45/45 | engineering |
| stop7_word | 0.0377137 | ic_online | 0.897593 | [0.823444, 0.989101] | 7.15643 | 45/45 | engineering |
| word_cover | 0.0347709 | ic_online | 0.827555 | [0.661756, 0.954295] | 7.7621 | 45/45 | engineering |

### Complete cold instructions

| Variant | Complete cold Ir | S = Ir / sqrt(r) | / IC reference | Descriptive 95% interval | / cold rho | / K floor | Verified / scheduled | Class |
|---|---|---|---|---|---|---|---|---|
| incumbent | 1.809e+06 | 2258.52 | 1 | reference | 0.706566 | 226124 | 45/45 | accounting |
| batch2_stop4_word_half | 1.82981e+06 | 2284.5 | 1.0115 | [0.861523, 1.25262] | 0.714694 | 457451 | 45/45 | engineering |
| ic_online | 1.84261e+06 | 2300.49 | 1.01858 | [1.00918, 1.0277] | 0.719695 | 230326 | 45/45 | accounting |
| rho | 2.56026e+06 | 3196.47 | 1.4153 | [1.25304, 1.61246] | 1 | not applicable | 45/45 | accounting |
| rho_online | 4.36748e+06 | 5452.77 | 2.41431 | [2.01265, 2.90425] | 1.70587 | not applicable | 45/45 | accounting |
| stop3_word_half | 1.84746e+06 | 2306.54 | 1.02126 | [0.822499, 1.36123] | 0.721588 | 615819 | 45/45 | engineering |
| stop4_word_cover | 1.89077e+06 | 2360.62 | 1.04521 | [0.900589, 1.26806] | 0.738507 | 472693 | 45/45 | engineering |
| stop4_word_half | 1.81071e+06 | 2260.66 | 1.00095 | [0.844349, 1.24743] | 0.707235 | 452677 | 45/45 | engineering |
| stop7_word | 1.76572e+06 | 2204.49 | 0.976077 | [0.908098, 1.02766] | 0.689663 | 252246 | 45/45 | engineering |
| word_cover | 2.88914e+06 | 3607.07 | 1.5971 | [1.361, 1.87544] | 1.12845 | 361142 | 45/45 | engineering |

### Complete cold native time

| Variant | Complete cold ms | / IC reference | Descriptive 95% interval | Verified / scheduled | Class |
|---|---|---|---|---|---|
| incumbent | 0.877279 | 1 | reference | 45/45 | accounting |
| batch2_stop4_word_half | 0.828381 | 0.944262 | [0.848629, 1.07376] | 45/45 | engineering |
| ic_online | 0.885736 | 1.00964 | [0.947445, 1.07437] | 45/45 | accounting |
| rho | 0.859795 | 0.980071 | [0.914178, 1.04512] | 45/45 | accounting |
| rho_online | 1.89826 | 2.1638 | [2.02623, 2.34818] | 45/45 | accounting |
| stop3_word_half | 0.837455 | 0.954606 | [0.806633, 1.171] | 45/45 | engineering |
| stop4_word_cover | 0.853137 | 0.972481 | [0.880068, 1.08584] | 45/45 | engineering |
| stop4_word_half | 0.824117 | 0.939401 | [0.837428, 1.07206] | 45/45 | engineering |
| stop7_word | 0.874938 | 0.997331 | [0.928945, 1.06854] | 45/45 | engineering |
| word_cover | 1.02076 | 1.16355 | [1.08785, 1.22205] | 45/45 | engineering |

### Actual base and matrix sizes

| Variant | Curve cell | Usable points B / folded columns K / final rank |
|---|---|---|
| incumbent | n17a1 | 272 / 8 / 8 |
| incumbent | n19a0 | 304 / 8 / 8 |
| incumbent | n23a0 | 368 / 8 / 8 |
| incumbent | n23a1 | 368 / 8 / 8 |
| incumbent | n31a0 | 496 / 8 / 8 |
| batch2_stop4_word_half | n17a1 | 136 / 4 / 4 |
| batch2_stop4_word_half | n19a0 | 152 / 4 / 4 |
| batch2_stop4_word_half | n23a0 | 184 / 4 / 4 |
| batch2_stop4_word_half | n23a1 | 184 / 4 / 4 |
| batch2_stop4_word_half | n31a0 | 248 / 4 / 4 |
| ic_online | n17a1 | 272 / 8 / 8 |
| ic_online | n19a0 | 304 / 8 / 8 |
| ic_online | n23a0 | 368 / 8 / 8 |
| ic_online | n23a1 | 368 / 8 / 8 |
| ic_online | n31a0 | 496 / 8 / 8 |
| stop3_word_half | n17a1 | 102 / 3 / 3 |
| stop3_word_half | n19a0 | 114 / 3 / 3 |
| stop3_word_half | n23a0 | 138 / 3 / 3 |
| stop3_word_half | n23a1 | 138 / 3 / 3 |
| stop3_word_half | n31a0 | 186 / 3 / 3 |
| stop4_word_cover | n17a1 | 136 / 4 / 4 |
| stop4_word_cover | n19a0 | 152 / 4 / 4 |
| stop4_word_cover | n23a0 | 184 / 4 / 4 |
| stop4_word_cover | n23a1 | 184 / 4 / 4 |
| stop4_word_cover | n31a0 | 248 / 4 / 4 |
| stop4_word_half | n17a1 | 136 / 4 / 4 |
| stop4_word_half | n19a0 | 152 / 4 / 4 |
| stop4_word_half | n23a0 | 184 / 4 / 4 |
| stop4_word_half | n23a1 | 184 / 4 / 4 |
| stop4_word_half | n31a0 | 248 / 4 / 4 |
| stop7_word | n17a1 | 238 / 7 / 7 |
| stop7_word | n19a0 | 266 / 7 / 7 |
| stop7_word | n23a0 | 322 / 7 / 7 |
| stop7_word | n23a1 | 322 / 7 / 7 |
| stop7_word | n31a0 | 434 / 7 / 7 |
| word_cover | n17a1 | 272 / 8 / 8 |
| word_cover | n19a0 | 304 / 8 / 8 |
| word_cover | n23a0 | 368 / 8 / 8 |
| word_cover | n23a1 | 368 / 8 / 8 |
| word_cover | n31a0 | 496 / 8 / 8 |

## Aa

Verified 30/30 pairs.

### Single-target online time

| Variant | Online ms | IC denominator | / IC reference | Descriptive 95% interval | Online rho / variant | Verified / scheduled | Class |
|---|---|---|---|---|---|---|---|
| incumbent | 0.0343183 | incumbent | 1 | reference | unknown | 15/15 | accounting |
| aa_control | 0.034089 | incumbent | 0.993318 | [0.972655, 1.01679] | unknown | 15/15 | accounting |

### Complete cold instructions

| Variant | Complete cold Ir | S = Ir / sqrt(r) | / IC reference | Descriptive 95% interval | / cold rho | / K floor | Verified / scheduled | Class |
|---|---|---|---|---|---|---|---|---|
| incumbent | 1.61006e+06 | 2010.15 | 1 | reference | unknown | 201257 | 15/15 | accounting |
| aa_control | 1.61007e+06 | 2010.16 | 1 | [1, 1.00001] | unknown | 201258 | 15/15 | accounting |

### Complete cold native time

| Variant | Complete cold ms | / IC reference | Descriptive 95% interval | Verified / scheduled | Class |
|---|---|---|---|---|---|
| incumbent | 0.798407 | 1 | reference | 15/15 | accounting |
| aa_control | 0.825431 | 1.03385 | [0.996806, 1.07129] | 15/15 | accounting |

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
| online_ns | 0.973335 | 1.03389 | 1.05058 |
| instructions | 0.990182 | 1.03495 | 1.03869 |
| cold_ns | 0.985697 | 1.02056 | 1.04795 |

### Replay promotion evidence

| Metric | Candidate / reference | Familywise upper | Largest cell ratio |
|---|---|---|---|
| online_ns | 0.97875 | 1.04037 | 1.0866 |
| instructions | 0.990183 | 1.03495 | 1.03869 |
| cold_ns | 0.98058 | 1.01437 | 1.01528 |

Decision reasons:

- No challenger passed every confirmation and replay threshold.

Complete run keys, admitted B/column counts, exclusive instruction ledgers, online phase clocks, peak process RSS and collection statistics are retained in `RUNS.csv` and `RESULTS.json`; raw certificates, sources, binaries and profiles are in the archive. Base-only memory was not instrumented and remains unknown. Yield intervals resample whole target/walk clusters, never individual dependent queries; few-target diagnostic intervals have weak coverage.

F4/F5/SAT and large sparse LA backends were not executed in this round. The tested mechanisms are pair-table PDP policies, Frobenius-orbit base construction and scalar-field Gaussian row kernels. No result here establishes a production-curve crossover or a globally fastest IC method.
