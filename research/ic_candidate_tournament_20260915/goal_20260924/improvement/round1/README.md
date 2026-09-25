# Bounded IC round 1

Round 1: retained. Selected challenger: stop6; retained/promoted winner: incumbent. 3243/3243 native/profile pairs independently verified by the measured round's frozen checker. Fresh Linux transport replay is the evidence-PR merge gate. All results are limited to the registered synthetic toy panel.

[Provenance, interpretation and reproduction commands](EVIDENCE.md).

Primary metric: one supplied target, after reusable preparation through scalar replay; fixture generation is outside both timed algorithms. Cold instruction and native process costs are supplementary promotion gates. Each table uses one cost unit. Values are equal-cell geometric means of three-process per-point medians, with no target amortization.

The IC reference is the qualified `pairinv` source. `rho` is the separately qualified cold-instruction reference; `rho_online` is the separately qualified online-time reference. Ratios to rho are descriptive. The K-instruction floor applies only to this full-rank collector; it is not a generic IC lower bound.

The class column labels engineering experiments and accounting controls. No asymptotic advance is claimed. Variant names are readable aliases; the machine-readable export retains every canonical candidate, workload and run ID. `stop3` uses the legacy adaptive orbit bound, which is not a universal three-column guarantee. Actual admitted bases and columns below are authoritative.

## Confirmation

Verified 864/864 pairs.

### Single-target online time

| Variant | Online ms | / IC reference | Descriptive 95% interval | Online rho / variant | Verified / scheduled | Class |
|---|---|---|---|---|---|---|
| incumbent | 0.0347468 | 1 | reference | 8.41491 | 216/216 | accounting |
| rho | 0.290317 | 8.35522 | [7.33193, 9.376] | 1.00714 | 216/216 | accounting |
| rho_online | 0.292392 | 8.41491 | [7.39792, 9.45168] | 1 | 216/216 | accounting |
| stop6 | 0.036997 | 1.06476 | [1.01009, 1.11961] | 7.90313 | 216/216 | engineering |

### Complete cold instructions

| Variant | Complete cold Ir | S = Ir / sqrt(r) | / IC reference | Descriptive 95% interval | / cold rho | / K floor | Verified / scheduled | Class |
|---|---|---|---|---|---|---|---|---|
| incumbent | 1.78276e+06 | 2790.95 | 1 | reference | 0.717195 | 222845 | 216/216 | accounting |
| rho | 2.48574e+06 | 3891.48 | 1.39432 | [1.24631, 1.55199] | 1 | not applicable | 216/216 | accounting |
| rho_online | 2.48564e+06 | 3891.33 | 1.39427 | [1.24627, 1.55192] | 0.999961 | not applicable | 216/216 | accounting |
| stop6 | 1.70579e+06 | 2670.45 | 0.956822 | [0.921452, 0.992762] | 0.686228 | 284298 | 216/216 | engineering |

### Complete cold native time

| Variant | Complete cold ms | / IC reference | Descriptive 95% interval | Verified / scheduled | Class |
|---|---|---|---|---|---|
| incumbent | 0.97376 | 1 | reference | 216/216 | accounting |
| rho | 0.939833 | 0.965159 | [0.92131, 1.01048] | 216/216 | accounting |
| rho_online | 0.93906 | 0.964365 | [0.919519, 1.0106] | 216/216 | accounting |
| stop6 | 0.939253 | 0.964563 | [0.942415, 0.98683] | 216/216 | engineering |

### Actual base and matrix sizes

| Variant | Curve cell | Usable points B / folded columns K / final rank |
|---|---|---|
| incumbent | n17a1 | 272 / 8 / 8 |
| incumbent | n19a0 | 304 / 8 / 8 |
| incumbent | n23a0 | 368 / 8 / 8 |
| incumbent | n23a1 | 368 / 8 / 8 |
| incumbent | n29a1 | 464 / 8 / 8 |
| incumbent | n31a0 | 496 / 8 / 8 |
| stop6 | n17a1 | 204 / 6 / 6 |
| stop6 | n19a0 | 228 / 6 / 6 |
| stop6 | n23a0 | 276 / 6 / 6 |
| stop6 | n23a1 | 276 / 6 / 6 |
| stop6 | n29a1 | 348 / 6 / 6 |
| stop6 | n31a0 | 372 / 6 / 6 |

## Replay

Verified 864/864 pairs.

### Single-target online time

| Variant | Online ms | / IC reference | Descriptive 95% interval | Online rho / variant | Verified / scheduled | Class |
|---|---|---|---|---|---|---|
| incumbent | 0.0348236 | 1 | reference | 8.37185 | 216/216 | accounting |
| rho | 0.29019 | 8.33315 | [7.38617, 9.34136] | 1.00464 | 216/216 | accounting |
| rho_online | 0.291538 | 8.37185 | [7.47242, 9.31891] | 1 | 216/216 | accounting |
| stop6 | 0.0373129 | 1.07148 | [1.00765, 1.13515] | 7.81333 | 216/216 | engineering |

### Complete cold instructions

| Variant | Complete cold Ir | S = Ir / sqrt(r) | / IC reference | Descriptive 95% interval | / cold rho | / K floor | Verified / scheduled | Class |
|---|---|---|---|---|---|---|---|---|
| incumbent | 1.78276e+06 | 2790.95 | 1 | reference | 0.717195 | 222845 | 216/216 | accounting |
| rho | 2.48574e+06 | 3891.49 | 1.39432 | [1.24632, 1.55199] | 1 | not applicable | 216/216 | accounting |
| rho_online | 2.48565e+06 | 3891.35 | 1.39427 | [1.24628, 1.55194] | 0.999965 | not applicable | 216/216 | accounting |
| stop6 | 1.70578e+06 | 2670.44 | 0.956821 | [0.92145, 0.99276] | 0.686227 | 284297 | 216/216 | engineering |

### Complete cold native time

| Variant | Complete cold ms | / IC reference | Descriptive 95% interval | Verified / scheduled | Class |
|---|---|---|---|---|---|
| incumbent | 0.971718 | 1 | reference | 216/216 | accounting |
| rho | 0.943677 | 0.971143 | [0.927014, 1.01771] | 216/216 | accounting |
| rho_online | 0.943521 | 0.970983 | [0.926287, 1.01869] | 216/216 | accounting |
| stop6 | 0.939902 | 0.967258 | [0.948695, 0.987191] | 216/216 | engineering |

### Actual base and matrix sizes

| Variant | Curve cell | Usable points B / folded columns K / final rank |
|---|---|---|
| incumbent | n17a1 | 272 / 8 / 8 |
| incumbent | n19a0 | 304 / 8 / 8 |
| incumbent | n23a0 | 368 / 8 / 8 |
| incumbent | n23a1 | 368 / 8 / 8 |
| incumbent | n29a1 | 464 / 8 / 8 |
| incumbent | n31a0 | 496 / 8 / 8 |
| stop6 | n17a1 | 204 / 6 / 6 |
| stop6 | n19a0 | 228 / 6 / 6 |
| stop6 | n23a0 | 276 / 6 / 6 |
| stop6 | n23a1 | 276 / 6 / 6 |
| stop6 | n29a1 | 348 / 6 / 6 |
| stop6 | n31a0 | 372 / 6 / 6 |

## Development

Verified 810/810 pairs.

### Single-target online time

| Variant | Online ms | / IC reference | Descriptive 95% interval | Online rho / variant | Verified / scheduled | Class |
|---|---|---|---|---|---|---|
| incumbent | 0.0352563 | 1 | reference | 8.62949 | 45/45 | accounting |
| compatibility | 0.0352661 | 1.00028 | [0.915036, 1.09709] | 8.62712 | 45/45 | accounting |
| rho | 0.302369 | 8.5763 | [7.34981, 10.0703] | 1.0062 | 45/45 | accounting |
| rho_online | 0.304244 | 8.62949 | [7.24042, 10.2083] | 1 | 45/45 | accounting |
| row_bounded | 0.0358234 | 1.01608 | [0.928955, 1.10844] | 8.4929 | 45/45 | engineering |
| row_suffix | 0.0346403 | 0.982528 | [0.892417, 1.08623] | 8.78295 | 45/45 | engineering |
| row_word | 0.0345295 | 0.979383 | [0.896681, 1.06213] | 8.81115 | 45/45 | engineering |
| stop2 | 0.0679139 | 1.92629 | [1.18149, 3.38307] | 4.47986 | 45/45 | engineering |
| stop2_word | 0.0647879 | 1.83762 | [1.11767, 3.26681] | 4.69601 | 45/45 | engineering |
| stop3 | 0.0691609 | 1.96166 | [1.09311, 3.81667] | 4.39908 | 45/45 | engineering |
| stop3_full | 0.0648977 | 1.84074 | [0.998979, 3.64114] | 4.68806 | 45/45 | engineering |
| stop3_word | 0.0668481 | 1.89606 | [1.07219, 3.5761] | 4.55128 | 45/45 | engineering |
| stop3_word_full | 0.063838 | 1.81068 | [0.976064, 3.58731] | 4.76588 | 45/45 | engineering |
| stop4 | 0.0486596 | 1.38017 | [0.897831, 2.25716] | 6.25251 | 45/45 | engineering |
| stop4_word | 0.0499396 | 1.41647 | [0.910688, 2.37107] | 6.09225 | 45/45 | engineering |
| stop4_word_full | 0.0508188 | 1.44141 | [0.896384, 2.48901] | 5.98684 | 45/45 | engineering |
| stop6 | 0.0387835 | 1.10004 | [0.924287, 1.31229] | 7.84469 | 45/45 | engineering |
| stop6_full | 0.0339077 | 0.961748 | [0.854766, 1.13779] | 8.97271 | 45/45 | engineering |

### Complete cold instructions

| Variant | Complete cold Ir | S = Ir / sqrt(r) | / IC reference | Descriptive 95% interval | / cold rho | / K floor | Verified / scheduled | Class |
|---|---|---|---|---|---|---|---|---|
| incumbent | 1.72112e+06 | 2148.81 | 1 | reference | 0.672651 | 215140 | 45/45 | accounting |
| compatibility | 1.72684e+06 | 2155.95 | 1.00332 | [1.00211, 1.00471] | 0.674886 | 215855 | 45/45 | accounting |
| rho | 2.55872e+06 | 3194.54 | 1.48666 | [1.1913, 1.73263] | 1 | not applicable | 45/45 | accounting |
| rho_online | 2.55868e+06 | 3194.5 | 1.48664 | [1.19127, 1.73263] | 0.999986 | not applicable | 45/45 | accounting |
| row_bounded | 1.72106e+06 | 2148.73 | 0.999963 | [0.999128, 1.00084] | 0.672626 | 215133 | 45/45 | engineering |
| row_suffix | 1.72443e+06 | 2152.94 | 1.00192 | [1.00078, 1.00307] | 0.673943 | 215554 | 45/45 | engineering |
| row_word | 1.71926e+06 | 2146.49 | 0.998918 | [0.998177, 0.999686] | 0.671923 | 214908 | 45/45 | engineering |
| stop2 | 2.33929e+06 | 2920.58 | 1.35916 | [0.909484, 2.1533] | 0.914241 | 1.16964e+06 | 45/45 | engineering |
| stop2_word | 2.33893e+06 | 2920.14 | 1.35896 | [0.909258, 2.1532] | 0.914103 | 1.16947e+06 | 45/45 | engineering |
| stop3 | 2.03632e+06 | 2542.33 | 1.18313 | [0.862259, 1.65424] | 0.795836 | 678773 | 45/45 | engineering |
| stop3_full | 2.10341e+06 | 2626.09 | 1.22211 | [0.913466, 1.65705] | 0.822055 | 701136 | 45/45 | engineering |
| stop3_word | 2.03549e+06 | 2541.29 | 1.18265 | [0.861813, 1.65388] | 0.795511 | 678496 | 45/45 | engineering |
| stop3_word_full | 2.10264e+06 | 2625.13 | 1.22166 | [0.912881, 1.65675] | 0.821754 | 700879 | 45/45 | engineering |
| stop4 | 1.80551e+06 | 2254.16 | 1.04903 | [0.832364, 1.36359] | 0.705629 | 451377 | 45/45 | engineering |
| stop4_word | 1.80389e+06 | 2252.14 | 1.04809 | [0.831268, 1.36309] | 0.704997 | 450972 | 45/45 | engineering |
| stop4_word_full | 2.04443e+06 | 2552.46 | 1.18785 | [1.00918, 1.44019] | 0.799005 | 511107 | 45/45 | engineering |
| stop6 | 1.80162e+06 | 2249.31 | 1.04677 | [0.950769, 1.17067] | 0.704111 | 300270 | 45/45 | engineering |
| stop6_full | 2.41789e+06 | 3018.72 | 1.40483 | [1.29982, 1.48749] | 0.944962 | 402982 | 45/45 | engineering |

### Complete cold native time

| Variant | Complete cold ms | / IC reference | Descriptive 95% interval | Verified / scheduled | Class |
|---|---|---|---|---|---|
| incumbent | 0.981844 | 1 | reference | 45/45 | accounting |
| compatibility | 0.984506 | 1.00271 | [0.958741, 1.05815] | 45/45 | accounting |
| rho | 0.946136 | 0.963631 | [0.898189, 1.02345] | 45/45 | accounting |
| rho_online | 0.945346 | 0.962827 | [0.898169, 1.026] | 45/45 | accounting |
| row_bounded | 0.984798 | 1.00301 | [0.94677, 1.0598] | 45/45 | engineering |
| row_suffix | 0.988162 | 1.00643 | [0.960356, 1.07065] | 45/45 | engineering |
| row_word | 0.956282 | 0.973965 | [0.944071, 1.00736] | 45/45 | engineering |
| stop2 | 1.0355 | 1.05465 | [0.886556, 1.30183] | 45/45 | engineering |
| stop2_word | 1.03288 | 1.05198 | [0.899399, 1.28675] | 45/45 | engineering |
| stop3 | 0.97891 | 0.997011 | [0.88293, 1.14783] | 45/45 | engineering |
| stop3_full | 0.988985 | 1.00727 | [0.911653, 1.13801] | 45/45 | engineering |
| stop3_word | 0.956279 | 0.973962 | [0.865902, 1.11229] | 45/45 | engineering |
| stop3_word_full | 0.956479 | 0.974165 | [0.886563, 1.09026] | 45/45 | engineering |
| stop4 | 0.918897 | 0.935889 | [0.859965, 1.02221] | 45/45 | engineering |
| stop4_word | 0.923877 | 0.940961 | [0.85349, 1.02574] | 45/45 | engineering |
| stop4_word_full | 0.948846 | 0.966392 | [0.882228, 1.05218] | 45/45 | engineering |
| stop6 | 0.967292 | 0.985179 | [0.933241, 1.04515] | 45/45 | engineering |
| stop6_full | 1.00556 | 1.02416 | [0.988143, 1.07122] | 45/45 | engineering |

Retained portfolio:

- `stop6_full`: single-target online-time leader.
- `row_word`: complete instruction-cost leader.
- `stop4`: complete native-time leader.
- `row_bounded`: non-dominated implementation family.
- `row_suffix`: non-dominated implementation family.
- `stop6`: predeclared exploration slot.

### Actual base and matrix sizes

| Variant | Curve cell | Usable points B / folded columns K / final rank |
|---|---|---|
| incumbent | n17a1 | 272 / 8 / 8 |
| incumbent | n19a0 | 304 / 8 / 8 |
| incumbent | n23a0 | 368 / 8 / 8 |
| incumbent | n23a1 | 368 / 8 / 8 |
| incumbent | n31a0 | 496 / 8 / 8 |
| compatibility | n17a1 | 272 / 8 / 8 |
| compatibility | n19a0 | 304 / 8 / 8 |
| compatibility | n23a0 | 368 / 8 / 8 |
| compatibility | n23a1 | 368 / 8 / 8 |
| compatibility | n31a0 | 496 / 8 / 8 |
| row_bounded | n17a1 | 272 / 8 / 8 |
| row_bounded | n19a0 | 304 / 8 / 8 |
| row_bounded | n23a0 | 368 / 8 / 8 |
| row_bounded | n23a1 | 368 / 8 / 8 |
| row_bounded | n31a0 | 496 / 8 / 8 |
| row_suffix | n17a1 | 272 / 8 / 8 |
| row_suffix | n19a0 | 304 / 8 / 8 |
| row_suffix | n23a0 | 368 / 8 / 8 |
| row_suffix | n23a1 | 368 / 8 / 8 |
| row_suffix | n31a0 | 496 / 8 / 8 |
| row_word | n17a1 | 272 / 8 / 8 |
| row_word | n19a0 | 304 / 8 / 8 |
| row_word | n23a0 | 368 / 8 / 8 |
| row_word | n23a1 | 368 / 8 / 8 |
| row_word | n31a0 | 496 / 8 / 8 |
| stop2 | n17a1 | 68 / 2 / 2 |
| stop2 | n19a0 | 76 / 2 / 2 |
| stop2 | n23a0 | 92 / 2 / 2 |
| stop2 | n23a1 | 92 / 2 / 2 |
| stop2 | n31a0 | 124 / 2 / 2 |
| stop2_word | n17a1 | 68 / 2 / 2 |
| stop2_word | n19a0 | 76 / 2 / 2 |
| stop2_word | n23a0 | 92 / 2 / 2 |
| stop2_word | n23a1 | 92 / 2 / 2 |
| stop2_word | n31a0 | 124 / 2 / 2 |
| stop3 | n17a1 | 102 / 3 / 3 |
| stop3 | n19a0 | 114 / 3 / 3 |
| stop3 | n23a0 | 138 / 3 / 3 |
| stop3 | n23a1 | 138 / 3 / 3 |
| stop3 | n31a0 | 186 / 3 / 3 |
| stop3_full | n17a1 | 102 / 3 / 3 |
| stop3_full | n19a0 | 114 / 3 / 3 |
| stop3_full | n23a0 | 138 / 3 / 3 |
| stop3_full | n23a1 | 138 / 3 / 3 |
| stop3_full | n31a0 | 186 / 3 / 3 |
| stop3_word | n17a1 | 102 / 3 / 3 |
| stop3_word | n19a0 | 114 / 3 / 3 |
| stop3_word | n23a0 | 138 / 3 / 3 |
| stop3_word | n23a1 | 138 / 3 / 3 |
| stop3_word | n31a0 | 186 / 3 / 3 |
| stop3_word_full | n17a1 | 102 / 3 / 3 |
| stop3_word_full | n19a0 | 114 / 3 / 3 |
| stop3_word_full | n23a0 | 138 / 3 / 3 |
| stop3_word_full | n23a1 | 138 / 3 / 3 |
| stop3_word_full | n31a0 | 186 / 3 / 3 |
| stop4 | n17a1 | 136 / 4 / 4 |
| stop4 | n19a0 | 152 / 4 / 4 |
| stop4 | n23a0 | 184 / 4 / 4 |
| stop4 | n23a1 | 184 / 4 / 4 |
| stop4 | n31a0 | 248 / 4 / 4 |
| stop4_word | n17a1 | 136 / 4 / 4 |
| stop4_word | n19a0 | 152 / 4 / 4 |
| stop4_word | n23a0 | 184 / 4 / 4 |
| stop4_word | n23a1 | 184 / 4 / 4 |
| stop4_word | n31a0 | 248 / 4 / 4 |
| stop4_word_full | n17a1 | 136 / 4 / 4 |
| stop4_word_full | n19a0 | 152 / 4 / 4 |
| stop4_word_full | n23a0 | 184 / 4 / 4 |
| stop4_word_full | n23a1 | 184 / 4 / 4 |
| stop4_word_full | n31a0 | 248 / 4 / 4 |
| stop6 | n17a1 | 204 / 6 / 6 |
| stop6 | n19a0 | 228 / 6 / 6 |
| stop6 | n23a0 | 276 / 6 / 6 |
| stop6 | n23a1 | 276 / 6 / 6 |
| stop6 | n31a0 | 372 / 6 / 6 |
| stop6_full | n17a1 | 204 / 6 / 6 |
| stop6_full | n19a0 | 228 / 6 / 6 |
| stop6_full | n23a0 | 276 / 6 / 6 |
| stop6_full | n23a1 | 276 / 6 / 6 |
| stop6_full | n31a0 | 372 / 6 / 6 |

## Smoke

Verified 270/270 pairs.

### Single-target online time

| Variant | Online ms | / IC reference | Descriptive 95% interval | Online rho / variant | Verified / scheduled | Class |
|---|---|---|---|---|---|---|
| incumbent | 0.0414595 | 1 | reference | 8.24864 | 15/15 | accounting |
| compatibility | 0.0466808 | 1.12594 | [0.998912, 1.42512] | 7.32603 | 15/15 | accounting |
| rho | 0.331739 | 8.00152 | [6.91821, 9.25446] | 1.03088 | 15/15 | accounting |
| rho_online | 0.341985 | 8.24864 | [6.93628, 9.8093] | 1 | 15/15 | accounting |
| row_bounded | 0.0418418 | 1.00922 | [0.988594, 1.03028] | 8.17327 | 15/15 | engineering |
| row_suffix | 0.0442462 | 1.06722 | [1.00365, 1.14663] | 7.72913 | 15/15 | engineering |
| row_word | 0.0434484 | 1.04797 | [0.959907, 1.1652] | 7.87105 | 15/15 | engineering |
| stop2 | 0.0384477 | 0.927355 | [0.682692, 1.29242] | 8.8948 | 15/15 | engineering |
| stop2_word | 0.0385728 | 0.930373 | [0.681697, 1.27986] | 8.86595 | 15/15 | engineering |
| stop3 | 0.0389066 | 0.938423 | [0.781481, 1.07914] | 8.7899 | 15/15 | engineering |
| stop3_full | 0.0388243 | 0.936439 | [0.784051, 1.08684] | 8.80852 | 15/15 | engineering |
| stop3_word | 0.0400013 | 0.964829 | [0.799571, 1.12555] | 8.54933 | 15/15 | engineering |
| stop3_word_full | 0.0382545 | 0.922696 | [0.789234, 1.05734] | 8.93972 | 15/15 | engineering |
| stop4 | 0.0393433 | 0.948958 | [0.793684, 1.12729] | 8.69232 | 15/15 | engineering |
| stop4_word | 0.0416154 | 1.00376 | [0.80422, 1.30243] | 8.21774 | 15/15 | engineering |
| stop4_word_full | 0.0400582 | 0.966199 | [0.771701, 1.17757] | 8.53721 | 15/15 | engineering |
| stop6 | 0.0414627 | 1.00008 | [0.705577, 1.46555] | 8.248 | 15/15 | engineering |
| stop6_full | 0.0400991 | 0.967185 | [0.6931, 1.36308] | 8.5285 | 15/15 | engineering |

### Complete cold instructions

| Variant | Complete cold Ir | S = Ir / sqrt(r) | / IC reference | Descriptive 95% interval | / cold rho | / K floor | Verified / scheduled | Class |
|---|---|---|---|---|---|---|---|---|
| incumbent | 1.94699e+06 | 2430.81 | 1 | reference | 0.689415 | 243374 | 15/15 | accounting |
| compatibility | 1.95332e+06 | 2438.71 | 1.00325 | [1.00209, 1.00445] | 0.691656 | 244165 | 15/15 | accounting |
| rho | 2.82412e+06 | 3525.9 | 1.45051 | [1.2103, 1.70404] | 1 | not applicable | 15/15 | accounting |
| rho_online | 2.82409e+06 | 3525.86 | 1.45049 | [1.21028, 1.70403] | 0.99999 | not applicable | 15/15 | accounting |
| row_bounded | 1.9457e+06 | 2429.19 | 0.999336 | [0.998698, 0.999974] | 0.688956 | 243212 | 15/15 | engineering |
| row_suffix | 1.95045e+06 | 2435.13 | 1.00178 | [1.00064, 1.00299] | 0.69064 | 243806 | 15/15 | engineering |
| row_word | 1.9436e+06 | 2426.57 | 0.998257 | [0.997788, 0.998726] | 0.688213 | 242950 | 15/15 | engineering |
| stop2 | 1.68904e+06 | 2108.76 | 0.867514 | [0.764624, 0.984249] | 0.598077 | 844521 | 15/15 | engineering |
| stop2_word | 1.68874e+06 | 2108.38 | 0.867357 | [0.76442, 0.984155] | 0.597968 | 844368 | 15/15 | engineering |
| stop3 | 1.60918e+06 | 2009.05 | 0.826497 | [0.701536, 1.0244] | 0.569799 | 536394 | 15/15 | engineering |
| stop3_full | 1.71499e+06 | 2141.15 | 0.880841 | [0.74517, 1.07333] | 0.607265 | 571663 | 15/15 | engineering |
| stop3_word | 1.60846e+06 | 2008.15 | 0.826127 | [0.701291, 1.02407] | 0.569544 | 536154 | 15/15 | engineering |
| stop3_word_full | 1.71429e+06 | 2140.28 | 0.880483 | [0.744926, 1.07293] | 0.607018 | 571431 | 15/15 | engineering |
| stop4 | 1.60745e+06 | 2006.89 | 0.825608 | [0.675293, 1.04043] | 0.569186 | 401863 | 15/15 | engineering |
| stop4_word | 1.6059e+06 | 2004.96 | 0.824813 | [0.674708, 1.03982] | 0.568638 | 401476 | 15/15 | engineering |
| stop4_word_full | 1.84741e+06 | 2306.48 | 0.948854 | [0.75375, 1.14322] | 0.654154 | 461852 | 15/15 | engineering |
| stop6 | 1.58274e+06 | 1976.04 | 0.812916 | [0.745035, 0.882935] | 0.560436 | 263790 | 15/15 | engineering |
| stop6_full | 2.28429e+06 | 2851.92 | 1.17324 | [0.939703, 1.42684] | 0.808848 | 380715 | 15/15 | engineering |

### Complete cold native time

| Variant | Complete cold ms | / IC reference | Descriptive 95% interval | Verified / scheduled | Class |
|---|---|---|---|---|---|
| incumbent | 1.00367 | 1 | reference | 15/15 | accounting |
| compatibility | 1.01525 | 1.01154 | [0.96277, 1.06977] | 15/15 | accounting |
| rho | 0.987279 | 0.983669 | [0.934957, 1.03761] | 15/15 | accounting |
| rho_online | 0.996115 | 0.992472 | [0.95043, 1.04289] | 15/15 | accounting |
| row_bounded | 1.01107 | 1.00738 | [0.961045, 1.0607] | 15/15 | engineering |
| row_suffix | 1.03722 | 1.03343 | [0.982347, 1.09666] | 15/15 | engineering |
| row_word | 1.02724 | 1.02348 | [0.980733, 1.0566] | 15/15 | engineering |
| stop2 | 0.868114 | 0.864939 | [0.826586, 0.901316] | 15/15 | engineering |
| stop2_word | 0.927407 | 0.924015 | [0.870228, 0.981128] | 15/15 | engineering |
| stop3 | 0.882391 | 0.879163 | [0.848425, 0.916289] | 15/15 | engineering |
| stop3_full | 0.899878 | 0.896587 | [0.856064, 0.947817] | 15/15 | engineering |
| stop3_word | 0.882449 | 0.879222 | [0.806756, 0.958197] | 15/15 | engineering |
| stop3_word_full | 0.877445 | 0.874236 | [0.821711, 0.943928] | 15/15 | engineering |
| stop4 | 0.895041 | 0.891768 | [0.859841, 0.935816] | 15/15 | engineering |
| stop4_word | 0.909145 | 0.90582 | [0.863603, 0.950101] | 15/15 | engineering |
| stop4_word_full | 0.922189 | 0.918817 | [0.867159, 0.988972] | 15/15 | engineering |
| stop6 | 0.904897 | 0.901587 | [0.876533, 0.926856] | 15/15 | engineering |
| stop6_full | 0.985016 | 0.981413 | [0.931481, 1.04672] | 15/15 | engineering |

### Actual base and matrix sizes

| Variant | Curve cell | Usable points B / folded columns K / final rank |
|---|---|---|
| incumbent | n17a1 | 272 / 8 / 8 |
| incumbent | n19a0 | 304 / 8 / 8 |
| incumbent | n23a0 | 368 / 8 / 8 |
| incumbent | n23a1 | 368 / 8 / 8 |
| incumbent | n31a0 | 496 / 8 / 8 |
| compatibility | n17a1 | 272 / 8 / 8 |
| compatibility | n19a0 | 304 / 8 / 8 |
| compatibility | n23a0 | 368 / 8 / 8 |
| compatibility | n23a1 | 368 / 8 / 8 |
| compatibility | n31a0 | 496 / 8 / 8 |
| row_bounded | n17a1 | 272 / 8 / 8 |
| row_bounded | n19a0 | 304 / 8 / 8 |
| row_bounded | n23a0 | 368 / 8 / 8 |
| row_bounded | n23a1 | 368 / 8 / 8 |
| row_bounded | n31a0 | 496 / 8 / 8 |
| row_suffix | n17a1 | 272 / 8 / 8 |
| row_suffix | n19a0 | 304 / 8 / 8 |
| row_suffix | n23a0 | 368 / 8 / 8 |
| row_suffix | n23a1 | 368 / 8 / 8 |
| row_suffix | n31a0 | 496 / 8 / 8 |
| row_word | n17a1 | 272 / 8 / 8 |
| row_word | n19a0 | 304 / 8 / 8 |
| row_word | n23a0 | 368 / 8 / 8 |
| row_word | n23a1 | 368 / 8 / 8 |
| row_word | n31a0 | 496 / 8 / 8 |
| stop2 | n17a1 | 68 / 2 / 2 |
| stop2 | n19a0 | 76 / 2 / 2 |
| stop2 | n23a0 | 92 / 2 / 2 |
| stop2 | n23a1 | 92 / 2 / 2 |
| stop2 | n31a0 | 124 / 2 / 2 |
| stop2_word | n17a1 | 68 / 2 / 2 |
| stop2_word | n19a0 | 76 / 2 / 2 |
| stop2_word | n23a0 | 92 / 2 / 2 |
| stop2_word | n23a1 | 92 / 2 / 2 |
| stop2_word | n31a0 | 124 / 2 / 2 |
| stop3 | n17a1 | 102 / 3 / 3 |
| stop3 | n19a0 | 114 / 3 / 3 |
| stop3 | n23a0 | 138 / 3 / 3 |
| stop3 | n23a1 | 138 / 3 / 3 |
| stop3 | n31a0 | 186 / 3 / 3 |
| stop3_full | n17a1 | 102 / 3 / 3 |
| stop3_full | n19a0 | 114 / 3 / 3 |
| stop3_full | n23a0 | 138 / 3 / 3 |
| stop3_full | n23a1 | 138 / 3 / 3 |
| stop3_full | n31a0 | 186 / 3 / 3 |
| stop3_word | n17a1 | 102 / 3 / 3 |
| stop3_word | n19a0 | 114 / 3 / 3 |
| stop3_word | n23a0 | 138 / 3 / 3 |
| stop3_word | n23a1 | 138 / 3 / 3 |
| stop3_word | n31a0 | 186 / 3 / 3 |
| stop3_word_full | n17a1 | 102 / 3 / 3 |
| stop3_word_full | n19a0 | 114 / 3 / 3 |
| stop3_word_full | n23a0 | 138 / 3 / 3 |
| stop3_word_full | n23a1 | 138 / 3 / 3 |
| stop3_word_full | n31a0 | 186 / 3 / 3 |
| stop4 | n17a1 | 136 / 4 / 4 |
| stop4 | n19a0 | 152 / 4 / 4 |
| stop4 | n23a0 | 184 / 4 / 4 |
| stop4 | n23a1 | 184 / 4 / 4 |
| stop4 | n31a0 | 248 / 4 / 4 |
| stop4_word | n17a1 | 136 / 4 / 4 |
| stop4_word | n19a0 | 152 / 4 / 4 |
| stop4_word | n23a0 | 184 / 4 / 4 |
| stop4_word | n23a1 | 184 / 4 / 4 |
| stop4_word | n31a0 | 248 / 4 / 4 |
| stop4_word_full | n17a1 | 136 / 4 / 4 |
| stop4_word_full | n19a0 | 152 / 4 / 4 |
| stop4_word_full | n23a0 | 184 / 4 / 4 |
| stop4_word_full | n23a1 | 184 / 4 / 4 |
| stop4_word_full | n31a0 | 248 / 4 / 4 |
| stop6 | n17a1 | 204 / 6 / 6 |
| stop6 | n19a0 | 228 / 6 / 6 |
| stop6 | n23a0 | 276 / 6 / 6 |
| stop6 | n23a1 | 276 / 6 / 6 |
| stop6 | n31a0 | 372 / 6 / 6 |
| stop6_full | n17a1 | 204 / 6 / 6 |
| stop6_full | n19a0 | 228 / 6 / 6 |
| stop6_full | n23a0 | 276 / 6 / 6 |
| stop6_full | n23a1 | 276 / 6 / 6 |
| stop6_full | n31a0 | 372 / 6 / 6 |

## Selection

Verified 405/405 pairs.

### Single-target online time

| Variant | Online ms | / IC reference | Descriptive 95% interval | Online rho / variant | Verified / scheduled | Class |
|---|---|---|---|---|---|---|
| incumbent | 0.036151 | 1 | reference | 8.45799 | 45/45 | accounting |
| rho | 0.305436 | 8.44889 | [7.15204, 9.79975] | 1.00108 | 45/45 | accounting |
| rho_online | 0.305765 | 8.45799 | [7.14888, 9.82322] | 1 | 45/45 | accounting |
| row_bounded | 0.0380344 | 1.0521 | [0.994404, 1.11908] | 8.03916 | 45/45 | engineering |
| row_suffix | 0.0373456 | 1.03304 | [0.967833, 1.09955] | 8.18745 | 45/45 | engineering |
| row_word | 0.0377286 | 1.04364 | [0.971463, 1.11319] | 8.10432 | 45/45 | engineering |
| stop4 | 0.0528517 | 1.46197 | [1.08023, 2.11107] | 5.78534 | 45/45 | engineering |
| stop6 | 0.0419274 | 1.15978 | [0.946978, 1.50778] | 7.29272 | 45/45 | engineering |
| stop6_full | 0.0368712 | 1.01992 | [0.896244, 1.17588] | 8.29277 | 45/45 | engineering |

### Complete cold instructions

| Variant | Complete cold Ir | S = Ir / sqrt(r) | / IC reference | Descriptive 95% interval | / cold rho | / K floor | Verified / scheduled | Class |
|---|---|---|---|---|---|---|---|---|
| incumbent | 1.81215e+06 | 2262.46 | 1 | reference | 0.688648 | 226519 | 45/45 | accounting |
| rho | 2.63146e+06 | 3285.36 | 1.45212 | [1.26781, 1.6097] | 1 | not applicable | 45/45 | accounting |
| rho_online | 2.63137e+06 | 3285.25 | 1.45207 | [1.2678, 1.60964] | 0.999966 | not applicable | 45/45 | accounting |
| row_bounded | 1.81158e+06 | 2261.74 | 0.999685 | [0.998993, 1.00034] | 0.688431 | 226447 | 45/45 | engineering |
| row_suffix | 1.81545e+06 | 2266.58 | 1.00182 | [1.00072, 1.00305] | 0.689902 | 226931 | 45/45 | engineering |
| row_word | 1.80959e+06 | 2259.26 | 0.998588 | [0.997831, 0.999242] | 0.687676 | 226199 | 45/45 | engineering |
| stop4 | 1.86225e+06 | 2325.01 | 1.02765 | [0.867977, 1.25418] | 0.707687 | 465562 | 45/45 | engineering |
| stop6 | 1.75971e+06 | 2196.99 | 0.971064 | [0.912361, 1.03234] | 0.668722 | 293286 | 45/45 | engineering |
| stop6_full | 2.3675e+06 | 2955.81 | 1.30646 | [1.15063, 1.46212] | 0.899692 | 394584 | 45/45 | engineering |

### Complete cold native time

| Variant | Complete cold ms | / IC reference | Descriptive 95% interval | Verified / scheduled | Class |
|---|---|---|---|---|---|
| incumbent | 0.941602 | 1 | reference | 45/45 | accounting |
| rho | 0.944123 | 1.00268 | [0.945116, 1.05741] | 45/45 | accounting |
| rho_online | 0.942623 | 1.00108 | [0.94421, 1.05222] | 45/45 | accounting |
| row_bounded | 0.962518 | 1.02221 | [0.966415, 1.08933] | 45/45 | engineering |
| row_suffix | 0.948079 | 1.00688 | [0.960338, 1.0527] | 45/45 | engineering |
| row_word | 0.964427 | 1.02424 | [0.98869, 1.06165] | 45/45 | engineering |
| stop4 | 0.930759 | 0.988484 | [0.902456, 1.092] | 45/45 | engineering |
| stop6 | 0.917517 | 0.974421 | [0.928304, 1.02549] | 45/45 | engineering |
| stop6_full | 0.983063 | 1.04403 | [0.993541, 1.09937] | 45/45 | engineering |

### Actual base and matrix sizes

| Variant | Curve cell | Usable points B / folded columns K / final rank |
|---|---|---|
| incumbent | n17a1 | 272 / 8 / 8 |
| incumbent | n19a0 | 304 / 8 / 8 |
| incumbent | n23a0 | 368 / 8 / 8 |
| incumbent | n23a1 | 368 / 8 / 8 |
| incumbent | n31a0 | 496 / 8 / 8 |
| row_bounded | n17a1 | 272 / 8 / 8 |
| row_bounded | n19a0 | 304 / 8 / 8 |
| row_bounded | n23a0 | 368 / 8 / 8 |
| row_bounded | n23a1 | 368 / 8 / 8 |
| row_bounded | n31a0 | 496 / 8 / 8 |
| row_suffix | n17a1 | 272 / 8 / 8 |
| row_suffix | n19a0 | 304 / 8 / 8 |
| row_suffix | n23a0 | 368 / 8 / 8 |
| row_suffix | n23a1 | 368 / 8 / 8 |
| row_suffix | n31a0 | 496 / 8 / 8 |
| row_word | n17a1 | 272 / 8 / 8 |
| row_word | n19a0 | 304 / 8 / 8 |
| row_word | n23a0 | 368 / 8 / 8 |
| row_word | n23a1 | 368 / 8 / 8 |
| row_word | n31a0 | 496 / 8 / 8 |
| stop4 | n17a1 | 136 / 4 / 4 |
| stop4 | n19a0 | 152 / 4 / 4 |
| stop4 | n23a0 | 184 / 4 / 4 |
| stop4 | n23a1 | 184 / 4 / 4 |
| stop4 | n31a0 | 248 / 4 / 4 |
| stop6 | n17a1 | 204 / 6 / 6 |
| stop6 | n19a0 | 228 / 6 / 6 |
| stop6 | n23a0 | 276 / 6 / 6 |
| stop6 | n23a1 | 276 / 6 / 6 |
| stop6 | n31a0 | 372 / 6 / 6 |
| stop6_full | n17a1 | 204 / 6 / 6 |
| stop6_full | n19a0 | 228 / 6 / 6 |
| stop6_full | n23a0 | 276 / 6 / 6 |
| stop6_full | n23a1 | 276 / 6 / 6 |
| stop6_full | n31a0 | 372 / 6 / 6 |

## Aa

Verified 30/30 pairs.

### Single-target online time

| Variant | Online ms | / IC reference | Descriptive 95% interval | Online rho / variant | Verified / scheduled | Class |
|---|---|---|---|---|---|---|
| incumbent | 0.0328128 | 1 | reference | unknown | 15/15 | accounting |
| aa_control | 0.033444 | 1.01924 | [0.98013, 1.06977] | unknown | 15/15 | accounting |

### Complete cold instructions

| Variant | Complete cold Ir | S = Ir / sqrt(r) | / IC reference | Descriptive 95% interval | / cold rho | / K floor | Verified / scheduled | Class |
|---|---|---|---|---|---|---|---|---|
| incumbent | 1.79772e+06 | 2244.44 | 1 | reference | unknown | 224715 | 15/15 | accounting |
| aa_control | 1.79772e+06 | 2244.44 | 1 | [1, 1] | unknown | 224715 | 15/15 | accounting |

### Complete cold native time

| Variant | Complete cold ms | / IC reference | Descriptive 95% interval | Verified / scheduled | Class |
|---|---|---|---|---|---|
| incumbent | 0.953169 | 1 | reference | 15/15 | accounting |
| aa_control | 0.995065 | 1.04395 | [1.00605, 1.08036] | 15/15 | accounting |

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
| online_ns | 1.06476 | 1.14271 | 1.28251 |
| instructions | 0.956822 | 1.00812 | 1.09327 |
| cold_ns | 0.964563 | 0.995827 | 0.991782 |

### Replay promotion evidence

| Metric | Candidate / reference | Familywise upper | Largest cell ratio |
|---|---|---|---|
| online_ns | 1.07148 | 1.15926 | 1.24018 |
| instructions | 0.956821 | 1.00812 | 1.09327 |
| cold_ns | 0.967258 | 0.995708 | 1.01595 |

Decision reasons:

- No challenger passed every confirmation and replay threshold.

Complete run keys, admitted B/column counts, exclusive instruction ledgers, online phase clocks, peak process RSS and collection statistics are retained in `RUNS.csv` and `RESULTS.json`; raw certificates, sources, binaries and profiles are in the archive. Base-only memory was not instrumented and remains unknown. Yield intervals resample whole target/walk clusters, never individual dependent queries; few-target diagnostic intervals have weak coverage.

F4/F5/SAT and large sparse LA backends were not executed in this round. The tested mechanisms are pair-table PDP policies, Frobenius-orbit base construction and scalar-field Gaussian row kernels. No result here establishes a production-curve crossover or a globally fastest IC method.
