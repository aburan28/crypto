| part | bits | kernel | timer | A median | A2 median | B median | A/B (ratio of medians) | A/B per round [min, max] | A/A2 per round [min, max] | verdict | pinned | contended (A, A2, B) |
|---|---:|---|---|---:|---:|---:|---:|---|---|---|---|---|
| bsgs1 | 40 | chain | total_s | 0.1090 | 0.1077 | 0.0858 | 1.270 | [0.959, 1.492] | [0.842, 1.481] | no measurable change | identical | (0, 0, 0) |
| bsgs1 | 40 | chain | build_s | 0.0693 | 0.0671 | 0.0514 | 1.348 | [0.972, 1.554] | [0.771, 1.469] | no measurable change | identical | (0, 0, 0) |
| bsgs1 | 40 | chain | search_s | 0.0384 | 0.0386 | 0.0332 | 1.157 | [0.939, 1.431] | [0.857, 1.502] | no measurable change | identical | (0, 0, 0) |
| bsgs1 | 40 | paired | total_s | 0.0565 | 0.0638 | 0.0417 | 1.354 | [1.026, 1.540] | [0.825, 1.067] | no measurable change | identical | (1, 2, 2) |
| bsgs1 | 40 | paired | build_s | 0.0346 | 0.0381 | 0.0249 | 1.390 | [1.021, 1.559] | [0.852, 1.075] | no measurable change | identical | (1, 2, 2) |
| bsgs1 | 40 | paired | search_s | 0.0226 | 0.0251 | 0.0169 | 1.340 | [1.034, 1.512] | [0.788, 1.057] | no measurable change | identical | (1, 2, 2) |
| bsgs1 | 44 | chain | total_s | 0.3482 | 0.3253 | 0.2414 | 1.442 | [1.287, 1.617] | [1.024, 1.222] | improves | identical | (0, 1, 1) |
| bsgs1 | 44 | chain | build_s | 0.3339 | 0.3123 | 0.2315 | 1.442 | [1.287, 1.616] | [1.025, 1.229] | improves | identical | (0, 1, 1) |
| bsgs1 | 44 | chain | search_s | 0.0142 | 0.0128 | 0.0103 | 1.372 | [1.196, 1.632] | [1.008, 1.294] | no measurable change | identical | (0, 1, 1) |
| bsgs1 | 44 | paired | total_s | 0.2265 | 0.2269 | 0.1142 | 1.984 | [1.320, 2.652] | [0.902, 1.108] | improves | identical | (0, 0, 1) |
| bsgs1 | 44 | paired | build_s | 0.2176 | 0.2153 | 0.1086 | 2.003 | [1.318, 2.678] | [0.913, 1.112] | improves | identical | (0, 0, 1) |
| bsgs1 | 44 | paired | search_s | 0.0098 | 0.0106 | 0.0052 | 1.885 | [1.388, 2.130] | [0.703, 1.107] | improves | identical | (0, 0, 1) |
| bsgs1 | 48 | chain | total_s | 1.4351 | 1.4714 | 1.1003 | 1.304 | [1.243, 1.401] | [0.932, 0.996] | improves | identical | (2, 1, 1) |
| bsgs1 | 48 | chain | build_s | 1.4303 | 1.4665 | 1.0971 | 1.304 | [1.242, 1.401] | [0.932, 0.996] | improves | identical | (2, 1, 1) |
| bsgs1 | 48 | chain | search_s | 0.0048 | 0.0049 | 0.0034 | 1.391 | [1.353, 1.500] | [0.885, 1.021] | improves | identical | (2, 1, 1) |
| bsgs1 | 48 | paired | total_s | 1.0429 | 1.0561 | 0.5543 | 1.881 | [1.739, 1.957] | [0.956, 1.028] | improves | identical | (0, 2, 0) |
| bsgs1 | 48 | paired | build_s | 1.0384 | 1.0523 | 0.5522 | 1.880 | [1.737, 1.957] | [0.957, 1.028] | improves | identical | (0, 2, 0) |
| bsgs1 | 48 | paired | search_s | 0.0040 | 0.0039 | 0.0020 | 2.000 | [1.850, 3.294] | [0.537, 1.436] | improves | identical | (0, 2, 0) |
| bsgs1 | 52 | chain | total_s | 6.5163 | 6.4385 | 4.9373 | 1.320 | [1.241, 1.402] | [0.999, 1.015] | improves | identical | (0, 0, 0) |
| bsgs1 | 52 | chain | build_s | 6.5145 | 6.4366 | 4.9359 | 1.320 | [1.240, 1.402] | [0.999, 1.015] | improves | identical | (0, 0, 0) |
| bsgs1 | 52 | chain | search_s | 0.0019 | 0.0019 | 0.0014 | 1.357 | [1.286, 1.769] | [0.613, 1.353] | no measurable change | identical | (0, 0, 0) |
| bsgs1 | 52 | paired | total_s | 4.6114 | 4.5582 | 2.6462 | 1.743 | [1.667, 1.847] | [0.947, 1.038] | improves | identical | (0, 0, 0) |
| bsgs1 | 52 | paired | build_s | 4.6099 | 4.5567 | 2.6453 | 1.743 | [1.667, 1.847] | [0.947, 1.038] | improves | identical | (0, 0, 0) |
| bsgs1 | 52 | paired | search_s | 0.0015 | 0.0015 | 0.0009 | 1.667 | [1.500, 1.889] | [0.875, 1.133] | improves | identical | (0, 0, 0) |
| bsgs3 | 48 | chain | total_s | 0.5810 | 0.5514 | 0.4350 | 1.336 | [1.232, 1.440] | [0.942, 1.165] | improves | identical | (1, 2, 0) |
| bsgs3 | 48 | chain | build_s | 0.5739 | 0.5461 | 0.4303 | 1.334 | [1.234, 1.447] | [0.940, 1.167] | improves | identical | (1, 2, 0) |
| bsgs3 | 48 | chain | search_s | 0.0060 | 0.0055 | 0.0038 | 1.592 | [0.912, 2.286] | [0.627, 1.346] | no measurable change | identical | (1, 2, 0) |
| bsgs3 | 48 | paired | total_s | 0.4152 | 0.4604 | 0.2767 | 1.501 | [1.386, 1.853] | [0.855, 1.122] | improves | identical | (0, 1, 2) |
| bsgs3 | 48 | paired | build_s | 0.4107 | 0.4558 | 0.2744 | 1.497 | [1.385, 1.845] | [0.853, 1.120] | improves | identical | (0, 1, 2) |
| bsgs3 | 48 | paired | search_s | 0.0044 | 0.0043 | 0.0025 | 1.760 | [1.538, 2.857] | [0.848, 1.304] | improves | identical | (0, 1, 2) |
| rho1 | 32 | - | solve_s | 0.0312 | 0.0309 | 0.0053 | 5.836 | [3.577, 7.177] | [0.708, 1.317] | improves | identical | (0, 0, 0) |
| rho1 | 32 | - | setup_s | 0.0295 | 0.0299 | 0.0043 | 6.857 | [2.490, 8.890] | [0.716, 1.334] | improves | identical | (0, 0, 0) |
| rho1 | 36 | - | solve_s | 0.0532 | 0.0502 | 0.0142 | 3.745 | [2.505, 3.912] | [0.698, 1.139] | improves | identical | (0, 0, 1) |
| rho1 | 36 | - | setup_s | 0.0417 | 0.0407 | 0.0052 | 7.961 | [5.937, 10.321] | [0.892, 1.334] | improves | identical | (0, 0, 1) |
| rho1 | 40 | - | solve_s | 0.1071 | 0.1044 | 0.0588 | 1.822 | [1.448, 2.075] | [0.675, 1.389] | improves | identical | (0, 0, 1) |
| rho1 | 40 | - | setup_s | 0.0521 | 0.0527 | 0.0057 | 9.105 | [7.272, 11.092] | [0.771, 1.489] | improves | identical | (0, 0, 1) |
| rho1 | 44 | - | solve_s | 0.4518 | 0.4401 | 0.3894 | 1.160 | [0.979, 1.457] | [0.954, 1.206] | no measurable change | identical | (0, 0, 0) |
| rho1 | 44 | - | setup_s | 0.0700 | 0.0690 | 0.0064 | 10.959 | [9.597, 13.659] | [0.964, 1.500] | improves | identical | (0, 0, 0) |
| rho1 | 48 | - | solve_s | 1.1729 | 1.1583 | 1.0332 | 1.135 | [1.070, 1.270] | [0.942, 1.045] | improves | identical | (0, 0, 0) |
| rho1 | 48 | - | setup_s | 0.0908 | 0.0888 | 0.0070 | 13.023 | [9.972, 18.749] | [0.811, 1.168] | improves | identical | (0, 0, 0) |
| rho1 | 52 | - | solve_s | 1.7896 | 1.6843 | 1.6480 | 1.086 | [0.902, 1.151] | [0.830, 1.136] | no measurable change | identical | (0, 0, 0) |
| rho1 | 52 | - | setup_s | 0.1192 | 0.1092 | 0.0074 | 16.094 | [13.658, 25.896] | [0.897, 1.373] | improves | identical | (0, 0, 0) |
| rho3 | 32 | - | solve_s | 0.0283 | 0.0271 | 0.0109 | 2.595 | [1.978, 4.169] | [0.708, 1.359] | improves | identical | (1, 0, 0) |
| rho3 | 32 | - | setup_s | 0.0294 | 0.0300 | 0.0043 | 6.774 | [5.313, 9.146] | [0.715, 1.098] | improves | identical | (1, 0, 0) |
| rho3 | 36 | - | solve_s | 0.0714 | 0.0731 | 0.0364 | 1.961 | [1.627, 2.097] | [0.799, 1.271] | improves | identical | (1, 1, 0) |
| rho3 | 36 | - | setup_s | 0.0830 | 0.0809 | 0.0103 | 8.029 | [6.499, 8.524] | [0.879, 1.121] | improves | identical | (1, 1, 0) |
| rho3 | 40 | - | solve_s | 0.0728 | 0.0732 | 0.0304 | 2.393 | [1.997, 3.793] | [0.661, 1.406] | improves | identical | (0, 0, 0) |
| rho3 | 40 | - | setup_s | 0.0968 | 0.1026 | 0.0114 | 8.500 | [6.877, 10.815] | [0.771, 1.135] | improves | identical | (0, 0, 0) |
| rho3 | 44 | - | solve_s | 0.1612 | 0.1552 | 0.0928 | 1.737 | [1.496, 2.166] | [0.981, 1.246] | improves | identical | (0, 0, 0) |
| rho3 | 44 | - | setup_s | 0.1393 | 0.1385 | 0.0127 | 10.996 | [9.052, 12.045] | [0.930, 1.313] | improves | identical | (0, 0, 0) |
| rho3 | 48 | - | solve_s | 0.2142 | 0.2119 | 0.1233 | 1.737 | [1.546, 2.195] | [0.947, 1.218] | improves | identical | (1, 0, 0) |
| rho3 | 48 | - | setup_s | 0.1872 | 0.1770 | 0.0139 | 13.506 | [11.184, 16.549] | [0.827, 1.317] | improves | identical | (1, 0, 0) |
| rho3 | 52 | - | solve_s | 1.2046 | 1.2155 | 1.0852 | 1.110 | [1.014, 1.144] | [0.909, 1.035] | no measurable change | identical | (0, 0, 0) |
| rho3 | 52 | - | setup_s | 0.2186 | 0.2105 | 0.0148 | 14.783 | [11.676, 18.252] | [0.840, 1.113] | improves | identical | (0, 0, 0) |

| curve role | bits | PARI elllog median s (n) | native A median s | native B median s | B / PARI | A / PARI | PARI answers verified |
|---|---:|---:|---:|---:|---:|---:|---|
| bsgs | 40 | 1.006 (3) | 0.0565 | 0.0417 | 0.0415 | 0.0561 | true |
| bsgs | 44 | 4.462 (3) | 0.2265 | 0.1142 | 0.0256 | 0.0508 | true |
| bsgs | 48 | 18.495 (3) | 1.0429 | 0.5543 | 0.0300 | 0.0564 | true |
| bsgs | 52 | 27.749 (3) | 4.6114 | 2.6462 | 0.0954 | 0.1662 | true |
| rho | 32 | 0.032 (3) | 0.0312 | 0.0053 | 0.1669 | 0.9740 | true |
| rho | 36 | 0.141 (3) | 0.0532 | 0.0142 | 0.1007 | 0.3770 | true |
| rho | 40 | 0.679 (3) | 0.1071 | 0.0588 | 0.0866 | 0.1577 | true |
| rho | 44 | 3.162 (3) | 0.4518 | 0.3894 | 0.1231 | 0.1429 | true |
| rho | 48 | 11.754 (3) | 1.1729 | 1.0332 | 0.0879 | 0.0998 | true |
| rho | 52 | 80.788 (3) | 1.7896 | 1.6480 | 0.0204 | 0.0222 | true |
growth bsgs 40→52 bits, log2(time ratio)/bit: PARI 0.399, A 0.529, B 0.499
growth rho 32→52 bits, log2(time ratio)/bit: PARI 0.565, A 0.292, B 0.413

field-operation counts: 120 B rows compared with A, 0 mismatches
