# ANALYSIS (written by `ic_algebraic_floor analyze`)

Cost per relation in wall-clock ns (log2), one core, same binary and host. `measured` = every target finished; `lower_bound` = some target censored at its budget or the cell stopped early; intervals are 95% bootstrap percentiles.

## m = 3

### Cells and the floor

| n | d | \|F\| points | p (decomposable) | decomps/target | exhaustive = MITM | floor amortized log2 ns [95%] | floor online log2 ns |
|--:|--:|--:|--:|--:|--:|--:|--:|
| 13 | 5 | 33 | 0.494 | 0.770 | 12/12 | 14.47 [14.46, 14.48] | 11.56 |
| 17 | 6 | 63 | 0.274 | 0.331 | 12/12 | 14.40 [14.34, 14.46] | 13.33 |
| 19 | 7 | 117 | 0.404 | 0.537 | 12/12 | 14.76 [14.73, 14.80] | 13.43 |
| 23 | 8 | 237 | 0.247 | 0.284 | 12/12 | 15.47 [15.40, 15.56] | 14.92 |
| 29 | 10 | 1017 | 0.279 | 0.323 | 0/0 | 17.38 [17.31, 17.46] | 16.86 |
| 31 | 11 | 2101 | 0.479 | 0.678 | 0/0 | 17.82 [17.77, 17.90] | 16.56 |

Floor slope (log2 ns per unit n, amortized): **0.207** [0.204, 0.211]; online 0.289; counted (group additions, amortized) 0.337.

### Arms

| arm | n | d | status | targets run (yes/no) | censored | log2 ns per relation [95%] | ÷ floor | counted per relation |
|---|--:|--:|---|--:|--:|--:|--:|--:|
| f4-chain | 13 | 5 | measured | 6/6 | 0 | 23.56 [23.35, 23.77] | 5.471e2 | 6.224e5 word_xors |
| f4-chain | 17 | 6 | measured | 6/6 | 0 | 26.12 [25.94, 26.28] | 3.366e3 | 7.611e6 word_xors |
| f4-chain | 19 | 7 | measured | 6/6 | 0 | 27.28 [27.04, 27.53] | 5.872e3 | 1.590e7 word_xors |
| f4-chain | 23 | 8 | measured | 6/6 | 0 | 29.55 [29.39, 29.71] | 1.727e4 | 1.363e8 word_xors |
| f4-chain | 29 | 10 | measured | 6/6 | 0 | 32.66 [32.49, 32.85] | 3.989e4 | 1.644e9 word_xors |
| f4-chain | 31 | 11 | measured | 6/6 | 0 | 33.08 [32.90, 33.26] | 3.924e4 | 2.668e9 word_xors |

**f4-chain**: measured cells 6, lower-bound cells 0; slope on measured cells 0.534 (floor on the same cells 0.207), Δ = arm − floor 95% [0.315, 0.338]; lower-bound fit slope 0.534 (not a valid bound), Δ 95% [0.315, 0.338]. **Verdict: grows faster than the floor.**

| arm | n | d | status | targets run (yes/no) | censored | log2 ns per relation [95%] | ÷ floor | counted per relation |
|---|--:|--:|---|--:|--:|--:|--:|--:|
| f6-ic | 13 | 5 | measured | 6/6 | 0 | 22.71 [22.55, 22.86] | 3.033e2 | 5.119e5 word_xors |
| f6-ic | 17 | 6 | measured | 6/6 | 0 | 25.60 [25.39, 25.83] | 2.355e3 | 6.521e6 word_xors |
| f6-ic | 19 | 7 | measured | 6/6 | 0 | 26.38 [26.23, 26.52] | 3.144e3 | 1.261e7 word_xors |
| f6-ic | 23 | 8 | measured | 6/6 | 0 | 28.83 [28.68, 28.98] | 1.049e4 | 1.104e8 word_xors |
| f6-ic | 29 | 10 | measured | 6/6 | 0 | 32.68 [32.51, 32.83] | 4.023e4 | 1.632e9 word_xors |
| f6-ic | 31 | 11 | measured | 6/6 | 0 | 33.04 [32.89, 33.17] | 3.809e4 | 2.654e9 word_xors |

**f6-ic**: measured cells 6, lower-bound cells 0; slope on measured cells 0.583 (floor on the same cells 0.207), Δ = arm − floor 95% [0.366, 0.384]; lower-bound fit slope 0.583 (not a valid bound), Δ 95% [0.366, 0.384]. **Verdict: grows faster than the floor.**

| arm | n | d | status | targets run (yes/no) | censored | log2 ns per relation [95%] | ÷ floor | counted per relation |
|---|--:|--:|---|--:|--:|--:|--:|--:|
| s4-hybrid-f10 | 13 | 5 | measured | 6/6 | 0 | 27.49 [27.28, 27.71] | 8.342e3 | 1.476e7 word_xors |
| s4-hybrid-f10 | 17 | 6 | measured | 6/6 | 0 | 31.77 [31.59, 31.93] | 1.695e5 | 9.980e7 word_xors |
| s4-hybrid-f10 | 19 | 7 | measured | 6/6 | 0 | 35.66 [35.53, 35.80] | 1.957e6 | 6.372e8 word_xors |
| s4-hybrid-f10 | 23 | 8 | lower_bound | 6/6 | 6 | ≥ 38.56 | ≥ 8.910e6 | — |
| s4-hybrid-f10 | 29 | 10 | lower_bound | 2/2 | 4 | ≥ 38.81 | ≥ 2.834e6 | — |
| s4-hybrid-f10 | 31 | 11 | not_run | 0/0 | 0 | — | — | — |

**s4-hybrid-f10**: measured cells 3, lower-bound cells 2; slope on measured cells 1.320 (floor on the same cells 0.040), Δ = arm − floor 95% [1.240, 1.321]; lower-bound fit slope 0.720 (valid as a lower bound), Δ 95% [0.519, 0.541]. **Verdict: grows faster than the floor.**

| arm | n | d | status | targets run (yes/no) | censored | log2 ns per relation [95%] | ÷ floor | counted per relation |
|---|--:|--:|---|--:|--:|--:|--:|--:|
| s4-hybrid-f12 | 13 | 5 | measured | 6/6 | 0 | 30.12 [29.90, 30.29] | 5.140e4 | 6.533e8 word_xors |
| s4-hybrid-f12 | 17 | 6 | measured | 6/6 | 0 | 32.77 [32.60, 32.94] | 3.397e5 | 2.700e9 word_xors |
| s4-hybrid-f12 | 19 | 7 | measured | 6/6 | 0 | 34.24 [34.10, 34.38] | 7.316e5 | 1.126e9 word_xors |
| s4-hybrid-f12 | 23 | 8 | lower_bound | 6/6 | 7 | ≥ 38.64 | ≥ 9.392e6 | — |
| s4-hybrid-f12 | 29 | 10 | lower_bound | 2/2 | 4 | ≥ 38.81 | ≥ 2.834e6 | — |
| s4-hybrid-f12 | 31 | 11 | not_run | 0/0 | 0 | — | — | — |

**s4-hybrid-f12**: measured cells 3, lower-bound cells 2; slope on measured cells 0.684 (floor on the same cells 0.040), Δ = arm − floor 95% [0.608, 0.687]; lower-bound fit slope 0.584 (valid as a lower bound), Δ 95% [0.385, 0.405]. **Verdict: grows faster than the floor.**

| arm | n | d | status | targets run (yes/no) | censored | log2 ns per relation [95%] | ÷ floor | counted per relation |
|---|--:|--:|---|--:|--:|--:|--:|--:|
| sat | 13 | 5 | measured | 6/6 | 0 | 30.21 [30.02, 30.39] | 5.467e4 | 5.973e4 conflicts |
| sat | 17 | 6 | lower_bound | 6/6 | 3 | ≥ 37.73 | ≥ 1.056e7 | — |
| sat | 19 | 7 | lower_bound | 6/6 | 6 | ≥ 37.53 | ≥ 7.149e6 | — |
| sat | 23 | 8 | lower_bound | 2/2 | 4 | ≥ 38.82 | ≥ 1.069e7 | — |
| sat | 29 | 10 | not_run | 0/0 | 0 | — | — | — |
| sat | 31 | 11 | not_run | 0/0 | 0 | — | — | — |

**sat**: measured cells 1, lower-bound cells 3; slope on measured cells — (floor on the same cells —), Δ = arm − floor 95% —; lower-bound fit slope 0.825 (not a valid bound), Δ 95% [0.700, 0.744]. **Verdict: undetermined (fewer than 3 usable cells).**

## m = 4

### Cells and the floor

| n | d | \|F\| points | p (decomposable) | decomps/target | exhaustive = MITM | floor amortized log2 ns [95%] | floor online log2 ns |
|--:|--:|--:|--:|--:|--:|--:|--:|
| 11 | 3 | 3 | 0.009 | 0.017 | 12/12 | 17.99 [17.67, 18.71] | 16.87 |
| 13 | 4 | 19 | 0.562 | 1.008 | 12/12 | 14.95 [14.90, 15.02] | 14.18 |
| 17 | 5 | 29 | 0.181 | 0.223 | 12/12 | 17.04 [16.87, 17.24] | 16.89 |
| 19 | 5 | 35 | 0.093 | 0.103 | 12/12 | 18.51 [18.26, 18.81] | 18.44 |
| 23 | 6 | 53 | 0.026 | 0.026 | 12/12 | 21.27 [20.80, 21.86] | 21.26 |
| 29 | 8 | 253 | 0.281 | 0.333 | 0/0 | 21.78 [21.62, 21.94] | 21.78 |
| 31 | 8 | 261 | 0.082 | 0.085 | 0/0 | 23.77 [23.48, 24.09] | 23.76 |

Floor slope (log2 ns per unit n, amortized): **0.363** [0.342, 0.379]; online 0.412; counted (group additions, amortized) 0.499.

### Arms

| arm | n | d | status | targets run (yes/no) | censored | log2 ns per relation [95%] | ÷ floor | counted per relation |
|---|--:|--:|---|--:|--:|--:|--:|--:|
| f4-chain | 11 | 3 | measured | 6/6 | 0 | 28.70 [27.91, 29.94] | 1.672e3 | 1.831e7 word_xors |
| f4-chain | 13 | 4 | measured | 6/6 | 0 | 26.04 [25.81, 26.23] | 2.166e3 | 3.846e6 word_xors |
| f4-chain | 17 | 5 | measured | 6/6 | 0 | 30.38 [30.16, 30.63] | 1.036e4 | 1.519e8 word_xors |
| f4-chain | 19 | 5 | measured | 6/6 | 0 | 32.12 [31.85, 32.45] | 1.252e4 | 6.665e8 word_xors |
| f4-chain | 23 | 6 | unsupported | 1/0 | 0 | — | — | — |
| f4-chain | 29 | 8 | unsupported | 1/0 | 0 | — | — | — |
| f4-chain | 31 | 8 | unsupported | 1/0 | 0 | — | — | — |

**f4-chain**: measured cells 4, lower-bound cells 0; slope on measured cells 0.560 (floor on the same cells 0.156), Δ = arm − floor 95% [0.342, 0.453]; lower-bound fit slope 0.560 (not a valid bound), Δ 95% [0.342, 0.453]. **Verdict: grows faster than the floor.**

| arm | n | d | status | targets run (yes/no) | censored | log2 ns per relation [95%] | ÷ floor | counted per relation |
|---|--:|--:|---|--:|--:|--:|--:|--:|
| f6-ic | 11 | 3 | measured | 6/6 | 0 | 28.56 [27.78, 29.80] | 1.521e3 | 1.744e7 word_xors |
| f6-ic | 13 | 4 | measured | 6/6 | 0 | 25.88 [25.66, 26.11] | 1.948e3 | 3.846e6 word_xors |
| f6-ic | 17 | 5 | measured | 6/6 | 0 | 30.26 [30.03, 30.51] | 9.551e3 | 1.504e8 word_xors |
| f6-ic | 19 | 5 | measured | 6/6 | 0 | 32.11 [31.83, 32.43] | 1.236e4 | 6.630e8 word_xors |
| f6-ic | 23 | 6 | unsupported | 1/0 | 0 | — | — | — |
| f6-ic | 29 | 8 | unsupported | 1/0 | 0 | — | — | — |
| f6-ic | 31 | 8 | unsupported | 1/0 | 0 | — | — | — |

**f6-ic**: measured cells 4, lower-bound cells 0; slope on measured cells 0.573 (floor on the same cells 0.156), Δ = arm − floor 95% [0.354, 0.469]; lower-bound fit slope 0.573 (not a valid bound), Δ 95% [0.354, 0.469]. **Verdict: grows faster than the floor.**

| arm | n | d | status | targets run (yes/no) | censored | log2 ns per relation [95%] | ÷ floor | counted per relation |
|---|--:|--:|---|--:|--:|--:|--:|--:|
| s22-hybrid-g8 | 11 | 3 | measured | 6/6 | 0 | 35.49 [34.72, 36.69] | 1.858e5 | 3.061e8 word_xors |
| s22-hybrid-g8 | 13 | 4 | lower_bound | 6/6 | 6 | ≥ 36.98 | ≥ 4.263e6 | — |
| s22-hybrid-g8 | 17 | 5 | lower_bound | 2/2 | 4 | ≥ 39.27 | ≥ 4.929e6 | — |
| s22-hybrid-g8 | 19 | 5 | not_run | 0/0 | 0 | — | — | — |
| s22-hybrid-g8 | 23 | 6 | not_run | 0/0 | 0 | — | — | — |
| s22-hybrid-g8 | 29 | 8 | not_run | 0/0 | 0 | — | — | — |
| s22-hybrid-g8 | 31 | 8 | not_run | 0/0 | 0 | — | — | — |

**s22-hybrid-g8**: measured cells 1, lower-bound cells 2; slope on measured cells — (floor on the same cells —), Δ = arm − floor 95% —; lower-bound fit slope 0.622 (not a valid bound), Δ 95% [0.599, 0.753]. **Verdict: undetermined (fewer than 3 usable cells).**

| arm | n | d | status | targets run (yes/no) | censored | log2 ns per relation [95%] | ÷ floor | counted per relation |
|---|--:|--:|---|--:|--:|--:|--:|--:|
| sat | 11 | 3 | measured | 6/6 | 0 | 31.84 [31.02, 33.19] | 1.471e4 | 2.069e5 conflicts |
| sat | 13 | 4 | lower_bound | 6/6 | 1 | ≥ 35.40 | ≥ 1.424e6 | — |
| sat | 17 | 5 | lower_bound | 6/6 | 9 | ≥ 39.18 | ≥ 4.614e6 | — |
| sat | 19 | 5 | lower_bound | 2/2 | 4 | ≥ 40.24 | ≥ 3.460e6 | — |
| sat | 23 | 6 | not_run | 0/0 | 0 | — | — | — |
| sat | 29 | 8 | not_run | 0/0 | 0 | — | — | — |
| sat | 31 | 8 | not_run | 0/0 | 0 | — | — | — |

**sat**: measured cells 1, lower-bound cells 3; slope on measured cells — (floor on the same cells —), Δ = arm − floor 95% —; lower-bound fit slope 1.029 (not a valid bound), Δ 95% [0.799, 0.938]. **Verdict: undetermined (fewer than 3 usable cells).**

