# ANALYSIS (written by `ic_algebraic_floor analyze`)

Cost per relation in wall-clock ns (log2), one core, same binary and host. `measured` = every target finished; `lower_bound` = some target censored at its budget or the cell stopped early; intervals are 95% bootstrap percentiles.

## m = 3

### Cells and the floor

| n | d | \|F\| points | p (decomposable) | decomps/target | exhaustive = MITM | floor amortized log2 ns [95%] | floor online log2 ns |
|--:|--:|--:|--:|--:|--:|--:|--:|
| 13 | 5 | 33 | 0.494 | 0.770 | 12/12 | 11.56 [11.46, 11.66] | 11.56 |
| 17 | 6 | 63 | 0.274 | 0.331 | 12/12 | 13.33 [13.22, 13.45] | 13.33 |
| 19 | 7 | 117 | 0.404 | 0.537 | 12/12 | 13.43 [13.34, 13.51] | 13.43 |
| 23 | 8 | 237 | 0.247 | 0.284 | 12/12 | 14.92 [14.81, 15.04] | 14.92 |
| 29 | 10 | 1017 | 0.279 | 0.323 | 0/0 | 16.86 [16.76, 16.97] | 16.86 |
| 31 | 11 | 2101 | 0.479 | 0.678 | 0/0 | 16.56 [16.42, 16.73] | 16.56 |

Floor slope (log2 ns per unit n, amortized): **0.289** [0.281, 0.297]; online 0.289; counted (group additions, amortized) 0.339.

### Arms

| arm | n | d | status | targets run (yes/no) | censored | log2 ns per relation [95%] | ÷ floor | counted per relation |
|---|--:|--:|---|--:|--:|--:|--:|--:|
| f4-chain | 13 | 5 | measured | 6/6 | 0 | 23.56 [23.35, 23.77] | 4.121e3 | 6.224e5 word_xors |
| f4-chain | 17 | 6 | measured | 6/6 | 0 | 26.12 [25.94, 26.28] | 7.042e3 | 7.611e6 word_xors |
| f4-chain | 19 | 7 | measured | 6/6 | 0 | 27.28 [27.04, 27.53] | 1.481e4 | 1.590e7 word_xors |
| f4-chain | 23 | 8 | measured | 6/6 | 0 | 29.55 [29.39, 29.71] | 2.535e4 | 1.363e8 word_xors |
| f4-chain | 29 | 10 | measured | 6/6 | 0 | 32.66 [32.49, 32.85] | 5.726e4 | 1.644e9 word_xors |
| f4-chain | 31 | 11 | measured | 6/6 | 0 | 33.08 [32.90, 33.26] | 9.426e4 | 2.668e9 word_xors |

**f4-chain**: measured cells 6, lower-bound cells 0; slope on measured cells 0.534 (floor on the same cells 0.289), Δ = arm − floor 95% [0.233, 0.257]; lower-bound fit slope 0.534 (not a valid bound), Δ 95% [0.233, 0.257]. **Verdict: grows faster than the floor.**

| arm | n | d | status | targets run (yes/no) | censored | log2 ns per relation [95%] | ÷ floor | counted per relation |
|---|--:|--:|---|--:|--:|--:|--:|--:|
| f6-ic | 13 | 5 | measured | 6/6 | 0 | 22.71 [22.55, 22.86] | 2.284e3 | 5.119e5 word_xors |
| f6-ic | 17 | 6 | measured | 6/6 | 0 | 25.60 [25.39, 25.83] | 4.926e3 | 6.521e6 word_xors |
| f6-ic | 19 | 7 | measured | 6/6 | 0 | 26.38 [26.23, 26.52] | 7.931e3 | 1.261e7 word_xors |
| f6-ic | 23 | 8 | measured | 6/6 | 0 | 28.83 [28.68, 28.98] | 1.539e4 | 1.104e8 word_xors |
| f6-ic | 29 | 10 | measured | 6/6 | 0 | 32.68 [32.51, 32.83] | 5.775e4 | 1.632e9 word_xors |
| f6-ic | 31 | 11 | measured | 6/6 | 0 | 33.04 [32.89, 33.17] | 9.150e4 | 2.654e9 word_xors |

**f6-ic**: measured cells 6, lower-bound cells 0; slope on measured cells 0.583 (floor on the same cells 0.289), Δ = arm − floor 95% [0.284, 0.304]; lower-bound fit slope 0.583 (not a valid bound), Δ 95% [0.284, 0.304]. **Verdict: grows faster than the floor.**

| arm | n | d | status | targets run (yes/no) | censored | log2 ns per relation [95%] | ÷ floor | counted per relation |
|---|--:|--:|---|--:|--:|--:|--:|--:|
| s4-hybrid-f10 | 13 | 5 | measured | 6/6 | 0 | 27.49 [27.28, 27.71] | 6.283e4 | 1.476e7 word_xors |
| s4-hybrid-f10 | 17 | 6 | measured | 6/6 | 0 | 31.77 [31.59, 31.93] | 3.545e5 | 9.980e7 word_xors |
| s4-hybrid-f10 | 19 | 7 | measured | 6/6 | 0 | 35.66 [35.53, 35.80] | 4.936e6 | 6.372e8 word_xors |
| s4-hybrid-f10 | 23 | 8 | lower_bound | 6/6 | 6 | ≥ 38.56 | ≥ 1.308e7 | — |
| s4-hybrid-f10 | 29 | 10 | lower_bound | 2/2 | 4 | ≥ 38.81 | ≥ 4.068e6 | — |
| s4-hybrid-f10 | 31 | 11 | not_run | 0/0 | 0 | — | — | — |

**s4-hybrid-f10**: measured cells 3, lower-bound cells 2; slope on measured cells 1.320 (floor on the same cells 0.331), Δ = arm − floor 95% [0.950, 1.026]; lower-bound fit slope 0.720 (valid as a lower bound), Δ 95% [0.386, 0.407]. **Verdict: grows faster than the floor.**

| arm | n | d | status | targets run (yes/no) | censored | log2 ns per relation [95%] | ÷ floor | counted per relation |
|---|--:|--:|---|--:|--:|--:|--:|--:|
| s4-hybrid-f12 | 13 | 5 | measured | 6/6 | 0 | 30.12 [29.90, 30.29] | 3.872e5 | 6.533e8 word_xors |
| s4-hybrid-f12 | 17 | 6 | measured | 6/6 | 0 | 32.77 [32.60, 32.94] | 7.107e5 | 2.700e9 word_xors |
| s4-hybrid-f12 | 19 | 7 | measured | 6/6 | 0 | 34.24 [34.10, 34.38] | 1.845e6 | 1.126e9 word_xors |
| s4-hybrid-f12 | 23 | 8 | lower_bound | 6/6 | 7 | ≥ 38.64 | ≥ 1.378e7 | — |
| s4-hybrid-f12 | 29 | 10 | lower_bound | 2/2 | 4 | ≥ 38.81 | ≥ 4.067e6 | — |
| s4-hybrid-f12 | 31 | 11 | not_run | 0/0 | 0 | — | — | — |

**s4-hybrid-f12**: measured cells 3, lower-bound cells 2; slope on measured cells 0.684 (floor on the same cells 0.331), Δ = arm − floor 95% [0.318, 0.392]; lower-bound fit slope 0.584 (valid as a lower bound), Δ 95% [0.251, 0.272]. **Verdict: grows faster than the floor.**

| arm | n | d | status | targets run (yes/no) | censored | log2 ns per relation [95%] | ÷ floor | counted per relation |
|---|--:|--:|---|--:|--:|--:|--:|--:|
| sat | 13 | 5 | measured | 6/6 | 0 | 30.21 [30.02, 30.39] | 4.117e5 | 5.973e4 conflicts |
| sat | 17 | 6 | lower_bound | 6/6 | 3 | ≥ 37.73 | ≥ 2.209e7 | — |
| sat | 19 | 7 | lower_bound | 6/6 | 6 | ≥ 37.53 | ≥ 1.803e7 | — |
| sat | 23 | 8 | lower_bound | 2/2 | 4 | ≥ 38.82 | ≥ 1.569e7 | — |
| sat | 29 | 10 | not_run | 0/0 | 0 | — | — | — |
| sat | 31 | 11 | not_run | 0/0 | 0 | — | — | — |

**sat**: measured cells 1, lower-bound cells 3; slope on measured cells — (floor on the same cells —), Δ = arm − floor 95% —; lower-bound fit slope 0.825 (not a valid bound), Δ 95% [0.480, 0.521]. **Verdict: undetermined (fewer than 3 usable cells).**

## m = 4

### Cells and the floor

| n | d | \|F\| points | p (decomposable) | decomps/target | exhaustive = MITM | floor amortized log2 ns [95%] | floor online log2 ns |
|--:|--:|--:|--:|--:|--:|--:|--:|
| 11 | 3 | 3 | 0.009 | 0.017 | 12/12 | 16.87 [16.04, 18.13] | 16.87 |
| 13 | 4 | 19 | 0.562 | 1.008 | 12/12 | 14.18 [14.08, 14.28] | 14.18 |
| 17 | 5 | 29 | 0.181 | 0.223 | 12/12 | 16.89 [16.70, 17.11] | 16.89 |
| 19 | 5 | 35 | 0.093 | 0.103 | 12/12 | 18.44 [18.18, 18.76] | 18.44 |
| 23 | 6 | 53 | 0.026 | 0.026 | 12/12 | 21.26 [20.80, 21.85] | 21.26 |
| 29 | 8 | 253 | 0.281 | 0.333 | 0/0 | 21.78 [21.62, 21.94] | 21.78 |
| 31 | 8 | 261 | 0.082 | 0.085 | 0/0 | 23.76 [23.48, 24.09] | 23.76 |

Floor slope (log2 ns per unit n, amortized): **0.412** [0.376, 0.437]; online 0.412; counted (group additions, amortized) 0.500.

### Arms

| arm | n | d | status | targets run (yes/no) | censored | log2 ns per relation [95%] | ÷ floor | counted per relation |
|---|--:|--:|---|--:|--:|--:|--:|--:|
| f4-chain | 11 | 3 | measured | 6/6 | 0 | 28.70 [27.91, 29.94] | 3.646e3 | 1.831e7 word_xors |
| f4-chain | 13 | 4 | measured | 6/6 | 0 | 26.04 [25.81, 26.23] | 3.718e3 | 3.846e6 word_xors |
| f4-chain | 17 | 5 | measured | 6/6 | 0 | 30.38 [30.16, 30.63] | 1.149e4 | 1.519e8 word_xors |
| f4-chain | 19 | 5 | measured | 6/6 | 0 | 32.12 [31.85, 32.45] | 1.312e4 | 6.665e8 word_xors |
| f4-chain | 23 | 6 | unsupported | 1/0 | 0 | — | — | — |
| f4-chain | 29 | 8 | unsupported | 1/0 | 0 | — | — | — |
| f4-chain | 31 | 8 | unsupported | 1/0 | 0 | — | — | — |

**f4-chain**: measured cells 4, lower-bound cells 0; slope on measured cells 0.560 (floor on the same cells 0.294), Δ = arm − floor 95% [0.239, 0.296]; lower-bound fit slope 0.560 (not a valid bound), Δ 95% [0.239, 0.296]. **Verdict: grows faster than the floor.**

| arm | n | d | status | targets run (yes/no) | censored | log2 ns per relation [95%] | ÷ floor | counted per relation |
|---|--:|--:|---|--:|--:|--:|--:|--:|
| f6-ic | 11 | 3 | measured | 6/6 | 0 | 28.56 [27.78, 29.80] | 3.316e3 | 1.744e7 word_xors |
| f6-ic | 13 | 4 | measured | 6/6 | 0 | 25.88 [25.66, 26.11] | 3.344e3 | 3.846e6 word_xors |
| f6-ic | 17 | 5 | measured | 6/6 | 0 | 30.26 [30.03, 30.51] | 1.059e4 | 1.504e8 word_xors |
| f6-ic | 19 | 5 | measured | 6/6 | 0 | 32.11 [31.83, 32.43] | 1.295e4 | 6.630e8 word_xors |
| f6-ic | 23 | 6 | unsupported | 1/0 | 0 | — | — | — |
| f6-ic | 29 | 8 | unsupported | 1/0 | 0 | — | — | — |
| f6-ic | 31 | 8 | unsupported | 1/0 | 0 | — | — | — |

**f6-ic**: measured cells 4, lower-bound cells 0; slope on measured cells 0.573 (floor on the same cells 0.294), Δ = arm − floor 95% [0.254, 0.309]; lower-bound fit slope 0.573 (not a valid bound), Δ 95% [0.254, 0.309]. **Verdict: grows faster than the floor.**

| arm | n | d | status | targets run (yes/no) | censored | log2 ns per relation [95%] | ÷ floor | counted per relation |
|---|--:|--:|---|--:|--:|--:|--:|--:|
| s22-hybrid-g8 | 11 | 3 | measured | 6/6 | 0 | 35.49 [34.72, 36.69] | 4.050e5 | 3.061e8 word_xors |
| s22-hybrid-g8 | 13 | 4 | lower_bound | 6/6 | 6 | ≥ 36.98 | ≥ 7.318e6 | — |
| s22-hybrid-g8 | 17 | 5 | lower_bound | 2/2 | 4 | ≥ 39.27 | ≥ 5.466e6 | — |
| s22-hybrid-g8 | 19 | 5 | not_run | 0/0 | 0 | — | — | — |
| s22-hybrid-g8 | 23 | 6 | not_run | 0/0 | 0 | — | — | — |
| s22-hybrid-g8 | 29 | 8 | not_run | 0/0 | 0 | — | — | — |
| s22-hybrid-g8 | 31 | 8 | not_run | 0/0 | 0 | — | — | — |

**s22-hybrid-g8**: measured cells 1, lower-bound cells 2; slope on measured cells — (floor on the same cells —), Δ = arm − floor 95% —; lower-bound fit slope 0.622 (not a valid bound), Δ 95% [0.496, 0.557]. **Verdict: undetermined (fewer than 3 usable cells).**

| arm | n | d | status | targets run (yes/no) | censored | log2 ns per relation [95%] | ÷ floor | counted per relation |
|---|--:|--:|---|--:|--:|--:|--:|--:|
| sat | 11 | 3 | measured | 6/6 | 0 | 31.84 [31.02, 33.19] | 3.207e4 | 2.069e5 conflicts |
| sat | 13 | 4 | lower_bound | 6/6 | 1 | ≥ 35.40 | ≥ 2.445e6 | — |
| sat | 17 | 5 | lower_bound | 6/6 | 9 | ≥ 39.18 | ≥ 5.117e6 | — |
| sat | 19 | 5 | lower_bound | 2/2 | 4 | ≥ 40.24 | ≥ 3.627e6 | — |
| sat | 23 | 6 | not_run | 0/0 | 0 | — | — | — |
| sat | 29 | 8 | not_run | 0/0 | 0 | — | — | — |
| sat | 31 | 8 | not_run | 0/0 | 0 | — | — | — |

**sat**: measured cells 1, lower-bound cells 3; slope on measured cells — (floor on the same cells —), Δ = arm − floor 95% —; lower-bound fit slope 1.029 (not a valid bound), Δ 95% [0.687, 0.787]. **Verdict: undetermined (fewer than 3 usable cells).**

