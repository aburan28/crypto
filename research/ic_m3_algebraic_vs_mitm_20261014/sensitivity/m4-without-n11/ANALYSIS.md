# ANALYSIS (written by `ic_algebraic_floor analyze`)

Cost per relation in wall-clock ns (log2), one core, same binary and host. `measured` = every target finished; `lower_bound` = some target censored at its budget or the cell stopped early; intervals are 95% bootstrap percentiles.

## m = 4

### Cells and the floor

| n | d | \|F\| points | p (decomposable) | decomps/target | exhaustive = MITM | floor amortized log2 ns [95%] | floor online log2 ns |
|--:|--:|--:|--:|--:|--:|--:|--:|
| 13 | 4 | 19 | 0.562 | 1.008 | 12/12 | 14.95 [14.90, 15.02] | 14.18 |
| 17 | 5 | 29 | 0.181 | 0.223 | 12/12 | 17.04 [16.87, 17.23] | 16.89 |
| 19 | 5 | 35 | 0.093 | 0.103 | 12/12 | 18.51 [18.26, 18.79] | 18.44 |
| 23 | 6 | 53 | 0.026 | 0.026 | 12/12 | 21.27 [20.80, 21.87] | 21.26 |
| 29 | 8 | 253 | 0.281 | 0.333 | 0/0 | 21.78 [21.62, 21.94] | 21.78 |
| 31 | 8 | 261 | 0.082 | 0.085 | 0/0 | 23.77 [23.48, 24.08] | 23.76 |

Floor slope (log2 ns per unit n, amortized): **0.456** [0.443, 0.470]; online 0.489; counted (group additions, amortized) 0.539.

### Arms

| arm | n | d | status | targets run (yes/no) | censored | log2 ns per relation [95%] | ÷ floor | counted per relation |
|---|--:|--:|---|--:|--:|--:|--:|--:|
| f4-chain | 13 | 4 | measured | 6/6 | 0 | 26.04 [25.82, 26.23] | 2.166e3 | 3.846e6 word_xors |
| f4-chain | 17 | 5 | measured | 6/6 | 0 | 30.38 [30.16, 30.62] | 1.036e4 | 1.519e8 word_xors |
| f4-chain | 19 | 5 | measured | 6/6 | 0 | 32.12 [31.83, 32.44] | 1.252e4 | 6.665e8 word_xors |
| f4-chain | 23 | 6 | unsupported | 1/0 | 0 | — | — | — |
| f4-chain | 29 | 8 | unsupported | 1/0 | 0 | — | — | — |
| f4-chain | 31 | 8 | unsupported | 1/0 | 0 | — | — | — |

**f4-chain**: measured cells 3, lower-bound cells 0; slope on measured cells 1.025 (floor on the same cells 0.583), Δ = arm − floor 95% [0.412, 0.478]; lower-bound fit slope 1.025 (not a valid bound), Δ 95% [0.412, 0.478]. **Verdict: grows faster than the floor.**

| arm | n | d | status | targets run (yes/no) | censored | log2 ns per relation [95%] | ÷ floor | counted per relation |
|---|--:|--:|---|--:|--:|--:|--:|--:|
| f6-ic | 13 | 4 | measured | 6/6 | 0 | 25.88 [25.66, 26.11] | 1.948e3 | 3.846e6 word_xors |
| f6-ic | 17 | 5 | measured | 6/6 | 0 | 30.26 [30.04, 30.51] | 9.551e3 | 1.504e8 word_xors |
| f6-ic | 19 | 5 | measured | 6/6 | 0 | 32.11 [31.82, 32.42] | 1.236e4 | 6.630e8 word_xors |
| f6-ic | 23 | 6 | unsupported | 1/0 | 0 | — | — | — |
| f6-ic | 29 | 8 | unsupported | 1/0 | 0 | — | — | — |
| f6-ic | 31 | 8 | unsupported | 1/0 | 0 | — | — | — |

**f6-ic**: measured cells 3, lower-bound cells 0; slope on measured cells 1.045 (floor on the same cells 0.583), Δ = arm − floor 95% [0.424, 0.499]; lower-bound fit slope 1.045 (not a valid bound), Δ 95% [0.424, 0.499]. **Verdict: grows faster than the floor.**

| arm | n | d | status | targets run (yes/no) | censored | log2 ns per relation [95%] | ÷ floor | counted per relation |
|---|--:|--:|---|--:|--:|--:|--:|--:|
| s22-hybrid-g8 | 13 | 4 | lower_bound | 6/6 | 6 | ≥ 36.98 | ≥ 4.263e6 | — |
| s22-hybrid-g8 | 17 | 5 | lower_bound | 2/2 | 4 | ≥ 39.27 | ≥ 4.929e6 | — |
| s22-hybrid-g8 | 19 | 5 | not_run | 0/0 | 0 | — | — | — |
| s22-hybrid-g8 | 23 | 6 | not_run | 0/0 | 0 | — | — | — |
| s22-hybrid-g8 | 29 | 8 | not_run | 0/0 | 0 | — | — | — |
| s22-hybrid-g8 | 31 | 8 | not_run | 0/0 | 0 | — | — | — |

**s22-hybrid-g8**: measured cells 0, lower-bound cells 2; slope on measured cells — (floor on the same cells —), Δ = arm − floor 95% —; lower-bound fit slope 0.574 (not a valid bound), Δ 95% [-0.005, 0.107]. **Verdict: undetermined (fewer than 3 usable cells).**

| arm | n | d | status | targets run (yes/no) | censored | log2 ns per relation [95%] | ÷ floor | counted per relation |
|---|--:|--:|---|--:|--:|--:|--:|--:|
| sat | 13 | 4 | lower_bound | 6/6 | 1 | ≥ 35.40 | ≥ 1.424e6 | — |
| sat | 17 | 5 | lower_bound | 6/6 | 9 | ≥ 39.18 | ≥ 4.614e6 | — |
| sat | 19 | 5 | lower_bound | 2/2 | 4 | ≥ 40.24 | ≥ 3.460e6 | — |
| sat | 23 | 6 | not_run | 0/0 | 0 | — | — | — |
| sat | 29 | 8 | not_run | 0/0 | 0 | — | — | — |
| sat | 31 | 8 | not_run | 0/0 | 0 | — | — | — |

**sat**: measured cells 0, lower-bound cells 3; slope on measured cells — (floor on the same cells —), Δ = arm − floor 95% —; lower-bound fit slope 0.826 (not a valid bound), Δ 95% [0.135, 0.350]. **Verdict: undetermined (fewer than 3 usable cells).**

