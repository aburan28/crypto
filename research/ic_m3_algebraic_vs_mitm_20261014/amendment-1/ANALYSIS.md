# ANALYSIS (written by `ic_algebraic_floor analyze`)

Cost per relation in wall-clock ns (log2), one core, same binary and host. `measured` = every target finished; `lower_bound` = some target censored at its budget or the cell stopped early; intervals are 95% bootstrap percentiles.

## m = 4

### Cells and the floor

| n | d | \|F\| points | p (decomposable) | decomps/target | exhaustive = MITM | floor amortized log2 ns [95%] | floor online log2 ns |
|--:|--:|--:|--:|--:|--:|--:|--:|
| 11 | 3 | 3 | 0.009 | 0.017 | 12/12 | 18.09 [17.84, 18.69] | 16.77 |
| 13 | 4 | 19 | 0.562 | 1.008 | 12/12 | 14.99 [14.93, 15.07] | 14.34 |
| 17 | 5 | 29 | 0.181 | 0.223 | 12/12 | 16.96 [16.79, 17.16] | 16.83 |
| 19 | 5 | 35 | 0.093 | 0.103 | 12/12 | 18.37 [18.11, 18.67] | 18.33 |
| 23 | 6 | 53 | 0.026 | 0.026 | 12/12 | 21.06 [20.57, 21.73] | 21.06 |
| 29 | 8 | 253 | 0.281 | 0.333 | 0/0 | 21.68 [21.53, 21.83] | 21.67 |
| 31 | 8 | 261 | 0.082 | 0.085 | 0/0 | 23.48 [23.20, 23.79] | 23.48 |

Floor slope (log2 ns per unit n, amortized): **0.348** [0.328, 0.362]; online 0.399; counted (group additions, amortized) 0.499.

### Arms

| arm | n | d | status | targets run (yes/no) | censored | log2 ns per relation [95%] | ÷ floor | counted per relation |
|---|--:|--:|---|--:|--:|--:|--:|--:|
| wide-chain | 11 | 3 | measured | 6/6 | 0 | 27.20 [26.45, 28.42] | 5.501e2 | 9.101e2 search_nodes |
| wide-chain | 13 | 4 | measured | 6/6 | 0 | 26.05 [25.81, 26.27] | 2.128e3 | 3.495e2 search_nodes |
| wide-chain | 17 | 5 | measured | 6/6 | 0 | 29.85 [29.62, 30.10] | 7.596e3 | 3.836e3 search_nodes |
| wide-chain | 19 | 5 | measured | 6/6 | 0 | 32.27 [31.99, 32.59] | 1.529e4 | 1.349e4 search_nodes |
| wide-chain | 23 | 6 | measured | 6/6 | 0 | 35.96 [35.45, 36.63] | 3.052e4 | 1.509e5 search_nodes |
| wide-chain | 29 | 8 | lower_bound | 6/6 | 7 | ≥ 38.58 | ≥ 1.223e5 | — |
| wide-chain | 31 | 8 | lower_bound | 6/6 | 8 | ≥ 40.52 | ≥ 1.346e5 | — |

**wide-chain**: measured cells 5, lower-bound cells 2; slope on measured cells 0.805 (floor on the same cells 0.333), Δ = arm − floor 95% [0.433, 0.504]; lower-bound fit slope 0.724 (valid as a lower bound), Δ 95% [0.358, 0.390]. **Verdict: grows faster than the floor.**

