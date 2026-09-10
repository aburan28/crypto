# Stage 1: matched planted-PDP solver matrix

This stage defines every factor base from public field algebra. No run enumerates the target subgroup or constructs factor-base discrete-log labels. Each cell exports one source system to WDSat ANF, CryptoMiniSat CNF-XOR, and Magma Boolean F4; native SAT and direct meet-in-the-middle use the same target and factor-base predicate.

The table reports solver-kernel wall time. WDSat and CryptoMiniSat also have isolated total core-seconds and peak RSS in `stage-1-composed-summary.json`. Native SAT and MITM have exact internal wall clocks, while their process CPU/RSS are combined in the exporter receipt; that accounting split is a remaining gate, not silently estimated here.

| n | ell | base | vars/eqs | native result, conflicts, s | WDSat result, conflicts, s | CMS result, conflicts, s | MITM result, s |
|--:|--:|:--|--:|:--|:--|:--|:--|
| 15 | 5 | standard | 42/42 | sat, 2731, 0.056453 | sat, 8, 0.447321 | sat, 0, 0.010865 | sat, 0.003889 |
| 31 | 5 | standard | 42/58 | sat, 244, 0.009184 | sat, 6615, 0.244580 | sat, 6790, 0.074314 | sat, 0.004038 |
| 31 | 5 | ggmp | 46/62 | sat, 319, 0.004193 | sat, 12, 0.217487 | sat, 0, 0.012488 | sat, 0.000347 |
| 41 | 5 | standard | 42/68 | sat, 1981, 0.046369 | sat, 6392, 0.250815 | sat, 28099, 0.303335 | sat, 0.004637 |
| 59 | 9 | standard | 78/110 | unknown, 100000, 8.918730 | timeout_inconclusive, -, 120.005169 | timeout_inconclusive, -, 120.014098 | sat, 3.126873 |

| n | base | combined exporter core-s / peak MiB | WDSat core-s / peak MiB | WDSat build core-s / peak MiB | CMS core-s / peak MiB |
|--:|:--|--:|--:|--:|--:|
| 15 | standard | 0.067911 / 5.44 | 0.002128 / 3.30 | 0.623483 / 56.56 | 0.009723 / 4.52 |
| 31 | standard | 0.021561 / 5.22 | 0.037935 / 3.50 | 0.558952 / 56.38 | 0.072967 / 6.03 |
| 31 | ggmp | 0.009990 / 4.22 | 0.002812 / 7.02 | 0.538840 / 56.12 | 0.011274 / 5.78 |
| 41 | standard | 0.060848 / 6.42 | 0.042600 / 3.59 | 0.539471 / 56.64 | 0.301732 / 8.50 |
| 59 | standard | 12.091462 / 188.80 | 117.999105 / 13.59 | 0.512332 / 54.95 | 119.261571 / 171.03 |

At n=31 and n=41, every available solver returned a source-validated SAT model. At n=59, direct MITM found the planted decomposition in about 3.13 seconds; native SAT stopped at 100,000 conflicts, WDSat reached the 120-second watchdog, and CryptoMiniSat reached the same watchdog. Both capped solver outcomes are inconclusive.

The n=31 GGMP cell derives its five-dimensional factor base as the kernel of a linearised polynomial obtained from a factor of T^31-1. It is an admissible construction, but this single planted target is not a causal comparison with the standard-basis cell because the two cells use different factor-base predicates and targets.

Magma is not installed on this host, so the generated `.magma` inputs are retained but F4 is an operationally missing baseline. The n=67 cell exceeds the current n<64 field-bitmask implementation and asserts nothing about PDP hardness.

This clears only the one-instance planted-PDP scaling sub-gate through n=59. Unknown-scalar end-to-end runs, fully isolated process accounting for every arm, Magma F4, automorphism-optimized Pollard rho, repeated scaling, external reproduction, and novelty review remain open. The result is not a Koblitz index-calculus SOTA claim.

Primary comparisons: Trimoska-Ionica-Dequen, *A SAT-Based Approach for Index Calculus on Binary Elliptic Curves* (ePrint 2019/313); Galbraith-Granger-Merz-Petit, *On Index Calculus Algorithms for Subfield Curves* (ePrint 2020/1315).
