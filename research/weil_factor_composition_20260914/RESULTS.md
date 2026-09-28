# Weil-friendly factor-base composition: measured results

Implemented an exact component-coordinate S3 solver with optional linear projection. The saved run contains 864 processes, 912 matched native root-set comparisons (1,824 reference/candidate target executions), 960 complete component root-set checks, and 21,888 full-DLP relation-attempt cross-checks. All completed checks pass. Four new unit tests and the no-default-features library/ic build also pass. Timeouts and incomplete attempts remain in the tables. No calibrated total-cost or exponent advance is established.

The common operation-unit ledger stays explicit. Partial field/XOR/reduction counters cannot price the complete attack, so total operations, `S` and boundary ratios are null for every variant, including rho. Runtime tables below are secondary shared-host diagnostics.

| Variant | Class | Total operations | S | Rho ratio | Floor ratio | Verified workload |
|---|---|---|---|---|---|---|
| ambient_f4 | engineering diagnostic | null | null | null | null | 33 verified DLPs / 120 attempts |
| charts_f4 | engineering diagnostic | null | null | null | null | 480 verified stage queries; 177 witnesses |
| charts_linear | engineering diagnostic | null | null | null | null | 96 verified DLPs / 120 attempts |
| enumerate | engineering diagnostic | null | null | null | null | 96 verified DLPs / 120 attempts |
| pair_table | engineering diagnostic | null | null | null | null | 96 verified DLPs / 120 attempts |
| rho | engineering diagnostic | null | null | null | null | 60 verified DLPs / 60 attempts |

## Process outcomes

| Kind | Finished | Timed out |
|---|---|---|
| native | 84 | 0 |
| stage | 207 | 33 |
| dlp | 393 | 87 |
| rho | 60 | 0 |

A finished DLP process may be incomplete; only an independently verified recovered scalar counts as a solve. Controls with zero projected columns are disqualified factor bases, not speedup wins. Every timed-out process keeps its partial stdout and stderr.

## Actual factor-base census

Coverage counts distinct nonzero subgroup points reachable by two factor-base summands. Same-component coverage is the incomplete ablation. Dimensions and point counts are not interchangeable.

| Case | n/k | Seed family | ell | Seed | Components | Product ranks | Pair-weighted rank | Points | Projected columns | Exact coverage | Same-component coverage | Tensor bytes |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 7/1 | power_span | 2 | 17,937 | 7 | 3–4 | 3.750 | 29 | 2 | 70/70 | 42 | 672 |
| 1 | 9/1 | power_span | 3 | 17,937 | 9 | 5–7 | 6.200 | 3 | 0 | 0/126 | 0 | 1728 |
| 2 | 9/1 | scaled_subfield | 3 | 17,937 | 9 | 3–3 | 3.000 | 55 | 2 | 126/126 | 54 | 1728 |
| 3 | 9/1 | subfield_control | 3 | 17,937 | 1 | 3–3 | 3.000 | 3 | 0 | 0/126 | 0 | 192 |
| 4 | 11/1 | power_span | 3 | 17,937 | 11 | 5–9 | 8.000 | 23 | 1 | 198/990 | 44 | 2112 |
| 5 | 11/1 | random | 3 | 17 | 11 | 6–9 | 8.500 | 89 | 4 | 814/990 | 198 | 2112 |
| 5 | 11/1 | random | 3 | 937 | 11 | 6–9 | 8.333 | 45 | 2 | 352/990 | 66 | 2112 |
| 6 | 13/1 | power_span | 3 | 17,937 | 13 | 5–9 | 8.143 | 81 | 3 | 858/2002 | 182 | 2496 |
| 7 | 15/3 | power_span | 3 | 17,937 | 5 | 5–9 | 7.333 | 33 | 3 | 60/660 | 20 | 960 |
| 8 | 15/3 | scaled_subfield | 3 | 17,937 | 5 | 3–3 | 3.000 | 41 | 4 | 70/660 | 20 | 960 |
| 9 | 15/3 | subfield_control | 5 | 17,937 | 1 | 5–5 | 5.000 | 33 | 3 | 30/660 | 30 | 480 |

The degree-nine ordinary seed and unscaled-subfield control both have zero projected columns. Scaling the subfield produces two useful columns and full coverage of the 126 nonzero subgroup points. At degree fifteen over GF(8), the scaled GF(8) components keep every product rank at three; their coverage is 70/660, versus 60/660 for the ordinary seed. The complementary GF(32) control has three columns and coverage 30/660, so it avoids the zero-column trap but has lower yield. These families change the base and its counting boundary; they are not an equal-base performance ratio.

## Fresh solver workloads

Each completed process covers the same eight targets. Three repetitions use seed 937. Cold milliseconds include setup and first-witness queries; complete root enumeration is separate validation. Witnesses/reductions in timeout rows are partial observations, not comparable totals. F4 reduction counts are exact calls to differently sized systems, not common-cost operation units.

| Case | Variant | Complete processes | Witnesses observed | F4 reductions observed | Median cold ms / 8 targets |
|---|---|---|---|---|---|
| 0 | ambient_f4 | 3/3 | 24 | 1914 | 66.262 |
| 0 | charts_f4 | 3/3 | 24 | 228 | 1.081 |
| 0 | charts_linear | 3/3 | 24 | 57 | 0.728 |
| 0 | enumerate | 3/3 | 24 | 0 | 0.464 |
| 1 | ambient_f4 | 3/3 | 0 | 46368 | 4450.652 |
| 1 | charts_f4 | 3/3 | 0 | 1575 | 15.239 |
| 1 | charts_linear | 3/3 | 0 | 813 | 4.302 |
| 1 | enumerate | 3/3 | 0 | 0 | 0.659 |
| 2 | ambient_f4 | 3/3 | 24 | 9957 | 958.730 |
| 2 | charts_f4 | 3/3 | 24 | 279 | 2.734 |
| 2 | charts_linear | 3/3 | 24 | 42 | 1.167 |
| 2 | enumerate | 3/3 | 24 | 0 | 0.896 |
| 3 | ambient_f4 | 3/3 | 0 | 46368 | 4562.311 |
| 3 | charts_f4 | 3/3 | 0 | 24 | 0.611 |
| 3 | charts_linear | 3/3 | 0 | 0 | 0.483 |
| 3 | enumerate | 3/3 | 0 | 0 | 0.579 |
| 4 | ambient_f4 | 0/3 | 6 | 14553 | — |
| 4 | charts_f4 | 3/3 | 9 | 1761 | 18.959 |
| 4 | charts_linear | 3/3 | 9 | 846 | 5.204 |
| 4 | enumerate | 3/3 | 9 | 0 | 1.425 |
| 5 | ambient_f4 | 0/3 | 0 | 22839 | — |
| 5 | charts_f4 | 3/3 | 3 | 1794 | 23.395 |
| 5 | charts_linear | 3/3 | 3 | 1077 | 6.871 |
| 5 | enumerate | 3/3 | 3 | 0 | 2.370 |
| 6 | ambient_f4 | 0/3 | 3 | 1125 | — |
| 6 | charts_f4 | 3/3 | 12 | 2457 | 16.240 |
| 6 | charts_linear | 3/3 | 12 | 378 | 4.089 |
| 6 | enumerate | 3/3 | 12 | 0 | 3.396 |
| 7 | ambient_f4 | 0/3 | 0 | 0 | — |
| 7 | charts_f4 | 3/3 | 6 | 324 | 3.195 |
| 7 | charts_linear | 3/3 | 6 | 18 | 1.644 |
| 7 | enumerate | 3/3 | 6 | 0 | 2.429 |
| 8 | ambient_f4 | 0/3 | 0 | 0 | — |
| 8 | charts_f4 | 3/3 | 6 | 297 | 4.037 |
| 8 | charts_linear | 3/3 | 6 | 6 | 1.702 |
| 8 | enumerate | 3/3 | 6 | 0 | 3.038 |
| 9 | ambient_f4 | 0/3 | 0 | 0 | — |
| 9 | charts_f4 | 3/3 | 0 | 24 | 1.843 |
| 9 | charts_linear | 3/3 | 0 | 0 | 1.573 |
| 9 | enumerate | 3/3 | 0 | 0 | 2.790 |

Both chart variants finish all 60 stage processes across development and holdout seeds. The ambient implementation times out on 33/60. Enumeration remains stronger in several binary-curve cases; the projected adapter reduces its overhead on the GF(8)-defined cases. No finite speedup ratio is assigned to a timed-out baseline. Zero-yield holdouts cannot meet a cost-per-verified-relation gate.

## Full-DLP practicality

All cases and all four seed/scalar combinations are included. Every row has 12 attempts. Cold time is the median among verified completions only; use the completion and timeout columns to avoid survivor bias. Rho is run once per distinct curve (case 1 also covers cases 2–3; case 4 covers case 5; case 7 covers cases 8–9).

| Case | Variant | Verified / attempts | Timeouts | Finished incomplete | Median cold ms, verified only |
|---|---|---|---|---|---|
| 0 | ambient_f4 | 12/12 | 0 | 0 | 14.126 |
| 0 | charts_linear | 12/12 | 0 | 0 | 0.666 |
| 0 | enumerate | 12/12 | 0 | 0 | 0.539 |
| 0 | pair_table | 12/12 | 0 | 0 | 0.602 |
| 1 | ambient_f4 | 0/12 | 12 | 0 | — |
| 1 | charts_linear | 0/12 | 0 | 12 | — |
| 1 | enumerate | 0/12 | 0 | 12 | — |
| 1 | pair_table | 0/12 | 0 | 12 | — |
| 2 | ambient_f4 | 12/12 | 0 | 0 | 148.014 |
| 2 | charts_linear | 12/12 | 0 | 0 | 1.211 |
| 2 | enumerate | 12/12 | 0 | 0 | 1.046 |
| 2 | pair_table | 12/12 | 0 | 0 | 1.169 |
| 3 | ambient_f4 | 0/12 | 12 | 0 | — |
| 3 | charts_linear | 0/12 | 0 | 12 | — |
| 3 | enumerate | 0/12 | 0 | 12 | — |
| 3 | pair_table | 0/12 | 0 | 12 | — |
| 4 | ambient_f4 | 0/12 | 12 | 0 | — |
| 4 | charts_linear | 12/12 | 0 | 0 | 6.249 |
| 4 | enumerate | 12/12 | 0 | 0 | 2.556 |
| 4 | pair_table | 12/12 | 0 | 0 | 1.665 |
| 5 | ambient_f4 | 9/12 | 3 | 0 | 2582.867 |
| 5 | charts_linear | 12/12 | 0 | 0 | 3.371 |
| 5 | enumerate | 12/12 | 0 | 0 | 2.129 |
| 5 | pair_table | 12/12 | 0 | 0 | 1.970 |
| 6 | ambient_f4 | 0/12 | 12 | 0 | — |
| 6 | charts_linear | 12/12 | 0 | 0 | 4.605 |
| 6 | enumerate | 12/12 | 0 | 0 | 3.867 |
| 6 | pair_table | 12/12 | 0 | 0 | 2.776 |
| 7 | ambient_f4 | 0/12 | 12 | 0 | — |
| 7 | charts_linear | 12/12 | 0 | 0 | 3.690 |
| 7 | enumerate | 12/12 | 0 | 0 | 5.188 |
| 7 | pair_table | 12/12 | 0 | 0 | 3.173 |
| 8 | ambient_f4 | 0/12 | 12 | 0 | — |
| 8 | charts_linear | 12/12 | 0 | 0 | 4.529 |
| 8 | enumerate | 12/12 | 0 | 0 | 5.898 |
| 8 | pair_table | 12/12 | 0 | 0 | 3.843 |
| 9 | ambient_f4 | 0/12 | 12 | 0 | — |
| 9 | charts_linear | 12/12 | 0 | 0 | 6.684 |
| 9 | enumerate | 12/12 | 0 | 0 | 12.466 |
| 9 | pair_table | 12/12 | 0 | 0 | 6.477 |
| 0 | rho | 12/12 | 0 | 0 | 0.215 |
| 1 | rho | 12/12 | 0 | 0 | 0.301 |
| 4 | rho | 12/12 | 0 | 0 | 0.414 |
| 6 | rho | 12/12 | 0 | 0 | 0.431 |
| 7 | rho | 12/12 | 0 | 0 | 0.597 |

## Fresh paired timing ratios

Candidate is projected charts throughout. Ratios below one favor the candidate. Only matched verified DLP completions or complete eight-target stage workloads enter a pair; missing pairs block an unqualified end-to-end claim. These are descriptive paired-bootstrap intervals on a shared host, not calibrated attack-cost ratios. Repeated seeds are not independent algorithm instances.

| Workload | Case | Reference | Pairs | Candidate/reference time | 95% interval |
|---|---|---|---|---|---|
| stage | 0 | ambient_f4 | 3 | 0.0108 | [0.0091, 0.0127] |
| stage | 0 | enumerate | 3 | 1.5465 | [1.2943, 1.7995] |
| dlp | 0 | ambient_f4 | 6 | 0.0540 | [0.0432, 0.0687] |
| dlp | 0 | enumerate | 6 | 1.2704 | [1.1052, 1.5326] |
| dlp | 0 | pair_table | 6 | 1.1635 | [0.9786, 1.4969] |
| dlp | 0 | rho | 6 | 3.1930 | [2.7191, 3.7082] |
| stage | 1 | ambient_f4 | 3 | 0.0010 | [0.0010, 0.0010] |
| stage | 1 | enumerate | 3 | 6.6327 | [6.2592, 7.0614] |
| dlp | 1 | ambient_f4 | 0 | — | — |
| dlp | 1 | enumerate | 0 | — | — |
| dlp | 1 | pair_table | 0 | — | — |
| dlp | 1 | rho | 0 | — | — |
| stage | 2 | ambient_f4 | 3 | 0.0012 | [0.0012, 0.0013] |
| stage | 2 | enumerate | 3 | 1.3096 | [1.2647, 1.3641] |
| dlp | 2 | ambient_f4 | 6 | 0.0111 | [0.0105, 0.0117] |
| dlp | 2 | enumerate | 6 | 1.0914 | [0.9567, 1.1892] |
| dlp | 2 | pair_table | 6 | 0.9854 | [0.9433, 1.0297] |
| dlp | 2 | rho | 6 | 3.6637 | [3.2518, 4.0815] |
| stage | 3 | ambient_f4 | 3 | 0.0001 | [0.0001, 0.0001] |
| stage | 3 | enumerate | 3 | 0.8598 | [0.7990, 0.9042] |
| dlp | 3 | ambient_f4 | 0 | — | — |
| dlp | 3 | enumerate | 0 | — | — |
| dlp | 3 | pair_table | 0 | — | — |
| dlp | 3 | rho | 0 | — | — |
| stage | 4 | ambient_f4 | 0 | — | — |
| stage | 4 | enumerate | 3 | 3.6699 | [3.6103, 3.7498] |
| dlp | 4 | ambient_f4 | 0 | — | — |
| dlp | 4 | enumerate | 6 | 2.5348 | [1.8624, 3.3074] |
| dlp | 4 | pair_table | 6 | 3.7484 | [2.9884, 4.7460] |
| dlp | 4 | rho | 6 | 17.3262 | [9.7763, 28.9917] |
| stage | 5 | ambient_f4 | 0 | — | — |
| stage | 5 | enumerate | 3 | 2.8081 | [2.6274, 3.0560] |
| dlp | 5 | ambient_f4 | 6 | 0.0014 | [0.0012, 0.0016] |
| dlp | 5 | enumerate | 6 | 1.6833 | [1.4703, 1.8991] |
| dlp | 5 | pair_table | 6 | 1.5687 | [1.2756, 1.9060] |
| dlp | 5 | rho | 6 | 8.2000 | [6.6406, 10.4788] |
| stage | 6 | ambient_f4 | 0 | — | — |
| stage | 6 | enumerate | 3 | 1.1986 | [1.1660, 1.2489] |
| dlp | 6 | ambient_f4 | 0 | — | — |
| dlp | 6 | enumerate | 6 | 1.2067 | [1.1769, 1.2372] |
| dlp | 6 | pair_table | 6 | 1.9722 | [1.6724, 2.3002] |
| dlp | 6 | rho | 6 | 13.8208 | [10.5124, 18.3350] |
| stage | 7 | ambient_f4 | 0 | — | — |
| stage | 7 | enumerate | 3 | 0.6611 | [0.6513, 0.6789] |
| dlp | 7 | ambient_f4 | 0 | — | — |
| dlp | 7 | enumerate | 6 | 0.8258 | [0.8072, 0.8439] |
| dlp | 7 | pair_table | 6 | 1.1374 | [1.1044, 1.1623] |
| dlp | 7 | rho | 6 | 3.4762 | [2.2098, 4.5436] |
| stage | 8 | ambient_f4 | 0 | — | — |
| stage | 8 | enumerate | 3 | 0.5897 | [0.5360, 0.6204] |
| dlp | 8 | ambient_f4 | 0 | — | — |
| dlp | 8 | enumerate | 6 | 0.8399 | [0.7621, 0.9784] |
| dlp | 8 | pair_table | 6 | 1.1077 | [0.9233, 1.3580] |
| dlp | 8 | rho | 6 | 3.9843 | [2.9163, 4.8600] |
| stage | 9 | ambient_f4 | 0 | — | — |
| stage | 9 | enumerate | 3 | 0.5592 | [0.5359, 0.5717] |
| dlp | 9 | ambient_f4 | 0 | — | — |
| dlp | 9 | enumerate | 6 | 0.5196 | [0.5017, 0.5374] |
| dlp | 9 | pair_table | 6 | 1.0888 | [1.0507, 1.1356] |
| dlp | 9 | rho | 6 | 12.0836 | [7.7098, 17.3932] |

## Degree 131: the mixed-product obstruction

The field polynomial `x^131+x^13+x^2+x+1` passes the full irreducibility check. The following exact field-algebra diagnostics use seed 17; all offset lists and seed-937 controls are in `structure.json`. No proper subfield of dimension 8, 16, 32 or 44 exists at prime degree 131.

| ell | Family | Self-product rank | Mixed minimum | Mixed maximum | Pair-weighted rank |
|---|---|---|---|---|---|
| 8 | power_span | 15 | 22 | 64 | 62.197 |
| 8 | frobenius_span | 36 | 43 | 64 | 62.727 |
| 8 | random | 36 | 64 | 64 | 63.576 |
| 16 | power_span | 31 | 46 | 131 | 127.364 |
| 16 | frobenius_span | 129 | 131 | 131 | 130.970 |
| 16 | random | 131 | 131 | 131 | 131.000 |
| 32 | power_span | 63 | 94 | 131 | 129.409 |
| 32 | frobenius_span | 131 | 131 | 131 | 131.000 |
| 32 | random | 131 | 131 | 131 | 131.000 |
| 44 | power_span | 87 | 130 | 131 | 130.318 |
| 44 | frobenius_span | 131 | 131 | 131 | 131.000 |
| 44 | random | 131 | 131 | 131 | 131.000 |

At ell=44 the power seed has rank 87 only on its self-product. Relative offsets 1 and 130 have rank 130; every other nonzero offset has rank 131. The complete cover has 8,646 unordered component pairs and weighted rank 130.318. Nearby offsets already explain the loss algebraically: products with the squared seed span 130 consecutive powers, and offset two contains a full field basis. This rejects the naive assumption that the low self-product dimension survives Frobenius composition. It does not rule out every other seed construction or fast algorithm.

## Decision and remaining boundary

The frozen 20% stage gate is met on the positive-yield fresh cases with complete ambient baselines (cases 0 and 2): paired cold cost per witness is 0.01080 and 0.001213 times ambient F4. This is an encoding-stage improvement. Zero-column cases 1 and 3 are rejected as candidate factor bases, and ambient timeouts block finite matched ratios for cases 4–9.

The stronger reference gate is not met. No fresh full-DLP chart/pair-table interval establishes a chart win, and the fresh chart/rho time ratios range from 3.19 to 17.33 across the eight useful-base cases. Rho verifies all 60 attempts. These runtime ratios remain secondary; calibrated attack-cost ratios are still null.

Keep the complete-cover adapter opt-in. Prefer common-subfield scalar components as the concrete composite-degree experiment; admission still depends on projected columns, relation yield and total verified work. Mixing GF(8) and GF(32) abscissa spaces would fill all 15 dimensions in their mixed product and lose the shared-subfield benefit.

At prime degree 131, a distinct future experiment could deliberately collect only selected component combinations and price the lost relation yield; that would be an incomplete relation oracle and must report unknown for the unsearched cover. The present implementation does not make that substitution. Higher arity, degree-131 solving, chart-plan Redis serialization, calibrated operation conversions and attack-exponent claims remain outside this measured implementation.

Derivations and primary references are in [README.md](README.md). The frozen thresholds and negative controls are in [contract.json](contract.json). Raw data, source/binary hashes and every paired comparison are in [results/run-001](results/run-001).
