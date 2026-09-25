# Fixed quadratic coefficients accelerate complete generic search

The new representation preserves the complete search and reduces its measured cold cost. Primary holdouts pass all nine preregistered comparisons against the fastest of seven retained methods, with paired median ratios of **2.40–2.51x**. The unchanged confirmation gives **2.10–2.55x** paired medians but passes only seven of nine comparisons. Its two 16-variable lower confidence bounds fall below 2.0. The universal two-run gate is therefore **REJECTED**.

Every 20- and 24-variable comparison passes in both runs. That is a bounded observation on the frozen fixtures, not a new acceptance rule or a population-wide guarantee. All 6,144 recorded complete-solve observations finish with verified results.

## What changes

The solver compiles the original quadratic coefficients into immutable neighbor masks. Assigning a variable changes only the active-variable mask, residual affine coefficients and exact quadratic edge counts. XOR handles cancellations, edge counts expose degree drops, and exact coefficient/support comparison removes duplicate equations in the same order as the retained implementation. No monomial list is rebuilt at each node. The constructor is charged to each cold solve.

The contract and update equations are in [README.md](README.md), and the implementation is in [quadratic.rs](quadratic.rs). This restricted degree-at-most-two state is separate from the earlier general support-envelope product schedule. It needs no multiplier enumeration; a degree drop becomes an affine propagation event. General higher-degree inputs use the retained fallback. The earlier support-envelope timing failures are not replaced or reclassified by this result.

All search-only implementations produce identical outcomes/models, logical counters, generator ordering and trace digests. In particular, search-node counts and logical specialization input-term counts do not improve. This is **engineering**: a faster representation of the same work, with no reduction of the abstract search tree or evidence of a changed exponent.

## Cold complete costs at 24 variables

The table copies the confirmation evidence. Units are milliseconds per complete solve plus result validation, pooling two holdout seeds and eight rotated repetitions. The displayed ratio compares pooled medians against ordered specialization. Acceptance instead uses the pointwise fastest of all seven controls and paired intervals. Every row has correctness PASS and class engineering.

| Method | Planted (ms) | Cross-planted (ms) | Unplanted (ms) | Ordered / method, planted |
|---|---:|---:|---:|---:|
| Search-only | 83.948353 | 50.172584 | 172.670625 | 0.881 |
| Degree-3 flat kernel | 324.505584 | 110.399500 | 562.174501 | 0.228 |
| Degree-3 sparse bucket | 1039.781000 | 340.670584 | 2074.648312 | 0.071 |
| Degree-3 hybrid kernel | 432.203792 | 156.305584 | 920.137104 | 0.171 |
| Selective degree-2 flat | 98.751625 | 57.306479 | 201.832708 | 0.749 |
| Selective degree-2 one-word | 86.087291 | 49.644729 | 175.701604 | 0.859 |
| Ordered-specialization search | 73.976250 | 43.469708 | 151.862917 | 1.000 |
| Fixed-quadratic state | 30.557917 | 17.770312 | 62.343125 | 2.421 |

Primary-run fixed-quadratic times were 22.657812, 22.433333 and 25.485771 ms for the three families. The primary and confirmation have different holdout inputs, so their absolute times are not a before/after optimization comparison. Both full tables remain in [run_01/RESULT.md](run_01/RESULT.md) and [run_02/RESULT.md](run_02/RESULT.md).

## All preregistered comparisons

Ratios are fastest retained control / fixed-quadratic state. Intervals bootstrap the paired repeated timings on the two fixed fixtures in each size/family group. They do not estimate uncertainty over the wider distribution of Boolean systems. The implementation, controls and threshold were unchanged between runs.

| Variables | Family | Primary ratio [95% interval] | Confirmation ratio [95% interval] | Primary / confirmation gate |
|---:|---|---|---|---|
| 16 | planted | 2.458 [2.428, 2.519] | 2.439 [1.950, 2.609] | PASS / REJECTED |
| 16 | cross_planted | 2.457 [2.368, 2.539] | 2.102 [1.994, 2.237] | PASS / REJECTED |
| 16 | unplanted | 2.507 [2.472, 2.540] | 2.553 [2.517, 2.611] | PASS / PASS |
| 20 | planted | 2.423 [2.369, 2.459] | 2.435 [2.423, 2.479] | PASS / PASS |
| 20 | cross_planted | 2.472 [2.426, 2.501] | 2.441 [2.373, 2.501] | PASS / PASS |
| 20 | unplanted | 2.458 [2.437, 2.479] | 2.450 [2.417, 2.523] | PASS / PASS |
| 24 | planted | 2.412 [2.405, 2.431] | 2.409 [2.399, 2.426] | PASS / PASS |
| 24 | cross_planted | 2.422 [2.411, 2.447] | 2.448 [2.392, 2.464] | PASS / PASS |
| 24 | unplanted | 2.401 [2.367, 2.556] | 2.430 [2.412, 2.451] | PASS / PASS |

The failed 16-variable confirmation groups are retained. One planted fixture visits only 48 search nodes and the candidate takes about 54–60 microseconds in its eight repetitions. The cross-planted comparisons often select the one-word inference control as fastest. These observations identify short-solve and competing-policy costs to investigate; they do not establish that setup cost or timing noise alone caused either failure. No third confirmation was selected to replace these failures.

## Verification and evidence

The two grids contain **96 run cells**, **72 distinct generated systems** and **6,144 timed solver observations**. The 24 discovery fixtures recur in the second run. Each run also includes an untimed completed search reference per fixture. Primary outcomes are 38 SAT / 10 UNSAT fixtures per arm; confirmation outcomes are 35 SAT / 13 UNSAT. No measured result is censored.

SAT models are independently evaluated on the original equations. UNSAT outputs must agree with the completed retained search reference; this is a separate reference execution, not an external proof checker or unaffiliated reproduction. Small-system tests additionally compare with exhaustive Boolean enumeration. Trace hashes are diagnostics, not collision-resistant certificates.

Ten Rust tests cover exhaustive bounded solving, specialization, degree-two word kernels, trace equality, cancellations, duplicates, a second specialization, cubic fallback and UNKNOWN limits. Nineteen evidence tests replay both frozen manifests/reports, check unchanged confirmation source and fresh holdouts, reject altered models/traces/counters/profiles, reject selecting a weaker timing reference, and retain null completion cost under synthetic censoring. The new CI workflow runs the algebra and evidence checks without launching performance campaigns.

The maximum whole-worker peak RSS is 23,117,824 bytes, including fixtures, reference work and all eight arms. It is not a candidate-specific allocation measure. The same host ran both grids. Compilation, state setup, recursive solving, context destruction and result validation are inside cold arm timing; fixture generation, independent status-reference preparation and record formatting remain outside it and inside process receipts.

[RUN_LEDGER.json](RUN_LEDGER.json) pins both immutable runs and [SUMMARY.json](SUMMARY.json). The canonical scoreboard carries the same latest figures while preserving the predecessor studies. Full production-solver cost, calibrated operation ratio, complete index-calculus cost and rho ratio remain null. These bounded generic Boolean results do not establish a cryptanalytic crossover.

The next bounded experiment should isolate setup, trace enumeration and short-search costs, and compare any change against this new fastest representation. It must retain both failed confirmation groups as regressions and add fresh holdouts. Higher-degree coefficient updates require a new contract; this quadratic result does not validate them.
