# Matched legacy degree-83 IC and optimized rho

Campaign14October2026; native verification and publication15October2026.
This directory contains browsing copies of the study's inspected PDFs and
SVG/150dpi PNG figures. The [study](../../research/n83_matched_20261014/README.md)
is authoritative; [native statistics](../../research/n83_matched_20261014/records/postprocess-20261015T0018Z/statistics.json)
record every number below. Ratios are **IC/rho; below1 favors IC**.

The standalone addendum has17 pages; the comprehensive report has87, including
the preceding70 preserved pages. `COPY_MANIFEST.sha256` binds these copies.

| Actual complete-process workload | IC mean seconds | Rho mean seconds | IC/rho | Conditional95% interval |
|---|---:|---:|---:|---|
| Eight cold singleton pairs | 1159.349882 | 3.071649 | 377.435691 | [311.812847,499.010210] |
| Three genuinely reused32-target batch pairs | 1232.975217 | 78.957400 | 15.615702 | [14.743543,17.614848] |

Optimized rho is faster at both measured complete-process boundaries. Online
target means remain secondary: IC/rho0.716385, interval[0.460389,1.013127]
and crossed target/cohort sensitivity[0.452168,1.046078], both including1.
The IC clock excludes input/Q-validation/output that rho includes. Three
preparation clusters on one guest do not establish host-general uncertainty.
The model's1685-target crossover is unobserved;814/20000 draws are nonfinite,
and the extended95% upper endpoint is no eventual IC advantage.

Exact curve: `icv1-f2m83-t6151469093347-9eba6cd5`, a=1, tail149,
degree83 field,53-bit prime subgroup r8569786107849059. IC K600/rank29;
rho Fold/reference-u64/256 lanes/DP6/escape-off. Both use CPU3, one OS thread,
12GiB address-space limit and2100s whole-process budget. Cold begins with
materialized K600; its original construction and OS-cache flushing are excluded.

All32 processes completed,218 scalars and104 IC witnesses passed independent
arithmetic replay, and all five A/A work pairs matched. Main rho work includes
4488 internal capped walks, which remain charged; there was no process timeout
or reported solve failure. The memory-limit-hit indicator is unavailable.
The separately metered external audit is outside both methods' process costs.

- [Standalone addendum](n83_matched_ic_rho_evaluation.pdf)
- [Updated comprehensive PDF](unified_curve_basis_report_matched_ic_rho.pdf)
- [Postprocessing revision history](../../research/n83_matched_20261014/POSTPROCESSING.md)
- [Native evidence and inspected report source](../../research/n83_matched_20261014/reports/run-20261014T1940Z/)

Canonical progress and leaderboard entries are diagnostics. Calibrated S,
ecbench L2, independent-host replay and full IC rank/log transcript admission
remain unset. These figures do not promote an admitted best/floor/ledger row,
claim a global optimum, or multiply historical gains across different panels.
