# Fresh matched degree-83 IC versus optimized rho

Status: all32 native processes completed;218 scalar answers,104 IC witnesses
and five identical A/A work pairs independently verified. Frozen native
analysis passed every required row. Inspected PDF/evidence publication is
being completed. User authorized execution on14October2026.
Worker/controller/analysis source freeze:
`c6e8daafb5a59d27e1be88a2f61f0bd6ee83ef03`.

The [prospective protocol](PROTOCOL.md) freezes the exact legacy a=1/tail149
curve,32 new point-only targets, three genuinely reused-setup IC/rho cohort
pairs, eight actual singleton cold pairs, and five rho A/A pairs. IC uses the
retained K600 and rank seed29; rho uses the previously selected Fold/reference
u64/256-lane/DP6/escape-off policy. Both are single-threaded under one native
CPU3 reservation,12GiB address-space limit and2100s full-process budget.

The primary costs include setup, solving, recovery and built-in verification.
Cold measurements start from the retained K600 file; historical factor-base
selection/construction and OS-cache flushing are outside this boundary.
Independent external arithmetic/file audit is separately metered. Producer
online boundaries differ and remain secondary. No old IC time is divided by
a newly measured rho time.

| User requirement | Current status | Evidence or remaining gate |
|---|---|---|
| Larger unseen matched panel | Fixtures independently verified | [Generation and native replay](verification/README.md) |
| Freeze both algorithms before testing | Verified before main timing | [Protocol](PROTOCOL.md), [build receipt](provenance/build_receipt.json), [frozen matrix](provenance/frozen_campaign_manifest.json) |
| Actual cold-start comparison | Verified complete | Eight actual singleton pairs, [native statistics](records/postprocess-20261015T0018Z/statistics.json) |
| Genuine setup reuse for both methods | Verified complete | Three32-target cohort pairs, each method preparing once per process |
| Complete cost, memory and failures | Verified complete | [All32 raw process directories](records/campaign-20261014T1940Z/runs/), no solver failure/timeout; internal work retained |
| Means, medians, slow cases and uncertainty | Verified native analysis | [Input](records/postprocess-20261015T0018Z/statistics_input.json), [statistics](records/postprocess-20261015T0018Z/statistics.json); three preparation clusters and nonfinite modeled mass retained |
| Independent correctness and provenance | Verified complete with explicit postprocessing corrections | [V3 PASS receipt](records/postprocess-20261015T0018Z/independent_replay_v3_receipt.json), [sealed extraction receipt](records/postprocess-20261015T0018Z/collector_extraction_receipt.json), [revision history](POSTPROCESSING.md) |
| Comprehensive PDF and reviewable publication | In progress after verified native results | New addendum preserving prior70pages; extend PR1828 |

All ratios below are **IC/rho; below1 favors IC**. Native recorded primary
means and paired percentile95% intervals are:

| Complete-process workload | IC mean seconds | Rho mean seconds | IC/rho | Conditional95% interval |
|---|---:|---:|---:|---|
| Eight actual cold singleton pairs | 1159.349882 | 3.071649 | 377.435691 | [311.812847,499.010210] |
| Three genuine32-target batch pairs | 1232.975217 | 78.957400 | 15.615702 | [14.743543,17.614848] |

Optimized rho is faster in both observed complete-process workloads. These
are local results for two frozen policies on one guest, not global algorithm
optimality. The batch interval has only three preparation clusters. Cold
begins with materialized K600 and does not include its historical construction.

Secondary producer-online target means are1.760364 s IC and2.457289 s rho:
ratio0.716385, target-block95% interval[0.460389,1.013127], crossed
target/cohort sensitivity[0.452168,1.046078]. Both include parity; clock
boundaries differ. The model's1685-target threshold is unobserved. Its814
nonfinite draws out of20000 (0.0407) are retained; the extended95% upper endpoint
is no eventual IC advantage, not the finite-only quantile.

The original v2 checker and frozen collector were both attempted and their
metadata round-trip failures retained. A separately versioned, narrowly scoped
v3 verifier and SHA-sealed collect-only extractor passed without changing any
measured code, fixture, raw record, protocol or statistical estimand. See
[postprocessing corrections](POSTPROCESSING.md). No solver was rerun.

Previous reports and all frozen studies remain historical evidence. This
campaign is a local diagnostic, with calibrated S, L2/physical-host admission,
independent-host claims and global parameter optimality unset. A complete
matrix and independent verification are required for a full-workload winner.

Conductor is unavailable; the pre-edit scope attempt and exit127 are retained
under `checks/conductor.*`. Separately, current PR-CI repairs install the Lean
workflow's missing ripgrep dependency and make curve-name diff checking honor
the existing frozen-archive exemption without modifying archived bytes.
