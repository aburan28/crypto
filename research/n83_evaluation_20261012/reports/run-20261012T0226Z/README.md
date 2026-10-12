# Degree-83 retained IC / strong-rho evaluation

Completed diagnostic run, UTC start 2026-10-12 02:26. All 22 children completed verified; standalone arithmetic replay verified 22 answers and all six IC target witnesses. Two warm-ups are excluded from the five matched pairs and five A/A control pairs. No performance row is promoted.

The native paired geometric mean of IC/strong prepared-online cost is **1.382330**, 95% interval **[1.035774, 1.882300]**. Fresh-process wall from the retained base is **218.034569**, interval **[172.358384, 285.573727]**. These are reciprocals of the frozen strong/IC estimator, supplied by the native exporter; they are not ratios of arithmetic mean times. M4 has a lower IC online interval, and M5 is near a tie. The overall point-specific geometric mean gives greater IC cost.

Every pair repeats one fixed public point on one host. These intervals describe conditional seed/repeat variation, not independently sampled targets or hosts. Fresh-process cold starts with the retained base file, omits original representative selection/construction, and does not flush OS page cache. Both producers have timer residuals. Primary online values use complete intervals; outer process wall includes the frozen runner's polling/observer effects.

The report corrects the old 8.1813× packed-reference interpretation and wrong-modulus `cdcc5432` alias. The correct exact-field slug is `icv1-f2m83-t6151469093347-9eba6cd5`; polynomial, subgroup, G and Q remain the primary binding. The historical scalar/base/witness arithmetic remains valid. Full rank/log transcript, calibrated common operation-equivalent S, exact sealed ecbench L2 integration and a second-host admission are missing.

Files:

- [PDF](n83_strong_reference_evaluation.pdf) and [editable Markdown](report_source.md).
- [Structured source](report_content.json) and [native data copies](data/).
- [Recorded per-pair values](figures/paired_runtime.svg) and [native intervals](figures/paired_intervals.svg).
- [Three spaces](figures/three_spaces.svg) and [IC/rho coefficient transport](figures/coefficient_transport.svg).
- [Layout QA](preview/layout_qa.json) and [all-page contact sheet](preview/contact_sheet.jpg).

Rebuild with installed ReportLab, PyMuPDF, Pillow, Matplotlib and Inkscape:

```bash
python figures.py
python build_report.py
```

Python selects/formats native recorded values and lays out figures/PDF pages only; it computes no estimator, interval, arithmetic replay or benchmark result. The editable diagram SVGs describe existing source invariants without new algebraic discovery or timing claims. Final PDF inspection covers every page, graph and long identifier.

Authoritative executed binding: [records/frozen_panel_inputs.sha256](../../records/frozen_panel_inputs.sha256), corroborated by each start/execution receipt. The earlier preparation snapshot is not an executed identity. Source baseline: `a8194d616ed8c04489ec46177d369f1b1670360a`. Durable native records, source audit, validation logs and replay are in the [study](../../README.md).

Validation: final default-feature release library check **3947 passed, 0 failed, 112 ignored**; runner **7/7**, exporter **4/4**, IC example **10/10** and strong example **2/2**. Checks and independent replay occur separately from timing.

The study's canonical graph review distinguishes the a=0 / 81-bit capacity anchor from this a=1 retained diagnostic. This report is an accounting/reference correction; it makes no admitted frontier or ledger promotion. See the study README for graph/source checks and retained evidence limitations.
