# Registered m=3 refutation-degree law and native Macaulay solve diagnostics

[Source study](../../research/ic_semiregular_law_20261005/RESULTS.md) is the source of truth. These browsing copies were added during integration on 2026-10-12. PDFs and SVGs are copied byte-for-byte; PNGs render the same SVGs at 150 dpi.

68/68 units resolved; x4 meets its generic degree window on 32/32 draws, rr falls below on 8. Heavy degree-10 solves: 239k–247k rows ×199k columns, 42–51 min on 4 threads and about 90 min on 2; registered verdict neither

Host and isolation: Local 4-core, 16-GB container (two concurrent 2-thread lanes) plus three remote 4-core hosts; no isolation beyond one process per lane; wall times advisory.

Run started 2026-10-05T16:58Z and finished 2026-10-11T23:27Z, with restart and usage-limit interruptions retained in the source. Native macaulay_dense solves; historical Python orchestration/analysis remains frozen. The n=131 cost is a conditional model output, not a measured end-to-end cost.

Run reports:

- [run-20261005T1658Z](../../research/ic_semiregular_law_20261005/reports/run-20261005T1658Z/report.pdf)
