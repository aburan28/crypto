# 2026-10-16 — native GLV-chain benchmark: Bandersnatch, Tom-384, CP6-782 G1

Copies for browsing; the source of truth is
[`research/endo_native_bench_20261016`](../../research/endo_native_bench_20261016/README.md)
(protocol, counting and timing records, frozen inputs, figure source
`figures/plot.py`, report sources). The figures are drawn only from the
committed JSON records; the PDFs and graphs here are byte-identical copies.

- [`REPORT.pdf`](REPORT.pdf) — results report
  (source `research/endo_native_bench_20261016/REPORT.typ`)
- [`run-report.pdf`](run-report.pdf) — run report
  (source `research/endo_native_bench_20261016/reports/run-20261016T1759Z/report.typ`)
- `fig_counted_vs_model` — counted vs model ratio per curve, and every total
  arm's counted deviation from the model;
  `fig_timed` — declared GLV arm vs the reference: counted, median and min
  ratio with the A/A floor, and per-round wall times;
  `fig_grid` — width × table grid of the declared chain (counted M_eq and
  median µs); `fig_phi` — φ(P) stage diagnostics, counted M_eq vs median µs
  (SVG and 150-dpi PNG each).
