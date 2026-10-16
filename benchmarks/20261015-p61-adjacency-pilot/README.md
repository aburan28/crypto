# p61 selected-degree adjacency production calibration

This benchmark mirror records startup-inclusive generation and independent
replay for one-, 64-, and 4,096-vertex windows of the selected degree
`5/13/37/59` relation on class `p61-i0025` over `p = 2^61 - 1`.

The 4,096-vertex production window emitted 28,672 directed edges and 471,040
kernel coefficients in 191.551 seconds, then passed the distinct kernel-first
verifier in 296.931 seconds. It occupies 41,543,871 logical JSONL bytes and
8,552,556 gzip bytes. These are single observations on a shared, unisolated
macOS arm64 host while the p61 vertex census was active; they are descriptive
capacity measurements, not a comparative speed result.

The study source of truth is
[`research/p61_isogeny_class_20261014/reports/run-20261015T0436Z-adjacency-pilot/`](../../research/p61_isogeny_class_20261014/reports/run-20261015T0436Z-adjacency-pilot/).

- [PDF report](p61-adjacency-pilot-report.pdf)
- [Calibration graph SVG](adjacency-calibration.svg)
- [Calibration graph PNG](adjacency-calibration.png)
- [Evidence pipeline SVG](evidence-pipeline.svg)
- [Evidence pipeline PNG](evidence-pipeline.png)
