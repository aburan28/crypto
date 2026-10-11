# Inherited F4: eager rewrite of kept rows (2026-10-11)

Results report: [`REPORT.pdf`](REPORT.pdf) (source `REPORT.typ`, graph `totals.svg`).
Run report: [`reports/run-20261011T1657Z/report.pdf`](reports/run-20261011T1657Z/report.pdf),
raw probe output in `reports/run-20261011T1657Z/raw/`.

`f6_ic_probe`, one binary, `KIC_F4_EAGER=0` (lazy rows) against the default
(eager), three interleaved rounds, medians, isolation L0:

| cell | arm | lazy ms | eager ms | lazy / eager |
|:--|:--|--:|--:|--:|
| K₀/2²³ dim 8, 4 targets | inherited F4 | 1477.0 | 1058.4 | 1.40 |
| | F6-IC + index | 1034.8 | 639.5 | 1.62 |
| K₀/2³¹ dim 9, 1 target | inherited F4 | 1543.5 | 1041.3 | 1.48 |
| | F6-IC + index | 1443.6 | 946.1 | 1.53 |

Reductions, splits, gate counters and every verdict are identical between
modes. A constant-factor change to the Gröbner track; the |F|^(m−2)
enumeration and every IC-versus-rho comparison are untouched.
