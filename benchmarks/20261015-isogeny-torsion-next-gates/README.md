# 2026-10-15 isogeny and torsion next gates

Browse [results.pdf](results.pdf) for the short cross-gate report. The
[study report](../../research/isogeny_torsion_fibers_20261012/RESULTS_20261015.md)
links the frozen protocols, native sources, public inputs, raw results,
replay receipts and exact limitations. Every PDF, SVG and PNG in this
directory is a byte-identical copy of its study counterpart; `cmp` checked
the copies after generation. The source directories remain authoritative.

| Mirror stem | Source run report |
| --- | --- |
| `chi-exact` | [Exhaustive χ equation check](../../research/isogeny_torsion_fibers_20261012/chi_polynomial_20261015/reports/run-20261015T1635Z/REPORT.md) |
| `chi-solve` | [χ paired relation solve](../../research/isogeny_torsion_fibers_20261012/chi_polynomial_20261015/reports/run-20261015T1651Z/REPORT.md) |
| `different-j-selector` | [Different-j selector](../../research/isogeny_torsion_fibers_20261012/different_j_20261015/reports/run-20261015T1648Z/REPORT.md) |
| `different-j-solve` | [Paired different-j relation solve](../../research/isogeny_torsion_fibers_20261012/different_j_20261015/reports/run-20261015T1713Z/REPORT.md) |
| `different-j-fibers` | [Exhaustive rational-fiber audit](../../research/isogeny_torsion_fibers_20261012/different_j_20261015/reports/run-20261015T1738Z/REPORT.md) |
| `pairing-initial` | [Preserved initial batch attempt](../../research/isogeny_torsion_fibers_20261012/pairing_batch_20261015/reports/run-20261015T1700Z/REPORT.md) |
| `pairing-v3-v4` | [Phase and replay ledger](../../research/isogeny_torsion_fibers_20261012/pairing_batch_20261015/reports/run-20261015T1725Z/REPORT.md) |
| `pairing-v5` | [Extension-field rho control](../../research/isogeny_torsion_fibers_20261012/pairing_batch_20261015/reports/run-20261015T1729Z/REPORT.md) |
| `pairing-v6` | [Direct prime-field rho control](../../research/isogeny_torsion_fibers_20261012/pairing_batch_20261015/reports/run-20261015T1740Z/REPORT.md) |
| `pairing-v6b` | [Same-input repeat after release test exit](../../research/isogeny_torsion_fibers_20261012/pairing_batch_20261015/reports/run-20261015T1752Z/REPORT.md) |
| `pairing-v7` | [Direct source-group BSGS reuse control](../../research/isogeny_torsion_fibers_20261012/pairing_batch_20261015/reports/run-20261015T1805Z/REPORT.md) |

The benchmark index row is prepared but cannot be inserted while Conductor
task T-15 holds `benchmarks/README.md` exclusively. It must be added before
this entry is considered fully indexed; the mirror does not override that
required index.
