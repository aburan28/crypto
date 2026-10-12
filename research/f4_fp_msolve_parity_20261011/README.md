# `f4_fp` against msolve and Singular (2026-10-12)

Results report: [`REPORT.pdf`](REPORT.pdf) (source `REPORT.typ`; graphs `fp_vs_msolve.svg`, `ratios.svg` with PNG renders).
Run report: [`reports/run-20261011T2359Z/report.pdf`](reports/run-20261011T2359Z/report.pdf), raw records in `reports/run-20261011T2359Z/raw/`.
Reproduce with `run.sh` (needs release builds of `examples/gb_cas_export.rs` with `f4_fp.rs` from `main` and from this branch, plus msolve and Singular); graphs with `plot.py <run dir>`.

Gebauer–Möller pair criteria, a matrix-step final interreduction, packed monomials, interned symbolic
preprocessing with a proper hasher, and a parallel symbolic pass. Every reduced basis is unchanged, and
each is confirmed by Singular or msolve.

| case | main ÷ new (1 thr) | msolve ÷ new (1 thr) | msolve ÷ new (4 thr) |
|:--|--:|--:|--:|
| quad_n9_p65521 | 29.8 | 1.22 | 0.83 |
| quad_n10_p65521 | 54.1 | 1.00 | 0.88 |
| quad_n11_p65521 | not run | 0.98 | 1.02 |
| cubic_n6_p31 | 47.0 | 0.96 | 1.11 |

Shared host, isolation L0, 3 interleaved rounds. Level with msolve is not "fastest": Magma was not run.
