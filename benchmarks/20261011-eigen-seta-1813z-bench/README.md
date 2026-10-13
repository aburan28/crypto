# Eigen-Séta historical KEM and field benchmark

[Source study](../../research/eigen_seta_rust_20261011) and [run/result record](../../research/eigen_seta_rust_20261011/reports/run-20261011T1813Z-bench/report.pdf) remain the source of truth. The PDF and SVG files below are byte-identical browsing copies added during integration on 2026-10-12. PNGs render the corresponding exact SVG at 150 dpi; no experiment or native result generation was run. The [copy manifest](copies.json) records source paths and hashes.

L1/L3/L5 decapsulation 13.424–14.997 / 51.834–53.940 / 120.180–123.953 s; 15/15 trials correct

Host and isolation: Intel Xeon Processor @ 2.80 GHz, 4 logical cores, Linux x86_64, shared cloud VM; no CPU pinning or frequency control. Original host/isolation wording and measurement qualifications are preserved in the copied PDFs and [verbatim source report](SOURCE_REPORT.md).

Historical source-bound records retained during PR #1680 reconciliation. The archived/current aggregate reports retain their respective claims and source bindings; current correctness checks do not revalidate historical side-channel or security claims. See research/eigen_seta_rust_20261011/INTEGRATION_RECONCILIATION.md.

## Visual copies

- [field_ops.svg](visuals/research/eigen_seta_rust_20261011/reports/run-20261011T1813Z-bench/field_ops.svg) — `research/eigen_seta_rust_20261011/reports/run-20261011T1813Z-bench/field_ops.svg`
- [kem_times.svg](visuals/research/eigen_seta_rust_20261011/reports/run-20261011T1813Z-bench/kem_times.svg) — `research/eigen_seta_rust_20261011/reports/run-20261011T1813Z-bench/kem_times.svg`
- [report.pdf](visuals/research/eigen_seta_rust_20261011/reports/run-20261011T1813Z-bench/report.pdf) — `research/eigen_seta_rust_20261011/reports/run-20261011T1813Z-bench/report.pdf`
- [report.pdf](visuals/research/eigen_seta_rust_20261011/report.pdf) — `research/eigen_seta_rust_20261011/report.pdf`
