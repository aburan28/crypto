# Eigen-Séta historical Toy16 instruction-count and timing probe

[Source study](../../research/eigen_seta_rust_20261011) and [run/result record](../../research/eigen_seta_rust_20261011/reports/run-20261011T1827Z-sidechannel/report.pdf) remain the source of truth. The PDF and SVG files below are byte-identical browsing copies added during integration on 2026-10-12. PNGs render the corresponding exact SVG at 150 dpi; no experiment or native result generation was run. The [copy manifest](copies.json) records source paths and hashes.

Non-allocator secret-phase instructions 1,116,186,003 in all 12 runs; dudect |t| ≤ 0.84

Host and isolation: Intel Xeon Processor @ 2.80 GHz, 4 logical cores, Linux x86_64, shared cloud VM; shared, one-key wall timing diagnostic. Original host/isolation wording and measurement qualifications are preserved in the copied PDFs and [verbatim source report](SOURCE_REPORT.md).

Historical source-bound records retained during PR #1680 reconciliation. The archived/current aggregate reports retain their respective claims and source bindings; current correctness checks do not revalidate historical side-channel or security claims. See research/eigen_seta_rust_20261011/INTEGRATION_RECONCILIATION.md.

## Visual copies

- [ct_probe.svg](visuals/research/eigen_seta_rust_20261011/reports/run-20261011T1827Z-sidechannel/ct_probe.svg) — `research/eigen_seta_rust_20261011/reports/run-20261011T1827Z-sidechannel/ct_probe.svg`
- [report.pdf](visuals/research/eigen_seta_rust_20261011/reports/run-20261011T1827Z-sidechannel/report.pdf) — `research/eigen_seta_rust_20261011/reports/run-20261011T1827Z-sidechannel/report.pdf`
- [report.pdf](visuals/research/eigen_seta_rust_20261011/report.pdf) — `research/eigen_seta_rust_20261011/report.pdf`
