# Eigen-Séta L3/L5 authenticated PKE timing samples

[Source study](../../research/eigen_seta_rust_20261011) and [run/result record](../../research/eigen_seta_rust_20261011/reports/run-20261011T1958Z-pke-l3-l5/report.pdf) remain the source of truth. The PDF and SVG files below are byte-identical browsing copies added during integration on 2026-10-12. PNGs render the corresponding exact SVG at 150 dpi; no experiment or native result generation was run. The [copy manifest](copies.json) records source paths and hashes.

L3/L5 Encrypt 1.1051 / 2.2990 s; Decrypt 47.138 / 104.223 s; one sample per size

Host and isolation: AMD EPYC 7763 shared cloud VM, Linux x86-64, 5 logical CPUs; no isolation or performance comparison; concurrent runs with other local validation. Original host/isolation wording and measurement qualifications are preserved in the copied PDFs and [verbatim source report](SOURCE_REPORT.md).

Historical source-bound records retained during PR #1680 reconciliation. The archived/current aggregate reports retain their respective claims and source bindings; current correctness checks do not revalidate historical side-channel or security claims. See research/eigen_seta_rust_20261011/INTEGRATION_RECONCILIATION.md.

## Visual copies

- [report.pdf](visuals/research/eigen_seta_rust_20261011/reports/run-20261011T1958Z-pke-l3-l5/report.pdf) — `research/eigen_seta_rust_20261011/reports/run-20261011T1958Z-pke-l3-l5/report.pdf`
- [timings.svg](visuals/research/eigen_seta_rust_20261011/reports/run-20261011T1958Z-pke-l3-l5/timings.svg) — `research/eigen_seta_rust_20261011/reports/run-20261011T1958Z-pke-l3-l5/timings.svg`
- [report.pdf](visuals/research/eigen_seta_rust_20261011/report.pdf) — `research/eigen_seta_rust_20261011/report.pdf`
- [ct_probe.svg](visuals/research/eigen_seta_rust_20261011/figures/ct_probe.svg) — `research/eigen_seta_rust_20261011/figures/ct_probe.svg`
- [ct_probe_verify.svg](visuals/research/eigen_seta_rust_20261011/figures/ct_probe_verify.svg) — `research/eigen_seta_rust_20261011/figures/ct_probe_verify.svg`
- [f_filter.svg](visuals/research/eigen_seta_rust_20261011/figures/f_filter.svg) — `research/eigen_seta_rust_20261011/figures/f_filter.svg`
- [field_ops.svg](visuals/research/eigen_seta_rust_20261011/figures/field_ops.svg) — `research/eigen_seta_rust_20261011/figures/field_ops.svg`
- [kem_times.svg](visuals/research/eigen_seta_rust_20261011/figures/kem_times.svg) — `research/eigen_seta_rust_20261011/figures/kem_times.svg`
- [pke_roundtrips.svg](visuals/research/eigen_seta_rust_20261011/figures/pke_roundtrips.svg) — `research/eigen_seta_rust_20261011/figures/pke_roundtrips.svg`

These aggregate SVGs accompany the maintained results PDF and include separately identified historical cohorts; their original host/source qualifications remain in the source report. The later PKE record does not revalidate historical side-channel claims.
