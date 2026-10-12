# Koblitz m=3 symmetric-S4 solver comparison

[Source study](../../research/koblitz-isogeny) and [run/result record](../../research/koblitz-isogeny/reports/run-20261011T1955Z/report.pdf) remain the source of truth. The PDF and SVG files below are byte-identical browsing copies added during integration on 2026-10-12. PNGs render the corresponding exact SVG at 150 dpi; no experiment or native result generation was run. The [copy manifest](copies.json) records source paths and hashes.

n=17 per-call e3 0.000116 / 0.00141 / 0.0216 s at l=3/4/5; symmetric l=5 timed out (>300 s for four probes)

Host and isolation: 4-vCPU shared x86_64 cloud container, no core isolation; wall-clock values indicative only. Original host/isolation wording and measurement qualifications are preserved in the copied PDFs and [verbatim source report](SOURCE_REPORT.md).

Stage/presentation diagnostics on synthetic public curves; no end-to-end sub-rho claim.

## Visual copies

- [full_solves.svg](visuals/research/koblitz-isogeny/reports/run-20261011T1955Z/full_solves.svg) — `research/koblitz-isogeny/reports/run-20261011T1955Z/full_solves.svg`
- [per_call_n17.svg](visuals/research/koblitz-isogeny/reports/run-20261011T1955Z/per_call_n17.svg) — `research/koblitz-isogeny/reports/run-20261011T1955Z/per_call_n17.svg`
- [report.pdf](visuals/research/koblitz-isogeny/reports/run-20261011T1955Z/report.pdf) — `research/koblitz-isogeny/reports/run-20261011T1955Z/report.pdf`
