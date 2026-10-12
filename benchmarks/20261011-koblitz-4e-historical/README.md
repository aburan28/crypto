# Historical closed-form Koblitz 4E cost comparison

[Source study](../../research/isogeny_fibers_20261010/round4/16_2_fast_4e) and [run/result record](../../research/isogeny_fibers_20261010/round4/16_2_fast_4e) remain the source of truth. The PDF and SVG files below are byte-identical browsing copies added during integration on 2026-10-12. PNGs render the corresponding exact SVG at 150 dpi; no experiment or native result generation was run. The [copy manifest](copies.json) records source paths and hashes.

n=131/233 table check 0.124/0.122µs versus τ-NAF160/98µs (order-of-magnitude only); counted multiplications are primary

Host and isolation: x86_64,4 cores shared with about 6 other agents; single thread; no isolated timing. Original host/isolation wording and measurement qualifications are preserved in the copied PDFs and [verbatim source report](SOURCE_REPORT.md).

Antecedent cited by PR#1697; this benchmark uses the research example field backend and excludes on-curve checking. It does not validate timing of the integrated library implementation.

## Visual copies

- [costs.svg](visuals/research/isogeny_fibers_20261010/round4/16_2_fast_4e/costs.svg) — `research/isogeny_fibers_20261010/round4/16_2_fast_4e/costs.svg`
- [report.pdf](visuals/research/isogeny_fibers_20261010/round4/16_2_fast_4e/report.pdf) — `research/isogeny_fibers_20261010/round4/16_2_fast_4e/report.pdf`
