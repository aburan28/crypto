# Ecbench root-finding decomposition assessment

This directory holds browsing copies of the [study](../../research/ecbench_rootfinding_assessment_20261012/) visual record. The study directory is the source of truth for the frozen specs, sealed sessions, comparisons, audit receipts, source assessment, and figure source. The measurements use small public-target curves; each of the six method arms has 24 verified measured runs and all local executions passed full replay.

| Copy | Source |
| --- | --- |
| [Run report](run-report.pdf) | [`run_report.pdf`](../../research/ecbench_rootfinding_assessment_20261012/run_report.pdf) |
| [Results report](results.pdf) | [`report.pdf`](../../research/ecbench_rootfinding_assessment_20261012/report.pdf) |
| [Counted-cost graph](mean_s.svg) and [150 dpi PNG](mean_s.png) | [`mean_s.svg`](../../research/ecbench_rootfinding_assessment_20261012/mean_s.svg), rendered by `render_visuals.mjs` from sealed-record summaries |
| [Pipeline diagram](pipeline.svg) and [150 dpi PNG](pipeline.png) | [`pipeline.svg`](../../research/ecbench_rootfinding_assessment_20261012/pipeline.svg), rendered by `render_visuals.mjs` |

The PNGs are Typst renders of the matching SVGs at `--ppi 150` using the two `*_png.typ` wrappers in the study directory; `sips` sets their physical-resolution metadata to 150 dpi. The chart marks counted lower bounds; all wall grades are L0. The Koblitz pair-table/subtraction counted-$S$ quotient is 0.057047 (95% cluster interval 0.041159–0.092069), while the prime quotient is 0.018861 (0.018547–0.019242) with native work unpriced. The matched rho arms have lower counted $S$ on both curves; cross-method complete-cost and wall-time claims remain outside this measurement.
