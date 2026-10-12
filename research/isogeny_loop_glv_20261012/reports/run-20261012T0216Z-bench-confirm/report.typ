#set page(paper: "a4", margin: 1.8cm)
#set text(font: "DejaVu Sans", size: 9pt)
#set par(justify: false)
#set table(stroke: 0.4pt + rgb("#c3c2b7"), inset: 4pt)
#show table: set text(7.6pt)
#text(13pt, weight: "bold")[Run report: run-20261012T0216Z-bench-confirm]

*What ran.* `taskset -c 3 ./target/release/examples/gls_loop_bench --d -619 --pbits 127 --seed 20261012 --pairs 256 --rounds 30 --w 5 --json research/isogeny_loop_glv_20261012/bench_run2.json`. This is the planned confirmation of run-20261011T2333Z-bench (`bench_run1`). Source: branch `claude/isogeny-loop-glv-sweeper` at the round-3 commits, after helper moves and a clippy-only rename. rustc 1.97.0, release.

*Host.* A different machine from run 1: the container restarted onto an Intel Xeon \@ 2.80 GHz with 4 CPUs. Run 1 had a Xeon \@ 2.10 GHz. Load average 1.23 at start and 1.21 at end (`../../bench_run2_env.txt`); another agent was starting builds. Not isolated.

*Arms, unit, checks.* As in run 1. All correctness checks passed, and the sink value 2f8b53ca9ff88734 equals run 1's. Runs requested / completed: 1 / 1.

#table(columns: 6,
  [*arm*], [*median ns*], [*min ns*], [*wNAF / arm*], [*GLS-2 / arm*], [*run 1 GLS-2 / arm*],
  [wNAF], [211,895], [204,256], [1.000], [0.542], [0.572],
  [GLS-2], [114,823], [109,972], [1.845], [1.000], [1.000],
  [GLV-GLS-4 affine loop], [106,268], [101,386], [1.994], [1.081], [1.118],
  [GLV-GLS-4 projective loop], [105,273], [99,253], [2.013], [1.091], [1.147],
  [wNAF (A/A)], [214,312], [205,225], [0.989], [—], [—],
  [decompose4 (stage)], [5,306], [5,042], [—], [—], [—],
  [α′(P) affine (stage)], [9,556], [9,209], [—], [—], [—],
  [α′(P) projective (stage)], [7,574], [7,407], [—], [—], [—],
)
Ratios are derived from the recorded medians.

#figure(image("../../figures/bench_arms.svg", width: 100%), caption: [Per-round values for runs 1 and 2 at D = −619, 30 rounds × 256 scalars each. Run 2 is on a different host, so compare its absolute times only within the panel.])

*Observations only.* Interpretation is in `../../README.md` §3.
