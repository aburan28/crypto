#set document(title: "Factor-base designs: Amendment 2 run", author: "fb_designs_slope")
#set page(paper: "a4", margin: (x: 1.8cm, y: 1.8cm), numbering: "1")
#set text(size: 9.5pt)
#set par(justify: true)
#show raw: set text(size: 8pt)

= Factor-base designs: the A11 follow-up F1–F4 (run-20261015T2052Z)

*Kind:* solver-axis stage diagnostic (AGENTS.md §2). Observations only: nothing
here prices a pipeline, reports S or shows a speedup. The protocol is
#raw("PREREGISTRATION.md") Amendment 2 (F1–F4), with Amendment 3 (registered
before any F1–F4 row). The motivation is #raw("reports/review-20261015T2030Z/REVIEW.pdf"), point 4.

== Run

- *Commands* (#raw("B=target/release/examples/fb_designs_slope"), #raw("R=research/factor_base_designs_20261015/runs"), run serially in this order): \
  #raw("$B check-p2") \
  #raw("$B regress --rows $R/A11-m3.jsonl --arm A11 --m 3 --t 2") \
  #raw("$B run --prime p2 --arms A0,A11 --m 2,3 --file 'F1-{arm}-m{m}-p2' --out-dir $R") \
  #raw("$B run --arms N11 --m 2,3 --file 'F2-{arm}-m{m}' --out-dir $R") \
  #raw("$B run --arms A0,A11 --m 4 --file 'F4-{arm}-m{m}' --out-dir $R") \
  #raw("$B run --arms A0,A11 --m 3 --t-min 5 --t-max 5 --budget 3600 --grace 300 --file 'F3-{arm}-m3-t5' --out-dir $R").
  Defaults: 2 planted targets per cell, F4 budget 300 s and hard kill 420 s (F3: 3600 s and 3900 s), #raw("ulimit -v 13000000") per child process, t rising from 2 until the first cell with a limit.
- *Binary:* built once from commit #raw("2ca6973a6"). F1 rows record #raw("2ca6973a6…"). F2, F3 and F4 rows record #raw("f287673fe…"), a later commit that added only F1 rows (no change under #raw("examples/") or #raw("src/")). The same binary ran every series.
- *Host:* "vm", Intel Xeon @ 2.80 GHz, 4 cores, 15 GB RAM, a Claude Code cloud container. *Isolation:* one system at a time; F4 uses every core. Times are wall-clock practicality notes, not the metric.
- *Window:* 2026-10-15, 20:52Z to 23:40Z (UTC).
- *Seeds:* #raw("0x46424453 ⊕ prime≪48 ⊕ arm≪40 ⊕ m≪32 ⊕ t≪16"), p#sub[1] = 0, p#sub[2] = 1, N11 = arm 6. Target k from #raw("cell ⊕ (k+1)≪8").
- *Records:* #raw("runs/F1-<arm>-m<m>-p2.jsonl"), #raw("runs/F2-N11-m<m>.jsonl"), #raw("runs/F3-<arm>-m3-t5.jsonl"), #raw("runs/F4-<arm>-m4.jsonl"), each with its #raw(".log"); #raw("runs/f1-instrument-check-p2.txt"); #raw("runs/regression-A11-m3-t2.txt"); #raw("runs/analysis-amendment2.txt") (native #raw("analyze --amendment2")). Every number below is copied from these files.

== Checks before the rows

*p#sub[2] instrument check: PASS.* p#sub[2] = 15·2#super[27] + 1 is prime, and 2#super[27] ∥ p#sub[2] − 1. Since p#sub[2] > 2#super[16], #raw("f4_fp") takes its wide path at p#sub[2], as at p#sub[1]: one Barrett reduction per update, with operands at most (p−1) + (p−1)#super[2] < 2#super[64]. The 32-bit and lazy-accumulation paths apply only to p ≤ 2#super[16]. A0 and A11 at t = 2, m = 2 and 3 (target 0) passed every check. Each system vanishes at its planted solution, F4 finishes consistent with the planted solution on its basis, no pair lies above the degree bound, every input reduces to 0 modulo the basis, and a second run gives the same D.

*Regression on p#sub[1]: PASS.* A11, m = 3, t = 2 was rebuilt with the extended harness. Both rows are byte-identical to the committed rows on D and 26 structural fields (seeds, set, curve, target, system shape, F4 trace). The first run's #raw("analysis.txt") is also reproduced byte for byte.

== Graphs

#figure(image("d_vs_b.svg", width: 100%),
  caption: [Measured D (#raw("solving_degree_max")) against b = log#sub[2]|V| per summand, for the F1 (p#sub[2]), F4 (m = 4) and F2 (N11, beside the committed A11 rows on p#sub[1]) series. n = 2 planted systems per cell. The two coincide in every completed cell. Open markers: systems stopped by the 300 s budget, at their partial D (a lower bound). Run run-20261015T2052Z, plotted from the rows by #raw("plot.py").])

#figure(image("cols_vs_b.svg", width: 100%),
  caption: [Largest F4 matrix width (#raw("max_cols"), log#sub[2] scale) against b, for the same series. Open markers: the width reached by budget-stopped systems (a lower bound). Run run-20261015T2052Z, #raw("plot.py").])

== Decision rules (Amendment 2), as #raw("analyze --amendment2") prints them

*F1, does A11 replicate on p#sub[2] at m = 3? → REPLICATES* (prediction held). There are 3 shared completed cells:

#table(columns: 5, align: center, stroke: 0.4pt + gray,
  [*t*], [*D A0*], [*D A11*], [*max_cols A0*], [*max_cols A11*],
  [2], [5], [4], [368], [62],
  [3], [6], [4], [3 285], [781],
  [4], [6], [4], [12 314], [3 096],
  [5], [timeout (300 s)], [timeout (300 s)], [≥ 97 456], [≥ 40 486],
)

Observed alongside the rule: on every completed system of F1, the p#sub[2] row has the same D, #raw("max_cols"), #raw("max_rows") and #raw("steps") as the p#sub[1] row at the same (arm, m, t, target). This holds for A0 and A11 at m = 2 and 3, although the curves, sets and targets differ (checked with #raw("jq") on the row files). The only difference is in the budget-stopped A11 m = 2, t = 10 rows, which stopped at different points of their traces. At these sizes the F4 trace did not depend on the prime or the coefficients drawn.

*F2, shape or structure? → SHAPE, NOT STRUCTURE* (prediction held). There are 11 completed N11 cells (m = 2: t = 2–9; m = 3: t = 2–4). Each is compared with the committed A11 cell, and N11 − A11 = 0 in D in all 11. Widths are equal at m = 2 in every cell. At m = 3 they differ: 62 vs 62, 745 vs 781 and 3 300 vs 3 096 (N11 vs A11). Both series stopped at the same cells (m = 2, t = 10; m = 3, t = 5).

*F3, the m = 3, t = 5 cell on p#sub[1], reported as recorded.*
- *A11:* both systems completed, with D = 5 and 5, #raw("max_cols") 40 486 and 40 486, in 564.7 s and 571.2 s. That width equals the width the 300 s run had reached when it stopped. With this cell as a fourth point, A11's m = 3 p#sub[1] fit becomes σ = 0.3000 over 4 cells, which reads *weak* (0.25 < σ < 0.5). Without it, σ was 0.0000 over 3 cells (*candidate*). The prediction "both A11 systems complete within 3600 s" held.
- *A0:* both systems aborted on a memory allocation failure under #raw("ulimit -v 13000000"), at 1542.3 s and 1548.6 s (SIGABRT). They recorded no D and no width. The cell is not completed, so A0's m = 3 fit stays σ = 0.5000 over 3 cells.

*F4, m = 4 on p#sub[1]: both arms are inconclusive* (fewer than 3 completed cells). Both series stopped at t = 4 on the 300 s budget.

#table(columns: 6, align: center, stroke: 0.4pt + gray,
  [*arm*], [*t = 2: D, width*], [*t = 3: D, width*], [*t = 4*], [*σ (2 cells)*], [*vars / eqs at t = 2, 3, 4*],
  [A0], [6, 1 911], [7, 21 619], [timeout, D ≥ 7, width ≥ 102 516], [1.0000], [10/11, 14/15, 18/19],
  [A11], [3, 34], [5, 2 159], [timeout, D ≥ 6, width ≥ 88 740], [2.0000], [21/22, 28/29, 35/36],
)

A11 at m = 4: n_eqs − n_vars = 1 on all 6 rows (#raw("norm_chain_excess")). lead(G) = 1 was checked symbolically in every built system (Amendment 3, item 3). The prediction "A0's σ ≥ 0.5" is *not decided* (2 cells). Over the two completed cells, A11's D is 3 below A0's at t = 2 and 2 below at t = 3, and its width is smaller at both.

== Secondary: Δlog#sub[2] max_cols / Δb (descriptive, not a verdict)

These are least-squares slopes over completed cells, from #raw("analysis-amendment2.txt").

#table(columns: 4, align: (left, center, center, center), stroke: 0.4pt + gray,
  [*series*], [*cells*], [*Δlog#sub[2] max_cols/Δb*], [*σ*],
  [A0 m=2, p#sub[1] / p#sub[2]], [8 / 8], [1.261 / 1.261], [0.1905 / 0.1905],
  [A11 m=2, p#sub[1] / p#sub[2]], [8 / 8], [1.314 / 1.314], [0.1786 / 0.1786],
  [N11 m=2, p#sub[1]], [8], [1.314], [0.1786],
  [A0 m=3, p#sub[1] / p#sub[2]], [3 / 3], [2.532 / 2.532], [0.5000 / 0.5000],
  [A11 m=3, p#sub[1] / p#sub[2]], [3 / 3], [2.821 / 2.821], [0.0000 / 0.0000],
  [A11 m=3, p#sub[1] with F3 (t = 2–5)], [4], [3.004], [0.3000],
  [N11 m=3, p#sub[1]], [3], [2.867], [0.0000],
  [A0 m=4, p#sub[1]], [2], [3.500], [1.0000],
  [A11 m=4, p#sub[1]], [2], [5.989], [2.0000],
)

At m = 3, A11's width per bit grows at least as fast as A0's (2.821 vs 2.532 over t = 2–4). Its widths are lower than A0's at every shared t. At m = 4, A11's width is lower at t = 2 and 3, and its per-bit growth over those two cells is larger (5.989 vs 3.500).

== Predictions scorecard (Amendment 2)

#table(columns: 3, align: (left, left, left), stroke: 0.4pt + gray,
  [*prediction*], [*observed*], [*held?*],
  [F1: A11 replicates], [replicates (3 shared cells; the p#sub[2] trace equals the p#sub[1] trace)], [yes],
  [F2: shape, not structure], [N11 = A11 in D in all 11 compared cells], [yes],
  [F3: both A11 systems complete within 3600 s (A0's may not)], [A11: both completed, D = 5, about 570 s; A0: both aborted on memory at about 1545 s], [yes (the A0 clause is not scored)],
  [F4: A0's σ ≥ 0.5], [σ = 1.0 on 2 cells: inconclusive], [not decided],
  [F4: A11, no prediction], [σ = 2.0 on 2 cells: inconclusive], [—],
)

== Deviations and limits

- *Amendment 3* was registered before any F1–F4 row. Item 3 corrects Amendment 2: at m = 4, G is monic (N = Q#super[2]\(X#super[3] + AX + B) − P#super[2] has leading term X#super[5]), not of leading coefficient 1 − p#sub[2]#super[2]. So the "vanishing lead" redraw cannot trigger. The 4 × 4 system for (q#sub[0], P) is replanted when it is singular. Items 1 and 2 fix N11's coefficient stream and how a cell's width is formed (the mean over its 2 systems).
- F3's A0 systems ended on the memory limit, not the time budget. Amendment 2 left the memory limit unchanged, so the cell is recorded as aborted, with no partial D (the process died before printing one).
- The F3 scorecard line in #raw("analyze") was corrected after the F2 and F4 rows were recorded, so that a missing cell reads "not decided" instead of "did not hold". This was analysis code only, and the run binary was not rebuilt while F3 ran.
- The run had no isolation beyond one system at a time, on a shared 4-core container. Time-limit stops depend on the host.
