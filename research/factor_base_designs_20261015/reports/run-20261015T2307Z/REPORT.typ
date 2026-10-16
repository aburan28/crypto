#set document(title: "Factor-base designs: A11 replication on p2", author: "fb_designs_slope")
#set page(paper: "a4", margin: (x: 1.8cm, y: 1.8cm), numbering: "1")
#set text(size: 9.5pt)
#set par(justify: true)
#show raw: set text(size: 8pt)

= Factor-base designs: A11 replication on p#sub[2] = 2013265921 (run-20261015T2307Z)

*Kind:* solver-axis stage diagnostic (AGENTS.md §2). Observations only: nothing here prices a pipeline, reports S or shows a speedup. *Protocol:* PREREGISTRATION.md §5, with Amendment 2 F1 (registered on #raw("research/factor-base-designs-lean") in commit #raw("e0afa082d"), before any p#sub[2] row) and Amendment 1 for everything F1 does not change. This is F1 only. F2–F4 were not run here.

== Run

- *Command:* #raw("target/release/examples/fb_designs_slope run --p 2013265921 --arms A0,A11 --m 2,3 --out-dir research/factor_base_designs_20261015/runs/p2")
  (defaults: t from 2, at most 27 = v#sub[2]\(p#sub[2] − 1), 2 planted targets per cell, F4 budget 300 s, hard kill 420 s, #raw("ulimit -v 13000000") per child).
- *Commit measured:* #raw("d5f675260867f2d025e080e65d8d7027431bc845"), recorded in every row.
- *Host:* "vm", Intel Xeon @ 2.10 GHz, 4 cores, a Claude Code cloud container. This is a different container and CPU from the p#sub[1] run (2.80 GHz), so times are not comparable across the two runs. *Isolation:* one system at a time, no other work running.
- *Window:* 2026-10-15 23:07Z to 2026-10-16 about 00:08Z (UTC).
- *Seeds (Amendment 2):* #raw("0x46424453 ⊕ prime≪48 ⊕ arm≪40 ⊕ m≪32 ⊕ t≪16") with p#sub[2] = 1. Target k is planted from #raw("cell ⊕ (k+1)≪8"). With p#sub[1] = 0, the p#sub[1] seeds are unchanged: A0 and A11 at m = 3, t = 3 were rebuilt and match their committed rows (x#sub[R], D, max_cols).
- *Records (after reconciliation, see below):* every system this run measured is in the F1 rows of record #raw("runs/F1-{A0,A11}-m{2,3}-p2.jsonl"), except A11 m = 2, t = 10 and t = 11, which are in #raw("runs/F1-annex-p2/A11-m2-p2-t10-t11.jsonl"). The annex also holds this run's logs and the row-by-row cross-check (#raw("crosscheck.txt")). The rows as first committed are at #raw("runs/p2/") in commit #raw("cb750c828"), with #raw("analysis.txt") and #raw("instrument-check-p2.txt").

== F1 pre-check

p#sub[2] < 2#super[31], and #raw("f4_fp") reduces modulo p < 2#super[32]. Eight planted systems at t = 2 were built before any p#sub[2] row: A0 and A11, m = 2 and 3, both targets. Each passed the set, presentation and vanishing checks, and each kept its planted solution in F4's basis (#raw("planted_ok") true, not inconsistent). The check was recorded in #raw("runs/p2/instrument-check-p2.txt") (commit #raw("cb750c828")). The parent lane's own F1 pre-check is #raw("runs/f1-instrument-check-p2.txt").

== Graphs

#figure(image("d_vs_b.svg", width: 100%),
  caption: [Measured D against b on p#sub[2], A0 and A11, m = 2 (left) and m = 3 (right). Both targets of each cell are plotted (n = 2 per cell; they coincide). Open markers are systems stopped by the 300 s F4 budget, at their partial D. Rows from #raw("runs/F1-*-p2.jsonl") and the annex, per #raw("plot.py").])

#figure(image("width_vs_b.svg", width: 100%),
  caption: [Widest F4 matrix (#raw("max_cols"), log#sub[2] scale) against b, A0 and A11 on p#sub[1] (dashed, hollow) and p#sub[2] (solid). The points are offset ±0.08 in b so both primes stay visible. Completed systems only.])

== F1 decision rule (registered for m = 3 on p#sub[2])

"A11 replicates if A11's cell D ≤ A0's cell D − 1 and A11's #raw("max_cols") < A0's at every shared completed t, with at least 3 shared completed cells." The cell #raw("max_cols") here is the mean over the cell's two systems.

#table(columns: 6, align: center, stroke: 0.4pt + gray,
  [*t*], [*A11 D*], [*A0 D*], [*A11 max_cols*], [*A0 max_cols*], [*condition*],
  [2], [4], [5], [62], [368], [met],
  [3], [4], [6], [781], [3285], [met],
  [4], [4], [6], [3096], [12314], [met],
)

*Reading: replicates*, on 3 shared completed cells, the minimum the rule accepts. At t = 5 both arms' systems hit the F4 budget, so that cell is missing for both, with partial D of 5 (A11) and 6 (A0) when stopped.

σ on p#sub[2] (from #raw("analysis.txt")): A0 m = 2 is 0.1905 (8 cells), A0 m = 3 is 0.5000 (3), A11 m = 2 is 0.1500 (9), and A11 m = 3 is 0.0000 (3). The registered width-growth slope Δlog#sub[2] #raw("max_cols")/Δb (descriptive) is 1.261, 2.532, 1.215 and 2.821, in the same order.

== Observations

- *p#sub[1] and p#sub[2] agree system for system.* 44 planted systems completed on both primes at the same (arm, m, t, target): A0 and A11, m = 2 and 3. All 44 have the same D, and all 44 have the same #raw("max_cols"). The two primes share no curve, set or target. At these sizes, D and matrix width do not vary with the prime or the instance.
- *m = 2 (descriptive; the rule is not registered here).* A11's D is one below A0's at t = 2–4 and 6–9, and equal at t = 5. A11 is narrower at t = 2–4 and 8–9, and wider at t = 5–7. The widest of those is t = 5: 1431 vs 606 columns.
- *Further cells.* A11 at m = 2 completed one more cell on p#sub[2] (t = 10, D = 4, about 257 s per system) than on p#sub[1], where t = 10 timed out at 300 s on a different CPU. A0 at m = 2 stops at t = 10 on both primes.
- *The two targets of a cell are not always distinct systems.* At small t only 2–4 elements of V lie on the curve, so both targets can draw the same abscissae. Then x#sub[R] coincides, and for the tower arms the system is identical. On p#sub[2] this happened in A11 at t = 2 (m = 2 and 3). On p#sub[1], eight cells at t ≤ 4 across arms have it. D is unaffected, but those cells hold one distinct system, not two. The draw rule is the registered one and was not changed.

== Every system (p#sub[2])

Target 0 / target 1. Time is F4 wall seconds. Status is ok or timeout (F4 budget of 300 s). There were no memory aborts and no planted-solution violations.

#set text(size: 8pt)
#table(columns: 9, align: center, stroke: 0.3pt + gray,
  [*arm*], [*m*], [*t*], [*vars*], [*eqs*], [*D*], [*max_cols*], [*time (s)*], [*status*],
  [A0], [2], [2], [4], [5], [4 / 4], [37 / 37], [0.0 / 0.0], [ok, ok],
  [A0], [2], [3], [6], [7], [4 / 4], [127 / 127], [0.0 / 0.0], [ok, ok],
  [A0], [2], [4], [8], [9], [4 / 4], [280 / 280], [0.0 / 0.0], [ok, ok],
  [A0], [2], [5], [10], [11], [4 / 4], [606 / 606], [0.0 / 0.0], [ok, ok],
  [A0], [2], [6], [12], [13], [5 / 5], [2451 / 2451], [0.2 / 0.2], [ok, ok],
  [A0], [2], [7], [14], [15], [5 / 5], [5031 / 5031], [1.1 / 1.2], [ok, ok],
  [A0], [2], [8], [16], [17], [5 / 5], [9119 / 9119], [9.3 / 9.0], [ok, ok],
  [A0], [2], [9], [18], [19], [5 / 5], [14922 / 14922], [59.6 / 59.4], [ok, ok],
  [A0], [2], [10], [20], [21], [5 / 5], [61361 / 61361], [301.1 / 301.5], [timeout, timeout],
  [A0], [3], [2], [7], [8], [5 / 5], [368 / 368], [0.0 / 0.0], [ok, ok],
  [A0], [3], [3], [10], [11], [6 / 6], [3285 / 3285], [0.3 / 0.3], [ok, ok],
  [A0], [3], [4], [13], [14], [6 / 6], [12314 / 12314], [9.5 / 9.7], [ok, ok],
  [A0], [3], [5], [16], [17], [6 / 6], [97456 / 97456], [301.0 / 301.3], [timeout, timeout],
  [A11], [2], [2], [9], [10], [3 / 3], [31 / 31], [0.0 / 0.0], [ok, ok],
  [A11], [2], [3], [12], [13], [3 / 3], [78 / 78], [0.0 / 0.0], [ok, ok],
  [A11], [2], [4], [15], [16], [3 / 3], [208 / 208], [0.0 / 0.0], [ok, ok],
  [A11], [2], [5], [18], [19], [4 / 4], [1431 / 1431], [0.0 / 0.1], [ok, ok],
  [A11], [2], [6], [21], [22], [4 / 4], [2673 / 2673], [0.5 / 0.6], [ok, ok],
  [A11], [2], [7], [24], [25], [4 / 4], [5282 / 5282], [2.4 / 2.3], [ok, ok],
  [A11], [2], [8], [27], [28], [4 / 4], [8763 / 8763], [10.6 / 10.4], [ok, ok],
  [A11], [2], [9], [30], [31], [4 / 4], [13579 / 13579], [53.9 / 53.7], [ok, ok],
  [A11], [2], [10], [33], [34], [4 / 4], [22118 / 22118], [259.5 / 255.7], [ok, ok],
  [A11], [2], [11], [36], [37], [4 / 4], [184920 / 184920], [302.0 / 301.9], [timeout, timeout],
  [A11], [3], [2], [15], [16], [4 / 4], [62 / 62], [0.0 / 0.0], [ok, ok],
  [A11], [3], [3], [20], [21], [4 / 4], [781 / 781], [0.0 / 0.1], [ok, ok],
  [A11], [3], [4], [25], [26], [4 / 4], [3096 / 3096], [1.3 / 1.3], [ok, ok],
  [A11], [3], [5], [30], [31], [5 / 5], [40486 / 40486], [300.8 / 301.6], [timeout, timeout],
)
#set text(size: 9.5pt)

== Requested and completed

4 series, 27 cells, 54 planted systems: 46 completed and 8 stopped at the F4 budget (A0 m = 2 t = 10, A0 m = 3 t = 5, A11 m = 2 t = 11, A11 m = 3 t = 5, both targets each). The commit that landed these rows (#raw("c195dfa64")) states "50 systems, 42 completed" in its message. That was a miscount; these figures are correct.

== Reconciliation with the parent lane (2026-10-16)

The parent lane ran F1 on its own, earlier (commits #raw("f287673fe") and #raw("fa60aba73")), with the same seed rule. #raw("fb_designs_slope compare-rows") compared the two runs system by system (#raw("runs/F1-annex-p2/crosscheck.txt")). Of 52 shared systems, 50 are identical in construction and F4 output. The other 2, A11 m = 2 t = 10, have identical construction, but the parent's run timed out at 300 s and this run completed (D = 4 in both). This run also has A11 m = 2 t = 11 (timeouts).

The parent's files are the F1 rows of record. This run's duplicate rows were removed, and its extra rows were kept in the annex. Every number in this report is unchanged, and the graphs rebuild byte-identical from the new sources. The F1 verdict (m = 3: replicates) rests on m = 3 rows that are identical in both runs. The σ values above use this run's rows. With the rows of record, A11 m = 2 on p#sub[2] reads σ = 0.1786 on 8 cells, instead of 0.1500 on 9; both are "candidate".

== Deviations

None from Amendment 2 F1. Not run in this lane: F2 (the N11 null), F3 (the m = 3, t = 5 cell at 3600 s) and F4 (m = 4). The parent lane ran all three (#raw("runs/analysis-amendment2.txt")). F2 reads "shape, not structure": the A11-shaped null has A11's D in every completed cell. So the D difference this run replicates comes from the shape of A11's equations, not from the curve.
