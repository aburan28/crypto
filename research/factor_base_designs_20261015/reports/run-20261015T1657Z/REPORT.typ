#set document(title: "Factor-base designs: slope run on p1", author: "fb_designs_slope")
#set page(paper: "a4", margin: (x: 1.8cm, y: 1.8cm), numbering: "1")
#set text(size: 9.5pt)
#set par(justify: true)
#show raw: set text(size: 8pt)

= Factor-base designs: D against b on p#sub[1] = 786433 (run-20261015T1657Z)

*Kind:* solver-axis stage diagnostic (AGENTS.md §2). Observations only: nothing
here prices a pipeline, reports S or shows a speedup. The protocol is
#raw("research/factor_base_designs_20261015/PREREGISTRATION.md") §1–§3 with
Amendment 1 (registered before the first data row).

== Run

- *Command:* #raw("target/release/examples/fb_designs_slope run --arms A0,N,A1z,A1g,A4,A11 --m 2,3 --out-dir research/factor_base_designs_20261015/runs")
  (defaults: t = 2 … 18, 2 planted targets per cell, F4 budget 300 s, hard kill at 420 s, #raw("ulimit -v 13000000") per child process).
- *Commit measured:* #raw("a698ba0399836c56cba11a622737d6c42f418e11") (recorded in every row; clean for #raw("examples/"), #raw("src/") and Cargo files).
- *Host:* "vm", Intel Xeon @ 2.80 GHz, 4 cores, 15 GB RAM, a Claude Code cloud container. *Isolation:* none beyond one system at a time. F4 uses every core, and the instrument checks ran before the study, not alongside it. The times below are wall-clock practicality notes, not the metric.
- *Window:* 2026-10-15, 16:57Z to about 19:42Z (UTC).
- *Instrument:* #raw("f4_fp::f4"), grevlex, degree bound n + d#sub[max] + 6, D = #raw("solving_degree_max"), with Gebauer–Möller criteria on (the current default; Amendment 1 item 1).
- *Seeds:* cell seed #raw("0x46424453 ⊕ arm≪40 ⊕ m≪32 ⊕ t≪16") draws the set and the curve. Target k is planted from #raw("cell ⊕ (k+1)≪8").
- *Records:* #raw("runs/<arm>-m<m>.jsonl") and #raw(".log") beside them, #raw("runs/analysis.txt") (native #raw("analyze")), #raw("runs/instrument-check.txt") and #raw("runs/instrument-check-gm0.txt").

== Instrument check (§5 stop rule)

A0 was rebuilt with the pilot's seed rule and draw order (#raw("repro")): 20 planted systems, m = 2 with t ≤ 7 and m = 3 with t ≤ 4. On every system, D, x#sub[R], the curve, the tower, the unknown and equation counts, #raw("inconsistent") and #raw("planted_ok") equal every committed pilot row in #raw("research/pkm_tower_pilot_20260924/runs/").

#raw("degree_reached") is lower than the pilot's on 16 of the 20 systems (the certification tail). With #raw("F4_FP_GM=0"), all 20 systems match on every field, #raw("degree_reached") included. The stop rule does not fire.

This study's own A0 series (FBDS seeds, so other towers, curves and targets) has the same D as the pilot's planted rows at every (m, t) both have: m = 2, t = 2–7 and m = 3, t = 2–3. At m = 3, t = 4, every pilot row timed out.

== Graph

#figure(image("d_vs_b.svg", width: 100%),
  caption: [Measured D (#raw("solving_degree_max")) against b = log#sub[2]|V| per summand, per arm, m = 2 (left) and m = 3 (right). Each completed system is plotted, both targets per cell (n = 2 per cell; the two coincide in every cell). Open markers: systems stopped by the 300 s F4 budget, at their partial D (a lower bound). Crosses: memory aborts, no D. Run run-20261015T1657Z, rows from #raw("runs/*.jsonl"), plotted by #raw("plot.py").])

== σ and the preregistered verdict

From #raw("runs/analysis.txt"). σ is the least-squares slope of the cell D (the mean of the cell's 2 systems) on b, over completed cells.

#table(columns: 5, align: (left, center, center, center, left), stroke: 0.4pt + gray,
  [*arm*], [*m*], [*completed cells (t)*], [*σ*], [*verdict (§1 thresholds)*],
  [A0], [2], [8 (2–9)], [0.1905], [candidate],
  [A0], [3], [3 (2–4)], [0.5000], [no middle ground],
  [N], [2], [8 (2–9)], [0.1905], [candidate],
  [N], [3], [3 (2–4)], [0.5000], [no middle ground],
  [A1z], [2], [5 (2–6)], [0.0000], [candidate],
  [A1z], [3], [3 (2–4)], [0.0000], [candidate],
  [A1g], [2], [5 (2–6)], [0.0000], [candidate],
  [A1g], [3], [3 (2–4)], [0.0000], [candidate],
  [A4], [2], [6 (2–7)], [0.2286], [candidate],
  [A4], [3], [3 (2–4)], [0.5000], [no middle ground],
  [A11], [2], [8 (2–9)], [0.1786], [candidate],
  [A11], [3], [3 (2–4)], [0.0000], [candidate],
)

Observations, stated as recorded:
- The slope rule labels the incumbent A0 itself "candidate" at m = 2, where D runs 4, 4, 4, 4, 5, 5, 5, 5 over t = 2–9. At m = 3, A0's three completed cells read 5, 6, 6.
- N has the same D as A0 at every completed cell, at both arities.
- A1z and A1g have D = 8 at every completed cell. 8 is also the highest input degree of their systems: x is quadratic in the unknowns, so S#sub[3] has degree 8. σ = 0 there means D stays at the input degree over b = 2–6 (m = 2) and b = 2–4 (m = 3). At equal t, A1z − A0 is +2 to +4, outside the predicted ±1. A1g − A1z is 0 at every shared cell, within the predicted ±1.
- A4 had no collisions at any t: |V| = 2#super[t], so b = t. At equal t, A4 − A0 is +1 in every shared cell but one (m = 3, t = 3: 0). Against A0 at ⌈t/2⌉, A4 is +1 to +2.
- A11's D is at or below every other arm's at every shared t: 3–4 at m = 2, where it ties A0 and N at t = 5, and 4 at m = 3. Its highest input degree is 2 (m = 2) and 3 (m = 3), with 3 + 3t and 5 + 5t unknowns.
- §5 says an arm reading "candidate" on p#sub[1] is replicated on p#sub[2] = 2013265921 before it is described as more than a single-prime reading. No p#sub[2] run was made here. At m = 3, every series has exactly 3 completed cells, the minimum the thresholds accept.

== Every system

D, time and status per target (target 0 / target 1). Time is F4 wall seconds, or process wall seconds for an abort. Status is ok, timeout (F4 budget of 300 s) or aborted (memory allocation failed under #raw("ulimit -v 13000000"), SIGABRT). Every series stopped at its first limit cell (Amendment 1 item 5). No planted-solution violation occurred.

#set text(size: 8pt)
#table(columns: 9, align: center, stroke: 0.3pt + gray,
  [*arm*], [*m*], [*t*], [*|V|*], [*vars*], [*eqs*], [*D*], [*time (s)*], [*status*],
  [A0], [2], [2], [4], [4], [5], [4 / 4], [0.0 / 0.0], [ok, ok],
  [A0], [2], [3], [8], [6], [7], [4 / 4], [0.0 / 0.0], [ok, ok],
  [A0], [2], [4], [16], [8], [9], [4 / 4], [0.0 / 0.0], [ok, ok],
  [A0], [2], [5], [32], [10], [11], [4 / 4], [0.1 / 0.1], [ok, ok],
  [A0], [2], [6], [64], [12], [13], [5 / 5], [0.3 / 0.3], [ok, ok],
  [A0], [2], [7], [128], [14], [15], [5 / 5], [2.0 / 2.0], [ok, ok],
  [A0], [2], [8], [256], [16], [17], [5 / 5], [12.9 / 12.9], [ok, ok],
  [A0], [2], [9], [512], [18], [19], [5 / 5], [92.2 / 92.0], [ok, ok],
  [A0], [2], [10], [1024], [20], [21], [5 / 5], [302.0 / 302.3], [timeout, timeout],
  [A0], [3], [2], [4], [7], [8], [5 / 5], [0.0 / 0.0], [ok, ok],
  [A0], [3], [3], [8], [10], [11], [6 / 6], [0.6 / 0.6], [ok, ok],
  [A0], [3], [4], [16], [13], [14], [6 / 6], [16.6 / 14.9], [ok, ok],
  [A0], [3], [5], [32], [16], [17], [6 / 6], [304.2 / 301.6], [timeout, timeout],
  [N], [2], [2], [4], [4], [5], [4 / 4], [0.0 / 0.0], [ok, ok],
  [N], [2], [3], [8], [6], [7], [4 / 4], [0.0 / 0.0], [ok, ok],
  [N], [2], [4], [16], [8], [9], [4 / 4], [0.0 / 0.0], [ok, ok],
  [N], [2], [5], [32], [10], [11], [4 / 4], [0.1 / 0.0], [ok, ok],
  [N], [2], [6], [64], [12], [13], [5 / 5], [0.3 / 0.4], [ok, ok],
  [N], [2], [7], [128], [14], [15], [5 / 5], [2.2 / 2.0], [ok, ok],
  [N], [2], [8], [256], [16], [17], [5 / 5], [14.8 / 15.0], [ok, ok],
  [N], [2], [9], [512], [18], [19], [5 / 5], [97.2 / 99.4], [ok, ok],
  [N], [2], [10], [1024], [20], [21], [5 / 5], [302.5 / 302.6], [timeout, timeout],
  [N], [3], [2], [4], [7], [8], [5 / 5], [0.0 / 0.0], [ok, ok],
  [N], [3], [3], [8], [10], [11], [6 / 6], [0.9 / 1.1], [ok, ok],
  [N], [3], [4], [16], [13], [14], [6 / 6], [24.7 / 24.5], [ok, ok],
  [N], [3], [5], [32], [16], [17], [6 / 6], [301.3 / 302.2], [timeout, timeout],
  [A1z], [2], [2], [4], [6], [7], [8 / 8], [0.0 / 0.0], [ok, ok],
  [A1z], [2], [3], [8], [10], [11], [8 / 8], [0.0 / 0.0], [ok, ok],
  [A1z], [2], [4], [16], [14], [15], [8 / 8], [0.1 / 0.1], [ok, ok],
  [A1z], [2], [5], [32], [18], [19], [8 / 8], [2.1 / 2.4], [ok, ok],
  [A1z], [2], [6], [64], [22], [23], [8 / 8], [110.2 / 106.0], [ok, ok],
  [A1z], [2], [7], [128], [26], [27], [– / –], [231.2 / 210.9], [aborted, aborted],
  [A1z], [3], [2], [4], [10], [11], [8 / 8], [0.0 / 0.0], [ok, ok],
  [A1z], [3], [3], [8], [16], [17], [8 / 8], [1.1 / 1.0], [ok, ok],
  [A1z], [3], [4], [16], [22], [23], [8 / 8], [237.5 / 247.9], [ok, ok],
  [A1z], [3], [5], [32], [28], [29], [8 / 8], [320.9 / 322.0], [timeout, timeout],
  [A1g], [2], [2], [4], [6], [7], [8 / 8], [0.0 / 0.0], [ok, ok],
  [A1g], [2], [3], [8], [10], [11], [8 / 8], [0.0 / 0.0], [ok, ok],
  [A1g], [2], [4], [16], [14], [15], [8 / 8], [0.1 / 0.2], [ok, ok],
  [A1g], [2], [5], [32], [18], [19], [8 / 8], [1.7 / 1.4], [ok, ok],
  [A1g], [2], [6], [64], [22], [23], [8 / 8], [69.8 / 69.0], [ok, ok],
  [A1g], [2], [7], [128], [26], [27], [8 / 8], [301.4 / 300.7], [timeout, timeout],
  [A1g], [3], [2], [4], [10], [11], [8 / 8], [0.0 / 0.0], [ok, ok],
  [A1g], [3], [3], [8], [16], [17], [8 / 8], [1.2 / 1.1], [ok, ok],
  [A1g], [3], [4], [16], [22], [23], [8 / 8], [281.3 / 266.4], [ok, ok],
  [A1g], [3], [5], [32], [28], [29], [– / –], [310.8 / 311.6], [aborted, aborted],
  [A4], [2], [2], [4], [6], [7], [5 / 5], [0.0 / 0.0], [ok, ok],
  [A4], [2], [3], [8], [8], [9], [5 / 5], [0.0 / 0.0], [ok, ok],
  [A4], [2], [4], [16], [10], [11], [5 / 5], [0.0 / 0.0], [ok, ok],
  [A4], [2], [5], [32], [12], [13], [5 / 5], [0.5 / 0.4], [ok, ok],
  [A4], [2], [6], [64], [14], [15], [6 / 6], [2.6 / 2.3], [ok, ok],
  [A4], [2], [7], [128], [16], [17], [6 / 6], [83.7 / 79.8], [ok, ok],
  [A4], [2], [8], [256], [18], [19], [6 / 6], [301.1 / 301.8], [timeout, timeout],
  [A4], [3], [2], [4], [10], [11], [6 / 6], [0.0 / 0.0], [ok, ok],
  [A4], [3], [3], [8], [13], [14], [6 / 6], [0.3 / 0.3], [ok, ok],
  [A4], [3], [4], [16], [16], [17], [7 / 7], [64.0 / 64.3], [ok, ok],
  [A4], [3], [5], [32], [19], [20], [8 / 8], [302.8 / 301.9], [timeout, timeout],
  [A11], [2], [2], [4], [9], [10], [3 / 3], [0.0 / 0.0], [ok, ok],
  [A11], [2], [3], [8], [12], [13], [3 / 3], [0.0 / 0.0], [ok, ok],
  [A11], [2], [4], [16], [15], [16], [3 / 3], [0.0 / 0.0], [ok, ok],
  [A11], [2], [5], [32], [18], [19], [4 / 4], [0.1 / 0.1], [ok, ok],
  [A11], [2], [6], [64], [21], [22], [4 / 4], [1.0 / 0.9], [ok, ok],
  [A11], [2], [7], [128], [24], [25], [4 / 4], [3.7 / 3.6], [ok, ok],
  [A11], [2], [8], [256], [27], [28], [4 / 4], [15.6 / 14.8], [ok, ok],
  [A11], [2], [9], [512], [30], [31], [4 / 4], [83.0 / 78.8], [ok, ok],
  [A11], [2], [10], [1024], [33], [34], [4 / 4], [301.0 / 301.0], [timeout, timeout],
  [A11], [3], [2], [4], [15], [16], [4 / 4], [0.0 / 0.0], [ok, ok],
  [A11], [3], [3], [8], [20], [21], [4 / 4], [0.1 / 0.1], [ok, ok],
  [A11], [3], [4], [16], [25], [26], [4 / 4], [1.8 / 1.8], [ok, ok],
  [A11], [3], [5], [32], [30], [31], [5 / 5], [300.8 / 300.8], [timeout, timeout],
)
#set text(size: 9.5pt)

== Requested and completed

12 series, 70 cells, 140 planted systems. 116 completed, 20 stopped at the F4 budget, 4 aborted for memory (A1z m = 2 t = 7, both targets; A1g m = 3 t = 5, both targets). The limit cells are missing data, and no D enters any fit from them.

== Deviations from the protocol as frozen

Amendment 1 records them all before the data: the instrument's Gebauer–Möller default, the "same cells" check done with the pilot's seed rule, the cell and target seed split, the cell-mean fit, the stop rule with a 420 s hard kill, N's chain at m = 3, and A11's sign convention and m = 3 form. Nothing else deviated. The replication on p#sub[2] that §5 requires for a "candidate" arm was not run in this lane.
