#set document(title: "Factor-base designs: Amendment 4 run", author: "fb_designs_slope")
#set page(paper: "a4", margin: (x: 1.8cm, y: 1.8cm), numbering: "1")
#set text(size: 9.5pt)
#set par(justify: true)
#show raw: set text(size: 8pt)

= Factor-base designs: cost score, third m = 4 cell, small-prime instrument, A11 ablation (run-20261016T0616Z)

*Kind:* solver-axis stage diagnostic (AGENTS.md §2). Observations only. Nothing here prices a pipeline, reports S or shows a speedup. The protocol is #raw("PREREGISTRATION.md") Amendment 4 (C1–C4), registered before any of its rows. The analysis is #raw("runs/analysis-amendment4.txt") (native #raw("analyze --amendment4")), and every number below is copied from it or from the row files.

== Run

- *Commands* (#raw("B=target/release/examples/fb_designs_slope"), #raw("R=research/factor_base_designs_20261015/runs"), run serially): \
  #raw("$B replay $R --out C1-rescore --groups p1,F1,F2,F3,F4") \
  #raw("$B run --arms H1,H2 --m 2,3 --file 'C4-{arm}-m{m}' --out-dir $R") \
  #raw("$B run --arms A0,A11 --m 4 --t-min 4 --t-max 4 --budget 3600 --grace 300 --file 'C2-{arm}-m4-t4' --out-dir $R") \
  #raw("$B replay $R --out C3-trace-p3 --groups p1,F2,F3,F4,C2,C4 --only-prime p1 --prime p3") \
  The C3 extension ran only if every compared system matched. It did not run (#raw("runs/C3-extension-not-run.txt")).
- *Binary:* built once from #raw("025746fd5"), the Amendment 4 harness commit. Later recorded commits (#raw("236e3bd58"), #raw("86958f03e")) only added rows, Lean files and README lines, with no change under #raw("examples/") or #raw("src/"). Before the run, the extended harness reproduced the committed regression, #raw("analysis.txt") and #raw("analysis-amendment2.txt") byte for byte. #raw("field_ops") was the same on two runs of the same system.
- *Host:* "vm", Intel Xeon @ 2.80 GHz, 4 cores, 15 GB RAM, a Claude Code cloud container. One system at a time, with no further isolation. *Window:* 2026-10-16, 06:16Z to 09:54Z (UTC).
- *Limits:* 300 s F4 budget and 420 s hard kill (C2: 3600 s and 3900 s), #raw("ulimit -v 13000000"). Seeds as in Amendments 1–3, plus prime index p#sub[3] = 2 and arms H1 = 7, H2 = 8.
- *Records:* #raw("runs/C1-rescore.jsonl"), #raw("runs/C4-H{1,2}-m{2,3}.jsonl"), #raw("runs/C2-{A0,A11}-m4-t4.jsonl"), #raw("runs/C3-trace-p3.jsonl"), each with its #raw(".log"); #raw("runs/analysis-amendment4.txt").
- *Lean* (outside the measurement): #raw("quartic_norm_factor") in #raw("lean/FactorBaseDesigns/Encodings.lean") proves Q#super[2]\(X#super[3] + AX + B) − P#super[2] = (X − x#sub[R])·G with G a monic quartic. That is Amendment 3's correction. #raw("lake build"), the axiom audit (standard axioms only) and #raw("leanchecker --fresh") pass. The logs are in #raw("proof_logs/").

== C1, cost score κ = Δlog#sub[2] field_ops / Δb against m − 1

The replay rebuilt and re-measured all 192 committed systems in completed cells (first run, F1, F2, F3, F4). *192 of 192 match their committed rows on D, #raw("max_cols"), #raw("max_rows") and #raw("steps").* #raw("field_ops") counts F#sub[p] multiplications in F4's row reduction only.

#figure(image("ops_vs_b.svg", width: 100%),
  caption: [#raw("field_ops") (cell mean of 2 systems, log#sub[2] scale) against b, completed cells on p#sub[1], m = 2, 3, 4. The dashed line has the break-even slope m − 1 bits per bit, drawn through the highest first point. It is a reference, not a fit. Rows from #raw("runs/C1-rescore.jsonl"), #raw("C2-*") and #raw("C4-*"). Run run-20261016T0616Z, #raw("plot.py").])

#table(columns: 6, align: (left, center, center, center, center, left), stroke: 0.4pt + gray,
  [*series*], [*cells*], [*κ*], [*m − 1*], [*σ*], [*C1 verdict*],
  [A0 m=2 p#sub[1]], [8], [3.593], [1], [0.1905], [unfavourable],
  [N m=2 p#sub[1]], [8], [3.596], [1], [0.1905], [unfavourable],
  [A1z m=2 p#sub[1]], [5], [6.050], [1], [0.0000], [unfavourable],
  [A1g m=2 p#sub[1]], [5], [5.809], [1], [0.0000], [unfavourable],
  [A4 m=2 p#sub[1]], [6], [4.902], [1], [0.2286], [unfavourable],
  [A11 m=2 p#sub[1]], [8], [3.995], [1], [0.1786], [unfavourable],
  [N11 m=2 p#sub[1]], [8], [4.014], [1], [0.1786], [unfavourable],
  [H1 m=2 p#sub[1]], [8], [3.535], [1], [0.2619], [unfavourable],
  [H2 m=2 p#sub[1]], [8], [3.923], [1], [0.0000], [unfavourable],
  [A0 m=3 p#sub[1]], [3], [7.468], [2], [0.5000], [unfavourable],
  [N m=3 p#sub[1]], [3], [7.772], [2], [0.5000], [unfavourable],
  [A1z m=3 p#sub[1]], [3], [8.905], [2], [0.0000], [unfavourable],
  [A1g m=3 p#sub[1]], [3], [9.189], [2], [0.0000], [unfavourable],
  [A4 m=3 p#sub[1]], [3], [8.028], [2], [0.5000], [unfavourable],
  [A11 m=3 p#sub[1] (with F3, t = 2–5)], [4], [8.769], [2], [0.3000], [unfavourable],
  [N11 m=3 p#sub[1]], [3], [8.582], [2], [0.0000], [unfavourable],
  [H1 m=3 p#sub[1]], [3], [7.660], [2], [0.5000], [unfavourable],
  [H2 m=3 p#sub[1]], [3], [9.158], [2], [0.0000], [unfavourable],
  [A0 / A11 m=2 p#sub[2]], [8 / 8], [3.595 / 3.995], [1], [0.1905 / 0.1786], [unfavourable],
  [A0 / A11 m=3 p#sub[2]], [3 / 3], [7.468 / 8.742], [2], [0.5000 / 0.0000], [unfavourable],
  [A0 m=4 p#sub[1]], [2], [10.773], [3], [1.0000], [inconclusive (2 cells)],
  [A11 m=4 p#sub[1] (with C2)], [3], [14.886], [3], [1.5000], [unfavourable],
)

Observations:
- Every decided series is cost-unfavourable. That is 23 series with at least 3 completed cells (22 at m ∈ {2, 3}, plus A11 at m = 4), so the prediction held. κ exceeds m − 1 by 2.5–5 bits per bit at m = 2 and by 5.5–7.2 at m = 3.
- σ and κ disagree on which arms look good. A1z, A1g, N11, H2 and A11 at m = 3 have σ = 0 or near it, yet their κ is at or above A0's (7.468).
- A11's κ is above A0's at m = 2 (3.995 vs 3.593) and at m = 3 (8.769 vs 7.468). Its field_ops are lower at small t. Per bit they grow faster.
- The p#sub[2] κ equals the p#sub[1] κ to within 0.03, consistent with the identical traces of Amendment 2.

== C2, the third m = 4 cell (p#sub[1], 3600 s budget)

The lab runner (isolab) has no CLI or credentials in this container, so the cell ran here under the same 13 GB limit. This was registered as a blocker in Amendment 4.
- *A11, t = 4:* both systems completed, each with D = 6 and #raw("max_cols") 88 740, with field_ops 8.630·10#super[11] in 1019.3 s and 1044.1 s. The cell joins F4's series. A11 at m = 4 now has 3 cells (D = 3, 5, 6): σ = 1.5000, which reads *no middle ground*, and κ = 14.886 against 3, *cost-unfavourable*.
- *A0, t = 4:* both systems aborted on memory under #raw("ulimit -v 13000000"), at 536.9 s and 512.8 s, with no D recorded. A0 at m = 4 stays at 2 cells: *inconclusive*.
- The prediction (A11 completes, A0 does not) held.

== C3, p#sub[3] = 40961 as an instrument: not admitted

194 completed p#sub[1] systems with t ≤ 13 (first run, F2, F3, F4, C2, C4) were rebuilt with p#sub[3] seeds. *192 of 194 match* on D, #raw("max_cols"), #raw("max_rows") and #raw("steps"). The 2 mismatches are both A1g, m = 3, t = 3. On p#sub[3] both targets have #raw("max_cols") 5 808 and #raw("max_rows") 5 830. On p#sub[1] the two targets have 5 288 / 5 127 and 5 783 / 5 805. D (8) and #raw("steps") (26) agree.

The rule requires every system to match, so p#sub[3] is *not admitted* and the extension (A0 m = 3, t = 5; A11 m = 3, t = 6; A0 and A11 m = 4, t = 4 on p#sub[3]) was not run. The prediction did not hold. Observed alongside: A1g's set {h·g#super[k]} has no fixed coset shape (ord(g) is drawn), and its two p#sub[1] targets already differ from each other in width. Every arm built on the coset g·μ#sub[2#super[t]] matched. #raw("field_ops") is not comparable between p#sub[3] and p#sub[1]: the narrow path counts dense row lengths, while the wide path counts non-zero updates.

== C4, the A11 ablation

#figure(image("ablation.svg", width: 92%),
  caption: [D, #raw("max_cols") and #raw("field_ops") against b at m = 2 (left) and m = 3 (right) for A0, A11 and the hybrids H1 (A11's coefficient equations with Kummer towers on the roots) and H2 (A0's S#sub[3] with A11's squaring chain). Completed cells, mean of 2 systems. Open markers: H1 and H2 systems stopped at the 300 s budget. Run run-20261016T0616Z, #raw("plot.py").])

At m = 3 (D; log#sub[2] field_ops; max_cols per cell):

#table(columns: 5, align: center, stroke: 0.4pt + gray,
  [*t*], [*A0*], [*A11*], [*H1*], [*H2*],
  [2], [5; 18.95; 368], [4; 12.19; 62], [5; 21.17; 710], [6; 14.05; 160],
  [3], [6; 27.83; 3 285], [4; 22.19; 781], [5; 28.89; 2 939], [6; 26.10; 2 853],
  [4], [6; 33.89; 12 314], [4; 29.67; 3 096], [6; 36.49; 23 053], [6; 32.36; 24 842],
  [5], [timeout], [5; F3, 3600 s budget], [timeout (D ≥ 6)], [timeout (D ≥ 6)],
)

*Rule outcome: mixed, both on D and on log#sub[2] field_ops* (3 shared t). The prediction "the curve relation carries it" did not hold.
- H1 is not cheaper than A0. Its field_ops are higher at every t, and its D is 5, 5, 6.
- H2 has A0's D (6) or more, and field_ops between A11's and A0's at every t. Its width is below A0's at t = 2 and 3 and above it at t = 4 (24 842 vs 12 314).
- Neither hybrid alone reaches A11's D or field_ops. In these rows, A11's lower cost appears only when the function relation and the chain are together. At m = 2 the four arms are within one unit of D of each other and close in field_ops (table in the analysis).

== Predictions scorecard (Amendment 4)

#table(columns: 3, align: (left, left, left), stroke: 0.4pt + gray,
  [*prediction*], [*observed*], [*held?*],
  [C1: every arm cost-unfavourable at m ∈ {2, 3}], [22 of 22 decided series unfavourable], [yes],
  [C2: A11's cell completes and A0's does not], [A11 completed (D = 6), A0 aborted on memory], [yes],
  [C3: p#sub[3] admitted], [192 of 194 match; A1g m = 3, t = 3 differs], [no],
  [C4: the curve relation carries it], [mixed on D and on field_ops], [no],
)

== Deviations and limits

- C2 ran on this container, not on a larger machine (blocker above). A0's t = 4 cell ended on the unchanged memory limit, so A0 at m = 4 still has 2 cells.
- #raw("field_ops") counts row-reduction multiplications only. Symbolic preprocessing, memory traffic and linear algebra over the relations are not counted. κ's break-even (m − 1) assumes relation collection dominates. The verdicts hold at the sizes measured and are not extrapolated.
- Within one cell, the two systems' #raw("field_ops") differ slightly (for example 2 080 268 927 and 2 080 135 234) while D and widths agree. The cell mean is used, as registered.
- The run had no isolation beyond one system at a time. Time-limit stops depend on the host.
