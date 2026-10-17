#set document(title: "Factor-base designs: review of Amendment 2", author: "factor-base designs study")
#set page(paper: "a4", margin: (x: 1.8cm, y: 1.8cm), numbering: "1")
#set text(size: 9.5pt)
#set par(justify: true)
#show raw: set text(size: 8pt)

= Review of Amendment 2 (F1–F4, run-20261015T2052Z)

*Kind:* solver-axis stage diagnostic (AGENTS.md §2). Interpretation of the
committed rows #raw("runs/F*.jsonl") and #raw("runs/analysis-amendment2.txt")
against #raw("PREREGISTRATION.md"), Amendments 2 and 3. It changes no row,
rule or verdict. Run report: #raw("reports/run-20261015T2052Z/REPORT.pdf").
Nothing here prices a pipeline, reports S or shows a speedup.

== Checks made on the rows, independently of the lane

- All 76 completed F-series systems satisfy their planted solution and are consistent.
- The p#sub[2] rows carry #raw("p = 2013265921") and their own curves (coefficients ~10#super[9]) and targets.
- Rebuilding the committed p#sub[1] cell A11 m = 3, t = 2 with the new harness gives a byte-identical row (#raw("runs/regression-A11-m3-t2.txt")).
- The p#sub[2] instrument check passed (#raw("runs/f1-instrument-check-p2.txt")).
- Amendment 3 corrected a wrong statement in Amendment 2. At m = 4, G is monic. It is not of leading coefficient 1 − p#sub[2]#super[2]. It was registered before any F1–F4 row.

== Decision rules (Amendment 2), as recorded

#table(columns: 3, align: (left, left, left), stroke: 0.4pt + gray,
  [*rule*], [*outcome*], [*prediction*],
  [F1: A11 replicates on p#sub[2], m = 3], [*replicates*. At t = 2, 3, 4: D 4 vs 5, 6, 6. Widths 62 / 781 / 3 096 vs 368 / 3 285 / 12 314], [held],
  [F2: shape or structure], [*shape, not structure*. N11's D equals A11's in all 11 shared cells. Widths are equal at m = 2 and within 7 % at m = 3], [held],
  [F3: m = 3, t = 5 on p#sub[1], 3600 s], [A11 completes at D = 5 (565 s, width 40 486). A0 aborts for memory at ≈ 1 545 s], [held],
  [F4: m = 4, σ], [*inconclusive*: 2 completed cells each (A0 σ = 1.0, A11 σ = 2.0). Both time out at t = 4], [not decided],
)

== What the outcomes mean

+ *The replication is close to automatic.* Every p#sub[2] series has the same matrix widths as its p#sub[1] counterpart, cell for cell, and the same D. At these sizes F4's work is fixed by the system's shape (support and degrees), not by the prime or the curve. That is also why N equalled A0 and N11 equals A11. F1 shows that A11's advantage does not depend on the prime. Given that cost is set by shape, it could hardly have come out otherwise.
+ *A11 is a different shape for the same factor base, not curve structure.* The norm-chain system is quadratic, with 3 + 3t unknowns at m = 2 and 5 + 5t at m = 3, instead of a degree-4 or -6 summation polynomial on 2t or 3t + 1 unknowns. Its lower D and narrower matrices come from that shape. F2 gives no sign of anything beyond it.
+ *The D slope of 0 at m = 3 does not hold.* With the F3 cell as a fourth point, A11's m = 3 series on p#sub[1] reads D = 4, 4, 4, 5: σ = 0.30, *weak* by §1. At t = 5, A11 finishes in 565 s where A0 runs out of memory. That is the largest per-cell gap measured.
+ *The growth per bit is no flatter, and at m = 4 it is steeper.* Δlog#sub[2] width per bit is 2.82 for A11 against 2.53 for A0 at m = 3, and 5.99 against 3.50 at m = 4. At m = 4, A11 starts far below A0 (34 vs 1 911 columns at t = 2) but both time out at t = 4, at a recorded width of 88 740 for A11 and 102 516 for A0. If A11 keeps that growth rate, its matrices are wider than A0's by t = 5. That is an extrapolation from 2 cells, not a measurement.

#figure(image("cols_vs_b_amendment2.svg", width: 100%),
  caption: [Largest F4 matrix width (log scale) against b for the F1 (p#sub[2]), F4 (m = 4) and F2 (N11 beside A11) series. Copied unchanged from the run report (#raw("reports/run-20261015T2052Z/cols_vs_b.svg"), plotted from the rows by its #raw("plot.py")). Open markers: stopped at the 300 s budget.])

== Reading

A11, the norm-divisibility chain, is a better *presentation* of the Kummer
fibre at small depth. It has lower D and matrices several times narrower,
on both primes. One cell finishes where the incumbent exhausts memory. No
measurement shows it changing the slope that the break-even argument
(README §1) needs, and m = 4 points the other way. Under the register's
rule it stays *open*, as an engineering presentation rather than a factor-base
advance. The question it leaves is whether the constant factor survives at
m = 4 for t ≥ 4. That needs a longer budget than the 300 s used here, and
its own amendment.
