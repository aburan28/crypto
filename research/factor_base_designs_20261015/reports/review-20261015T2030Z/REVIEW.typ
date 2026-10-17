#set document(title: "Factor-base designs: review of run-20261015T1657Z", author: "factor-base designs study")
#set page(paper: "a4", margin: (x: 1.8cm, y: 1.8cm), numbering: "1")
#set text(size: 9.5pt)
#set par(justify: true)
#show raw: set text(size: 8pt)

= Review of run-20261015T1657Z: what the preregistered verdicts do and do not show

*Kind:* solver-axis stage diagnostic (AGENTS.md §2). Interpretation of the
committed rows of #raw("runs/*.jsonl") against
#raw("PREREGISTRATION.md") (§1–§5, Amendment 1). It changes no row, verdict or
threshold. The run report is #raw("reports/run-20261015T1657Z/REPORT.pdf").
Nothing here prices a pipeline, reports S or shows a speedup.

== Verdicts as preregistered

Every series read *candidate* (σ ≤ 0.25 on ≥ 3 completed cells) except A0,
N and A4 at m = 3 (σ = 0.5, *no middle ground*). Those labels stand as the
protocol's output. The four observations below bear on what they mean.

== 1. At m = 2 the threshold does not discriminate

The incumbent A0 reads *candidate* at m = 2 (σ = 0.19, D = 4,4,4,4,5,5,5,5
over t = 2–9), and so does the null N. The σ ≤ 0.25 threshold came from the
m = 4, n = 131 cost model (README §1, #raw("BreakEven.lean")). At m = 2 every
presentation measured, structured or not, has a slope under it. A *candidate*
reading at m = 2 therefore carries no information about an arm.

== 2. The null equals the incumbent: D measures the presentation's shape

N draws a different curve, set and target in every cell (checked: #raw("curve"),
#raw("x_r") and #raw("planted_x") differ on all 26 row pairs) and gives
different F4 traces (#raw("steps") at m = 2, matrix sizes at m = 3). Its D
equals A0's in every cell at both arities. At these sizes D is set by the
system's support and degrees, not by the curve. So differences between arms
are differences between *presentations*, which is the question this study
asks. No arm other than A0 has a null of its own shape, so none of them can
be credited with structure beyond its shape.

== 3. D does not track cost for A1z and A1g

A1z and A1g sit at D = 8 in every completed cell. That is the input degree of
their systems (the substituted x is quadratic, so S#sub[3] has degree 8),
hence σ = 0. Their matrices are the widest of all arms at every t, and their
time grows fastest. At m = 2 and t = 6, the widest matrix has 42 652 columns
against A0's 2 451, and the time is 110 s against 0.3 s (target 0). Both
series hit the memory limit first. For this arm σ on D is a saturated metric.
The *candidate* label is the protocol's literal output, and these rows give
no reason to expect it to replicate as a cheaper design. The preregistered
prediction A1z = A0 ± 1 failed (+2 to +4). A1g = A1z ± 1 held (0 in every
cell): dropping the smoothness condition costs nothing in D here.

== 4. A11 at m = 3: a constant-factor shift, not a flatter slope

At m = 3, A11 is the one arm below the incumbent at equal b on both D and matrix width.

#table(columns: 7, align: center, stroke: 0.4pt + gray,
  [*m, t*], [*A0 D*], [*A11 D*], [*A0 width*], [*A11 width*], [*A0 time*], [*A11 time*],
  [3, 2], [5], [4], [368], [62], [0.0 s], [0.0 s],
  [3, 3], [6], [4], [3 285], [781], [0.6 s], [0.1 s],
  [3, 4], [6], [4], [12 314], [3 096], [16.5 s], [1.8 s],
  [3, 5], [≥ 6 (timeout)], [≥ 5 (timeout)], [≥ 97 456], [≥ 40 486], [300 s], [300 s],
  [2, 9], [5], [4], [14 922], [13 579], [92.0 s], [82.9 s],
)

Its width grows by a factor of ×12.6 and then ×4.0 per bit at m = 3, against
×8.9 and ×3.7 for A0. Its growth rate is no flatter: the curves in the graph
below run parallel, about 4× lower. At m = 2, A11's width meets A0's by
b = 6–9, and both time out at t = 10. At m = 3 the t = 5 cell timed out at a
partial D of 5, so A11's D is not flat beyond the 3 completed cells. Its
σ = 0 rests on the minimum number of cells the protocol accepts, on one prime,
with no null of its own shape.

#figure(image("cols_vs_b.svg", width: 100%),
  caption: [Largest F4 matrix width (#raw("max_cols"), log scale) against b, per arm, target 0 of each cell. Open markers: systems stopped at the 300 s budget, at the width reached (a lower bound). Rows from #raw("runs/*.jsonl"), plotted by #raw("plot_cols.py"). The D graph is in the run report.])

== Predictions scorecard (§3)

#table(columns: 3, align: (left, left, left), stroke: 0.4pt + gray,
  [*prediction*], [*observed*], [*held?*],
  [A1z: D within ±1 of A0, σ no smaller], [+2 to +4; σ = 0 (saturated at the input degree)], [no],
  [A1g = A1z ± 1], [0 in every shared cell], [yes],
  [A4 = A0(⌈t/2⌉) + ≤ 1 (depth)], [+1 to +2], [no],
  [A4 = A0(t) (Bezout)], [+1 in 8 of 9 shared cells, 0 in one (m = 3, t = 3)], [no; closer to Bezout + 1],
  [A11: no prediction], [lowest or tied-lowest D at every shared t (tied with A0 at m = 2, t = 5)], [—],
)

== What is owed before any further statement

+ §5 replication of every *candidate* on p#sub[2] = 2013265921. Only A11 at m = 3 is worth the cost, because points 1 and 3 empty the other labels.
+ A null of A11's shape: the same norm-chain support with random coefficients and a planted zero. Without it, A11's lower D cannot be told apart from the generic behaviour of a quadratic chain.
+ The A11 m = 3, t = 5 cell completed (a longer budget), giving a fourth point.
+ An m = 4 series, where the threshold was derived and A0's slope is about 1 (#raw("RESEARCH_DREG_MEASUREMENT.md")).

Each of these needs its own dated amendment before it runs.
