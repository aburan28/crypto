#set page(paper: "a4", margin: 1.7cm)
#set text(size: 9.3pt)
#set table(stroke: 0.4pt)
#show table: set text(size: 8.3pt)
= Prime-degree binary IC, three leads: all negative

2026-10-15. `research/pdp_shape_replay_20261015/`. The protocol was registered (`d0b27a7f4`) before the claim ladders. *Class: engineering, a stage diagnostic.*

F4Remake's idea for the roots of the inherited-F4 decomposition oracle (Koblitz curves over $FF_(2^n)$, $n$ prime). Rows that are dependent on earlier rows are recorded on the first root of a shape (layout, rows, degree, policy) and dropped from later roots of that shape. A replayed root spans a subspace, so the oracle stays exact and may split more. Recording is charged.

#table(columns: 6, align: (left, right, right, right, right, left),
  [ladder], [cells], [default], [replay], [default ÷ replay], [verdicts],
  [frozen (pilot)], [6], [45,517,794], [44,815,209], [1.016], [identical],
  [`chain`], [3], [7,103,844], [7,322,495], [*0.970*], [identical],
  [`chain-holdout`], [4], [43,501,345], [45,579,107], [*0.954*], [identical],
  [`holdout`], [2], [41,414,967], [42,137,329], [*0.983*], [identical],
)

Replay costs more on all 9 registered cells (0.753–0.985×) and splits more on 6 of them. A counting build showed why. The large rank-deficient roots rarely repeat their layout across targets, so they are recorded but not reused: on `K_0/2^15`, $m = 4$, 18 roots were recorded and 6 replayed. Where roots do repeat, only about 4 rows are dropped: on `K_1/2^15`, $m = 4$, 1,832 replays dropped 3.8 rows each.

All predictions held: exact on 9/9 cells, every ratio in $[0.95, 1.10]$, and no gain. By the registered decision this lead is a *recorded null*: `KIC_F4_ROOT_REPLAY` stays off and the oracle is unchanged. The algebraic oracle remains far behind the root-finding table oracle (4.1·10#super[7] word XORs against 592 multiplications per target at (23, 8)), which itself loses to rho, so no scoreboard row moves.

#figure(image("figures/ratio_by_cell.svg", width: 100%), caption: [Default ÷ replay word operations per cell. Values above 1 mean replay is cheaper.])

== Leads 2 and 3

*Lead 2* (a weighted or block grading of the symmetrised model) was already measured on 2026-09-12. The systems are multilinear in their blocks, so per-block bounds reproduce the plain matrix (1.00×) or enlarge it (0.21×). It was not rerun.

*Lead 3* (`decompose_by_e3`, exact $m = 3$ by enumerating $e_3$) now counts its field multiplications. Against the table oracle's $5|F| - 3$ multiplications per target it costs *553–34,901×* on 8 cells $(n, ell)$ with 16 uniform targets each. It is a sealed negative, and no ecbench engine is built (`PROTOCOL-lead3.md`, registered before the run).

#figure(image("figures/e3_vs_table.svg", width: 85%), caption: [$F_(2^n)$ multiplications per target: $e_3$ enumeration (counted) against the pair table.])

Environment: Claude Code cloud container, Xeon \@ 2.10 GHz, 4 vCPU (shared). Counts are deterministic. Conductor task T-3.
