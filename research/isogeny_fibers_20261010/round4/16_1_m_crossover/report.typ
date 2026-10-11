#set page(paper: "a4", margin: 1.8cm)
#set text(size: 9.5pt)
#set heading(numbering: none)

= Lead 16.1: where kernel saturation pays, the feasibility crossover in $m$

Round 4 of `RESEARCH_ISOGENY_FIBERS_ECDLP.md` §16.1. Toy Koblitz curves only:
`icv1-f2m17-tm101-00378d4e`, `icv1-f2m19-t797-b6cf2467` and
`icv1-f2m23-t5197-69e76b73`. No sub-rho claim and no claim about any deployed
curve. Full tables, exact commands and limits are in `RESULTS.md` in this
directory.

== Question and falsifier
The affine law (Prop. 15.2) charges a base closed under the rational
2-torsion kernel $M = ⟨T⟩$ a mitm cost factor $|M|^(m-2)$.
The $u$-frame symmetrised systems have lower degree. Is there a smallest
$(n, m)$ where the plain $x$-frame arm is infeasible under a shared counted
cap but the symmetrised arm is not, and does that arm beat the best feasible
arm at the same $n$? *Falsifier:* a crossover arm that beats every feasible
arm in total cost.

== What ran
- *Panel* (`examples/r4_16_1_m_crossover.rs`): the library oracles with
  `InheritedF4` (max degree 3) and the same `node_budget = 4096` for every
  solve. Each $(n, m)$ has one shared $V ∋ 1$ of dimension
  $ell = ceil((n + log_2 m!)/m)$, plus invariant-$V$ cells
  ($ell = 9$ at $n = 17$, $ell = 12$ at $n = 23$). There are 24 subgroup
  targets per curve, and every verdict is gated against exact enumeration
  (0 gate failures in all 25 cells). The unit is engine-counted word XORs, a
  lower bound.
- *Anchors* (`ic bench`): end-to-end with a planted logarithm, verified, and
  8 matched counted rho runs, on the invariant bases, under the same node
  budget. The $m = 3$ arms are capped at 300 trials and 10 min of contended
  wall time.

#figure(
  image("crossover.svg", width: 100%),
  caption: [Measured panel, 24 targets per $n$, node budget 4096. Top:
  $log_(10)$ of word XORs per relation found (hatched: not feasible, i.e.
  fewer than half the attempts solved, or none found). Bottom: solved
  fraction. The $n = 23$, $m = 3$ `x-direct` cells were not run (see
  Limits). At $n = 23$, $m = 4$ both chained layouts exceed the engine's 64
  unknowns. Data: `crossover_data.csv`, from `results/panel_*.json`.],
)

== Key measured numbers
#table(
  columns: 5,
  align: (left, left, right, right, left),
  [*cell*], [*arm*], [*XOR / relation*], [*solved*], [*feasible*],
  [$n=17$, $m=3$, $ell=9$ inv], [x-direct (S#sub[4], deg 6)], [— (1 target, 905 s)], [0/1], [no],
  [], [x-chained], [1.33e6], [24/24], [yes],
  [], [sym (S#sub[4], deg 4)], [2.26e7], [14/24], [yes],
  [$n=19$, $m=3$, $ell=8$], [x-direct], [7.12e9], [9/24], [no],
  [], [x-chained], [4.65e6], [24/24], [yes],
  [], [sym], [3.45e7], [6/24], [no],
  [$n=17$, $m=4$, $ell=6$], [x-chained], [2.36e7], [22/24], [yes],
  [], [sym-chained], [2.75e8], [1/24], [no],
  [$n=23$, $m=3$, $ell=12$ inv], [x-chained], [1.03e7], [23/24], [yes],
  [], [sym], [4.57e9], [1/24], [no],
)

End-to-end (`ic bench`, $S = "GAE" slash sqrt(r)$): at $n = 17$, rho 1.76,
mitm $m=2$ 17.1, x $m=2$ 757, sym $m=2$ 3 112, *sym $m=3$ 74 460*, and the
direct x $m=3$ arm gave no output in 10 min. At $n = 23$, rho 1.09, mitm
$m=2$ 143, *sym $m=3$ 938 966*, and the direct x $m=3$ arm gave no output in
10 min.

== Verdict: confirmed, the falsifier is not met
+ A crossover exists only against the direct S#sub[4] descent. The smallest
  cell is $(n = 17, m = 3)$ on the invariant $ell = 9$ base, where the
  symmetrised arm solves and recovers the logarithm and the direct arm does
  not.
+ Against the chained $x$-frame system, which is already in the library,
  there is no crossover anywhere in $n in {17, 19, 23}$, $m in {3, 4}$. The
  chain is feasible wherever the symmetrised arm is, and 2.4× to 440×
  cheaper per relation. At $m = 4$ the kernel-saturated chain is the
  infeasible one.
+ The crossover arm is 98× (x $m = 2$) and 4 360× (mitm $m = 2$) above the
  best feasible arms at $n = 17$, and far above rho. Kernel saturation buys
  feasibility only against one particular plain system, never total cost.

== Limits
The $x$-direct $m=3$ arm at $n = 23$ was not run in the panel: its first
attempt was killed after 30 min inside one solve. Its infeasibility rests on
the `ic` anchor (no output in 10 min) and on the 905 s solve at $n = 17$, so
it is a wall observation. No direct S#sub[5] exists in either frame, so
$m = 4$ compares chains only. With 24 targets per cell, hit rates carry noise
of about ±0.1. The panel and `ic` disagree on the sym $m=3$, $ell=12$ hit rate
(0.04 against 0.22), which I did not resolve. Word XORs exclude the Macaulay
build, and wall times are contended.
