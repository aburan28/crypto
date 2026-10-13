#set page(paper: "a4", margin: 1.7cm)
#set text(size: 9.3pt)
#set table(stroke: 0.4pt)
#show figure.caption: it => it.body
#show table: set text(size: 8.3pt)

= Gaudry's $k = 4$ index calculus with the $D_4$ oracle, end to end in ecbench

2026-10-13. `research/gaudry_quartic_d4_ecbench_20261013/`. Protocol registered before the session; audit OK (3 of 3 runs replayed identically). *Class: engineering.*

== Setup

$E(FF_(29^4))$ with rational 2-torsion, $\#E = 16 dot 44263$ (registered `icv1-fp5k4-tm926-7bc5b359`), one target. `ic.pipeline` with FGHR's $Y$-line base (columns ${plus.minus P, plus.minus (P + T)}$) and the `gaudry-d4` oracle ($S_5$ in the $D_4$ invariants, Macaulay degree 14); every phase charged, the oracle's $FF_p$ multiplications at the curve's 97 per addition. Reference: rho with the negation map.

#figure(image("figures/s_by_arm.svg", width: 100%), caption: [Left: $S$ per arm, with §11.16's model value for the IC arm. Right: the IC arm by phase.])

#table(columns: 4, align: (left, left, right, left),
  [arm], [method], [$S$], [isolation],
  [rho-neg], [`rho.negation`], [1.065], [L2],
  [ic-k4-d4], [`ic.pipeline`, `fghr-y-line`, `gaudry-d4`], [*$7.17 dot 10^8$*], [L1],
  [rho-neg-aa], [`rho.negation`], [1.065], [L2],
)

633 relation trials, 631 solves at $2.32 dot 10^10$ $FF_p$ multiplications each, 8 relations, rank 8 on 8 unknowns; relations are $approx 100 %$ of $S$. *Index calculus costs $6.7 dot 10^8 times$ rho*, within 2× of the model's $1.3 dot 10^9$; every registered prediction held. Three decompositions were lost to a summand outside the base (probably $Y = 0$); 9 of 631 systems were degenerate.

Scope: one small curve and one target; the method working and priced end to end where it was expected to lose by nine orders of magnitude. It confirms the model's relation-phase constant, not the projected crossover near $2^135$.

Environment: Claude Code cloud container, Xeon \@ 2.10 GHz, one reserved core, rustc 1.97.0. Conductor unreachable: no task card.
