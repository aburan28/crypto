#set page(paper: "a4", margin: 1.8cm)
#set text(size: 9.5pt)
#set table(stroke: 0.4pt)
#show figure.caption: it => it.body

= C4 falsifier: exponent-space structure of decompositions on PKM's Kummer tower

2026-10-13. `research/ic_exponent_falsifiers_20261013/`, program `examples/kummer_exponent_bias.rs`. *Class: falsification, negative.*

== Question

On PKM's tower $V = mu_(2^t)$ every base abscissa is $zeta^e$. A representation-technique search for $R = P_a + P_b + P_c$ needs a linear form of the exponents that the target predicts. Simplest form: is $e_a + e_b + e_c mod 2^k$ concentrated among the decompositions of one target?

== Method

Curve `icv1-fp20-tm728-d215a2a1` ($p = 946177$), $t = 8$, 246 signed base points. 2,000 uniform targets; every decomposition enumerated exactly and counted once (9,837 in all). Statistic: mean over targets of the largest residue-class share of the exponent sum mod $2^k$. Null: random exponent triples, same count per target. Positive control: the same detector on $F_p^*$ multiplication, where the sum is fixed by the target.

#figure(image("figures/c4_class_share.svg", width: 90%), caption: [Largest class share per $k$: the curve lies on the null; the positive control reads 1.])

#table(columns: 5, align: right,
  [$k$], [$1\/2^k$], [curve], [null], [positive control],
  [1], [0.500], [0.680], [0.680], [1.000],
  [2], [0.250], [0.508], [0.511], [1.000],
  [3], [0.125], [0.393], [0.396], [1.000],
  [4], [0.0625], [0.321], [0.321], [1.000],
  [5], [0.0313], [0.273], [0.275], [1.000],
  [6], [0.0156], [0.244], [0.248], [1.000],
)

== Decision

The exponent sum of a decomposition carries no information about the target (curve and null differ by under 1 % at every $k$). C4's simplest form is falsified on this base, as predicted: the elliptic group law does not respect the multiplicative structure of $mu_N$. Scope: one curve, one tower height, the plain sum.

Environment: Claude Code cloud container, Intel Xeon \@ 2.10 GHz, 4 vCPU, rustc 1.97.0. Conductor control plane unreachable: no task card.
