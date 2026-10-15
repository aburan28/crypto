#set page(paper: "a4", margin: 1.7cm)
#set text(size: 9.3pt)
#set table(stroke: 0.4pt)
#show figure.caption: it => it.body
#show table: set text(size: 8.3pt)

= $k = 4$ index calculus in characteristic 2: Frobenius-reduced D#sub[4] system, $C_4$ 481× lower, IC ÷ rho falling as $n^(-1\/4)$ (crossover extrapolated near $2^68$)

2026-10-15. `research/gaudry_quartic_d4_char2_20261015/`. Protocol registered (commit `666438375`) before any claim data; ecbench binary `92610bdf…` throughout. *Class: engineering:* a constant in one decomposition, plus a measured trend; the exponent is unchanged and no run reaches the crossover. Every run is listed in `EXPERIMENTS.md`.

== The system

With $T = (0, sqrt(b))$ and $t = beta\/(x + beta)$, $beta = b^(1\/4)$, adding $T$ is $t |-> t + 1$. The invariants are $w = t^2 + t$, $s_k = e_k (w)$ and $pi = sum t$ with $pi^2 + pi = s_1$. Only 21 of 105 monomials carry the interpolated relation, all even in $s_2, s_3, pi$, so the system is rewritten in $(s_2^2, s_3^2, s_4, pi^2)$: weights $(2,2,1,1)$, 64 solutions, Macaulay degree 11, a $384 times 448$ matrix (plain: $11116 times 11628$). Chains, sparsest pivots and replay with top-up come from the odd study (#1863).

== One decomposition ($C_4$, $FF_q$ multiplications; run report `run-20261015T0530Z-ladder`)

#table(columns: 3, align: (left, left, right),
  [configuration ($m = 8$, 10 planted)], [matrix], [$C_4$],
  [plain, degree 14, + sparsest + replay], [$11116 times 11628$], [$9.54 dot 10^7$],
  [weighted (2,2,2,1), degree 26, + s + r], [$3688 times 4200$], [$3.60 dot 10^6$],
  [*Frobenius (2,2,1,1), degree 11, + s + r*], [$384 times 448$], [*$1.98 dot 10^5$*],
)

Paired ratios (95 % bootstrap): plain/weighted 26.5× (26.4–26.6), weighted/Frobenius 18.1× (18.08–18.18), *plain/Frobenius 481×* (479–483). Across $m = 5..15$ the Frobenius $C_4$ is $1.13$–$2.50 dot 10^5$; only root finding grows. 512 closed solves, none disagreeing with MITM; the 8 planted failures ($m = 5, 6$) are positive-dimensional families (every decomposition contains $\{t, t + 1\}$), and no uniform target fails.

#figure(image("figures/c4_ablation.svg", width: 92%), caption: [$C_4$ by option at $m = 8$ and its split by phase.])

== End to end (ecbench; $S$ = group-addition equivalents ÷ $sqrt(r)$)

At $m = 6$ (session `ECBS1h03a1d62aba88`, 10/10 verified) the three forms use identical trials and relations; $S$ is 713,704 / 722,233 (plain), 38,542 / 31,128 (weighted), *4,501 / 3,654 (Frobenius)*, rho 1.017 / 1.314. `ecbench compare`: Frobenius/plain 0.0057 (0.0051–0.0063).

#table(columns: 5, align: (right, right, right, left, left),
  [$m$], [$log_2 n$], [$h$], [IC (Frobenius) $S$], [IC ÷ rho theory],
  [6], [24], [14], [4,501; 3,654], [5,080; 4,124],
  [8], [32], [6], [573.9; 459.6], [648; 519],
  [10], [40], [2], [71.4; 73.4], [80.6; 82.9],
  [12], [48], [10], [44.0; 42.6], [49.7; 48.1],
  [14], [56], [2], [5.37 (one target)], [*6.06*],
)

Sessions `ECBS1h4d63e1551dc1` (24/24) and `ECBS1hb9229e3a0f30` (3/3). IC is still 6.1× rho's theory (4.3× rho's measured run) at $n approx 2^56$. The fit of $log_2 (S_"IC" \/ 0.886)$ against $log_2 n$ for $m >= 8$ has slope *−0.248* ($n^(-1\/4)$ expected) and reaches rho at *$n^* approx 2^68.1$* (leave-one-out $2^66.7$–$2^71.6$): an extrapolation. The model with measured inputs gives $2^63.1$–$2^69.0$ for $h = 1, 2, 4$ ($r_infinity$ assumed). The uneven steps follow the cofactor $h$ (the ratio carries $sqrt(h)$); this reading is post hoc.

#figure(image("figures/crossover.svg", width: 80%), caption: [IC $S$ ÷ rho's theory $S$ per run, the fit ($m >= 8$) and the model. `raw/summary.json`.])

== Predictions

A1 (exact), A2 ($C_4 in [8 dot 10^4, 4 dot 10^5]$), A3 (closing), A4 (ablation), B1 (verified, replays identical, equal trials, order), B2 (ratio ranges at $m = 8, 10, 12, 14$), C1 (slope in $[-0.35, -0.15]$) and C2 ($n^*$ in $[2^60, 2^85]$): *all held*. The replay audit (`verify --replay-all`) reproduced all 37 runs of the three sessions.

*Scope.* Curves over $FF_(2^(4m))$ are not standard curves; the crossover is an extrapolation plus a model, and $m = 14$ is one target. Prime-degree binary fields are settled separately (#1885).

Environment: Claude Code cloud container, Xeon \@ 2.10 GHz, 4 vCPU, one core reserved per session. Conductor task T-2.
