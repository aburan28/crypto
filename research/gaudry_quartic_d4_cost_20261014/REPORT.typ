#set page(paper: "a4", margin: 1.7cm)
#set text(size: 9.3pt)
#set table(stroke: 0.4pt)
#show figure.caption: it => it.body
#show table: set text(size: 8.3pt)

= A cheaper $k = 4$ decomposition over $FF_(p^4)$: $C_4$ $2.5 dot 10^10 -> 6.9 dot 10^7$, $S$ ÷ 400, crossover (model) $2^135 -> 2^103$

2026-10-14. `research/gaudry_quartic_d4_cost_20261014/`. Protocol registered (commit `2be9ad3d7`) before any claim data; code frozen at `2b49f38d1`. *Class: engineering:* a constant in one decomposition. The exponent $q^(2-2\/n)$ and the asymptote past the crossover are unchanged.

== What changed

The `gaudry-d4` oracle (FGHR's $D_4$ system on a curve with rational 2-torsion) has four options. With all of them off it is the committed oracle, operation for operation (its committed counts were reproduced exactly).

- *weighted*: substitute $sigma_4 = pi^2$: four equations in $(sigma_1, sigma_2, sigma_3, pi)$ with weights $(2,2,2,1)$, weighted degree 8. The Hilbert series is $(1+t^2+t^4+t^6)^3 (1+t+dots+t^7)$, so the weighted Macaulay degree is 26: a $3979 times 4200$ matrix against $13412 times 11628$.
- *chains*: only 44 of the 512 products $pi b$ leave the standard monomials. $M_pi$ is block-companion along $pi$-chains, so the step needs 44 border normal forms, $det A(lambda)$ of a $44 times 44$ polynomial matrix for the characteristic polynomial, and its kernels for the eigenvectors.
- *sparsest*: the sparsest candidate row is each elimination column's pivot.
- *replay* (F4Remake): rows that reduced to zero are skipped. A rank shortfall is topped up from the held rows up to the recorded generic rank, which a specialised system cannot exceed, so the result is exact.

== One decomposition ($C_4$, $FF_p$ multiplications; mean, 95 % bootstrap)

#table(columns: 3, align: (left, right, left),
  [configuration ($p = 269$, 10 planted)], [$C_4$], [interval],
  [committed (plain, degree 14)], [$2.547 dot 10^10$], [2.544–2.550],
  [+ sparsest], [$9.67 dot 10^8$], [9.32–9.96],
  [+ sparsest + chains + replay], [$7.54 dot 10^8$], [7.49–7.59],
  [weighted (degree 26)], [$1.459 dot 10^9$], [1.433–1.483],
  [weighted + chains + sparsest], [$1.400 dot 10^8$], [1.392–1.407],
  [*weighted + chains + sparsest + replay*], [*$6.89 dot 10^7$*], [6.85–6.93],
)

Paired ratios: committed/new 369.8× (367.9–371.6). Sparsest pivots alone give 26.4×: most of the headline is pivot choice, which a competent eliminator on the old system would already have. The weighted form against the plain one at equal linear algebra gives *10.95×* (10.90–11.02). Across $p = 29, 269, 1009, 4049$ the new $C_4$ is $5.5$–$8.5 dot 10^7$, 311–430× below the committed oracle on paired targets. All 648 closed solves equal the meet-in-the-middle enumeration, and paired targets give identical solution sets. 1 % of uniform targets at $p = 29$ do not close (weighted quotient of dimension 509–511). They are charged and dropped: re-solving one would cost about 370 new solves.

#figure(image("figures/c4_ablation.svg", width: 100%), caption: [$C_4$ by option at $p = 269$ (left, log scale) and its split by phase (right). `raw/summary.json`.])

#figure(image("figures/c4_ladder.svg", width: 68%), caption: [$C_4$ against $p$; planted and uniform targets; the committed oracle on the same targets.])

== End to end (ecbench session `ECBS1h4687f1d786c2`, $E(FF_(29^4))$, $\#E = 16 dot 44263$)

#table(columns: 5, align: (left, right, right, left, left),
  [arm], [target 0 $S$], [target 1 $S$], [trials], [unsolved (charged)],
  [rho, negation map], [1.065], [1.412], [], [],
  [IC, committed oracle], [$6.571 dot 10^8$], [$7.228 dot 10^8$], [574 / 635], [9 / 8],
  [*IC, new oracle*], [*$1.653 dot 10^6$*], [*$1.780 dot 10^6$*], [574 / 635], [11 / 13],
)

New/committed: *398× and 406×*, with identical trials and relations. IC still costs $1.39 dot 10^6$× rho (1.26–1.55) at this size. All 8 runs verified. The IC runs graded L1, so wall time (221 s against 7,321 s) is descriptive.

#figure(image("figures/s_by_arm.svg", width: 68%), caption: [$S$ per arm and target, session `k4-p29-d4cost`.])

== Crossover (§11.16's model, extrapolated)

$S\/"rho" = sqrt(h) (181 (p\/4) C_4 \/ (126 p^2) + r_infinity \/ 4)$, with $C_4$ measured at $p = 4049$ and root finding extrapolated in $log p$. At $h = 2$ the crossover moves from $2^135.1$–$2^135.8$ to *$2^103.0$–$2^103.7$* (φ = 0.73 … 1); at $h = 1$, $2^100.7$–$2^101.2$. Past it the ratio still tends to 0.12–0.22× rho. This is a model, not a measured $S$.

#figure(image("figures/crossover.svg", width: 68%), caption: [The model with the committed and the new $C_4$.])

== Predictions and limits

All registered predictions held: A1 exact; A2 $C_4$ in range and ≥ 200× below; A3 ≤ 3 % unsolved; A4 26× and 11×; B2 identical trials; B3 398×, 406×; C1 $2^103$. B1: verified, with the replay audit in `audits/`.

What did not work: weighted degrees below 26; chains on the plain form, where $sigma_1$ has 260 chains and $det A$ costs more than the dense characteristic polynomial; replay without top-up, which fell back on 42 % of systems at $p = 29$; and solving the degenerate systems.

Not tried: F4 with weighted sugar, Wiedemann, Magma on the weighted system.

No `docs/bounds` domain exists for $FF_(p^4)$ curves, so this is a paired session, not a challenge verdict.

Environment: Claude Code cloud container, Xeon \@ 2.10 GHz, 4 vCPU, rustc 1.97.0. Conductor task T-1.
