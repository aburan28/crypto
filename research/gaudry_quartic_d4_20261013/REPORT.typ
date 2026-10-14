#set page(paper: "a4", margin: 1.7cm)
#set text(size: 9.3pt)
#set table(stroke: 0.4pt)
#show figure.caption: it => it.body
#show table: set text(size: 8.3pt)

= C8: the $D_4$ symmetry cuts $C_4$ 48×; the $k = 4$ crossover moves from $2^150$ to $2^135$

2026-10-13. `research/gaudry_quartic_d4_20261013/`; code `gaudry_quartic_d4.rs`. Protocol registered before the measured runs. *Class: engineering.*

== Construction

A curve over $FF_(p^4)$ with a rational 2-torsion point $T$ and FGHR's coordinate $Y$, $Y(P + T) = -Y(P)$; base ${P : Y(P) in FF_p}$. $S_5$ in $Y$ is invariant under $D_4 = (ZZ\/2)^3 ⋊ S_4$ (order 192). In $sigma_k = e_k (Y^2)$, $pi = Y_1 Y_2 Y_3 Y_4$ and $sigma_4 = pi^2$ the system is four quartics and $pi^2 - sigma_4$: Bézout 512 against 4,096, Macaulay degree 14 against 29. The invariance is checked by interpolation on fresh points.

== Result, $p = 269$

#table(columns: 4, align: (left, right, right, right),
  [], [$S_4$, $x$-line], [$D_4$, degree 14], [$D_4$, degree 15],
  [quotient dimension], [4,096], [512], [512],
  [matrix], [41,780 × 40,920], [13,412 × 11,628], [18,527 × 15,504],
  [$C_4$ ($FF_p$ mult.)], [$1.213$–$1.227 dot 10^12$], [*$2.535$–$2.552 dot 10^10$*], [$5.13 dot 10^10$],
  [wall a solve], [164–203 s], [4.2–5.6 s], [8.3–9.7 s],
)

$C_4$ falls *47.5–48.4×*. Decompositions equal meet in the middle on all 19 solved targets; one planted target of 20 is degenerate (quotient dimension 517) at degrees 14 and 15. A uniform target decomposes on the $Y$-line base at 0.133 (pooled, two curves, 4,000 targets each) of the $x$-line rate at equal base size: the 8 even $T$-translates of a decomposition hit the same target.

#figure(image("figures/c4_d4.svg", width: 100%), caption: [Left: $C_4$ per decomposition (median). Right: §11.16's model, $S\/"rho" = sqrt(h) (181 (p\/4) C_4 \/ (126 p^2) + r_infinity \/ 4)$, with the measured inputs.])

== Crossover

#table(columns: 4, align: (left, right, right, right),
  [], [$n^*$, $h = 1$], [$n^*$, $h = 2$], [limit past $n^*$],
  [$S_4$, $x$-line], [$2^149.4$ … $2^152.7$], [—], [0.34 … 0.63],
  [$D_4$, $Y$-line], [$2^132.9$ … $2^133.3$], [*$2^135.1$ … $2^135.8$*], [0.12 … 0.22 ($h = 2$)],
)

Rational 2-torsion forces $h >= 2$, so rho runs in a subgroup at least 2× smaller: the crossover moves 14–17 bits earlier. The model's exponents are extrapolated; no end-to-end $S$ is measured. Next: the $k = 4$ ladder in ecbench.

Environment: Claude Code cloud container, Xeon \@ 2.10 GHz, 4 vCPU, rustc 1.97.0. Conductor unreachable: no task card.
