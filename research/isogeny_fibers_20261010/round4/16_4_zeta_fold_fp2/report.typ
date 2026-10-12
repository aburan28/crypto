#set page(paper: "a4", margin: 1.8cm)
#set text(size: 9.5pt)
#set heading(numbering: "1.")

= Lead 16.4: the free automorphism fold of the Gaudry base over $FF_(p^2)$

Round 4, 2026-10-11. Example `examples/r4_16_4_zeta_fold_fp2.rs`. Raw records in this directory: `census_torsion.json`, `pipeline_{prime,torsion}.json`, `stream_{prime,torsion}.json`, `sweep_torsion_{10,12}.json` (plus `.txt` logs). Full tables, exact commands and limits are in `RESULTS.md`. Toy synthetic curves only. No claim about rho, sub-rho behaviour or deployed curves. Wall times are contended.

== Question and arms

On $E: y^2 = x^3 + b$ over $FF_(p^2)$, $p equiv 1 med 3$, the Weil-descent base $F = {P : x(P) in FF_p^*}$ is closed under $zeta(x,y) = (omega x, y)$. Prop. 15.2 predicts: (b) folding by $⟨-1, zeta⟩$ instead of (a) $⟨-1⟩$ costs $times 1\/3$, and (c) a base saturated by a rational 3-torsion point $T$ is neutral at $m = 2$ and loses $times 3$ at $m = 3$. Arm (c) is $F_c = phi^(-1){Q in E' : x(Q) in FF_p^*}$ for the 3-isogeny $phi: E -> E': y^2 = x^3 - 27 b$ with kernel $⟨T⟩$. It folds by $⟨-1, zeta, tau_T⟩$ (18 points per column). Through $phi$ it is the base $F'' = {Q in phi(E) : x(Q) in FF_p^*}$ on $E'$.

== Structural finding (derived; checked on every instance)

Rational 3-torsion on $y^2 = x^3 + b$ requires $b$ to be a square. Then Frobenius $pi = (zeta_6^j pi_1)^2$ is a square, so $E$ is a base change from $FF_p$ and $\#E = N_1 N_2$ with $N_i approx p$. *The census shape $\#E = 3r$ of §15.2 cannot occur over $FF_(p^2)$.* Arm (c) exists only on curves with $r <= (sqrt(p)+1)^2$. So two families were used. The prime family ($b$ a non-square, $h <= 16$, $r approx p^2$) runs arms (a) and (b). The torsion family ($b$ a square, $9 ∤ \#E$, $h approx p$) runs all three arms on one curve.

== Results

#figure(
  image("ratios.svg", width: 100%),
  caption: [Recorded ratios, one point per curve, against the Prop. 15.2 prediction (red). Left: exact census yields (16 torsion curves, $p approx 2^9$ to $2^(12)$, 1000 targets each), size-matched to $|F_c| = |F|$. Right: cost ratios. $U_b\/U_a$ and total (b)/(a) come from the repository's `LineOracle` (24 curves, $p$ up to $2^(16)$). Trials to full rank come from `full_rank_stream` (6 prime-family curves). Total (c)/(b) at $m=2$ uses the line oracle on 12 torsion curves. Total (c)/(b) at $m=3$ uses the example's sweep oracle (8 curves); the "on $E$" column is derived as $3 times$ the quotient sweep. Sources: the JSON files listed above.],
)

#table(
  columns: (auto, auto, auto, auto),
  align: (left, center, left, left),
  table.header([*comparison*], [*law*], [*measured (range over curves)*], [*status*]),
  [(b)/(a) total, $m=2$, line oracle], [1/3], [0.3333 on 24 of 24 curves (6 points/column, strict ζ-closure)], [confirmed, exact],
  [(a)/(b) relations to full rank], [3 × coupon factor], [2.79–4.19 (6 streams, log recovered on all)], [confirmed],
  [(c)/(b) total, $m=2$, line oracle, size-matched], [1], [0.54–1.11; 8 of 12 per-curve 95% intervals contain 1], [confirmed (neutral; Poisson saturation predicts about 0.85)],
  [(c)/(b) distinct rows, $m=3$], [1/9], [0.100–0.112 (16 curves)], [confirmed],
  [(c)/(b) rows per unknown, $m=3$], [1/3], [0.301–0.336 (16 curves)], [confirmed],
  [(c)/(b) total, $m=3$, sweep over $F_c$ on $E$], [3], [2.67–3.62 (derived: 3 × quotient sweep)], [reproduces the law],
  [(c)/(b) total, $m=3$, quotient sweep on $E'$], [3], [0.93–1.07 size-matched; 95% intervals exclude 3 on 8 of 8], [*law's ×3 refuted for this oracle*],
)

== Verdict

The ζ-fold is free and exact. The base is automatically ζ-closed, the columns drop by exactly 3, and the descent cost per call is unchanged. The repository supports it but had no driver for plain $j = 0$ curves over $FF_(p^2)$. The kernel-saturated arm follows the law's yield predictions exactly: the raw/distinct multiplicities are 3.00 and 9.0, and the size-matched yields are $1\/3$ and $1\/9$. Its $m = 3$ cost penalty of ×3 appears only for an oracle that enumerates over $F_c$. An algebraic oracle in the quotient coordinate $x' = x compose phi$ enumerates $F_c\/⟨T⟩$, and there the arm is measured neutral ($approx 1$). The $|M|^(m-2)$ term of Prop. 15.2 should therefore be stated per oracle.

Limits: toy sizes ($p <= 2^(16)$; census and sweep $p <= 2^(12)$). The $m=3$ oracle is the example's sweep and returns one pair per call. Arm (c) runs only on base-change curves. Torsion-family streams did not reach full rank, though all recovered the logarithm. Inversions are not priced.
