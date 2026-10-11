#set page(paper: "a4", margin: (x: 1.8cm, y: 1.8cm))
#set text(size: 9.5pt)
#set par(justify: true)
#set heading(numbering: "1.")

#align(center)[
  #text(size: 14pt, weight: "bold")[Lead 16.6: the affine saturation law on a genus-2 Jacobian]\
  Round 4 of `RESEARCH_ISOGENY_FIBERS_ECDLP.md`, 2026-10-11. Code: `examples/r4_16_6_genus2_census.rs`.
  Raw records: `research/isogeny_fibers_20261010/round4/16_6_genus2_law/raw/`.
]

= Question and verdict

*Question (§16.6).* Prop. 15.2 prices factor-base symmetry for the mitm
pipeline as cost factor $|M|^(m-2) slash (|A| dot i)$. Does that hold on the
Jacobian of $y^2 = x^5 + b$ over $FF_p$ ($p equiv 1 med (mod 10)$), where the
automorphism group $A = ⟨ zeta_5, iota ⟩$ has order 10 and rational
(2,2)-kernels exist? *Falsifier:* a measured factor outside the law beyond
sampling error.

*Verdict: confirmed, with one structural difference.*
+ All arms follow the law (measured, $p = 4111$; predicted in brackets). The
  m = 3 cost ratios are: A-closed 0.198--0.200 [0.2], μ = 2 1.89--1.98 [2],
  Richelot μ = 4 3.82 [4], non-isotropic μ = 4 3.82--3.92 [4], μ = 5 4.85 [5],
  combined A + M₅ 0.977 [1]. The class factors are 1.00, 0.255, 0.064 and
  0.041 [1, 0.25, 0.0625, 0.04]. The image arms reach $c_2$ ratio 16.0 [16]
  and 5.02 [5].
+ One leading-order deviation, at m = 2 on the A-closed base when A fixes the
  cofactor classes: measured +32%, +21% and +9% at $|F| = 50, 100, 200$. It
  matches the derived structured-pair term $1 + (h-1)(|A|-2) slash (2|F|)$ =
  1.32, 1.16 and 1.08. That term vanishes as $1 slash |F|$ and has the same
  form for $j = 0$ elliptic curves, so it is not a new regime.
+ Genus 2 differs in *realizability*. A base of degree-1 divisors is never
  closed under a nonzero translation, because $|(C - infinity) inter (C - infinity + t)| <= 2$.
  Kernel-saturated bases must therefore be $(mu - 1) slash mu$ degree-2
  divisors. That is fine for the mitm census, but unusable by a
  Gaudry/Nagao oracle. Automorphism folding stays free.

= Setup

Two families, three primes each, all validated: the exhaustive Mumford
count equals $L(1)$, $[\#J]D = 0$, associativity holds, Cantor agrees with
`jv_cover::Hyp`, and $zeta$ acts on $G$ as a 5th root of unity mod $r$.
*C2* has $-b$ a 5th power, so $J[2] subset J(FF_p)$ and $h = 16$: $p = 1051, 2131, 4111$
with $b = 1042, 32, 2172$. *C1* has $b$ square, $T_5 = (0, sqrt(b)) - infinity$
and $h = 5$: $b = 49, 4, 3984$. The census is the one of §15.2. $D_m$ is the
exact number of distinct nonzero classes of $G$ that are $m$-sums of base
elements. Columns are orbits of $[h]F$ under $A$. The cost is
$U dot r slash D_m dot L_m$. Each cell has 40 random bases; the figures show
ratios of means with 95% bootstrap intervals. Costs use the plain base's
nominal $|F| slash 2$ columns, because the ζ-stable x-set folds a few plain
pairs by accident (plain $U$ = 0.83 / 0.92 / 0.95 of nominal).

#figure(
  image("law_measured_vs_predicted.svg", width: 100%),
  caption: [Measured ÷ predicted for every arm and curve (1 = exactly the
  law). Points: ratio of means over 40 random bases; bars: 95% bootstrap CI;
  colour = p. The G arms fall below 1 in cost because the G-filtered base is
  partly A-closed by itself (U ratio 0.32–0.82), and their $c_3$ falls below
  1 at small $r$ (hit rate up to 0.54). Records: `raw/census_*.json`, data
  `plot_data.json`.],
)

#figure(
  image("A_pair_correction.svg", width: 72%),
  caption: [A-closed base, class-factor ratio at m = 2 against $|F|$. On C1
  (A fixes $H = ⟨ T_5 ⟩$) the pairs ${zeta^a P, zeta^b (-P)}$
  always land in $G$, which gives the derived term $1 + 16 slash |F|$ (dashed;
  `finite_size_A.json`). On C2 ζ fixes no nonzero cofactor class and the
  ratio is 1. Measured: `raw/census_c1_p4111_s{50,200}_Aonly.json` and the
  main runs.],
)

= Key measured values (p = 4111, ratio to plain, [prediction])

#table(
  columns: 7,
  stroke: 0.4pt,
  [*curve, |F|*], [*arm*], [*U*], [*$c$, m = 2 strict*], [*$c$, m = 3*], [*cost m = 2*], [*cost m = 3*],
  [C2, 160], [A], [0.217 [0.2]], [0.98 [1]], [0.999 [1]], [0.203 [0.2]], [0.200 [0.2]],
  [C2, 160], [M2], [0.523 [0.5]], [0.503 [0.5]], [0.255 [0.25]], [0.958 [1]], [1.888 [2]],
  [C2, 160], [M4 Richelot], [0.266 [0.25]], [0.232 [0.25]], [0.064 [0.0625]], [1.043 [1]], [3.825 [4]],
  [C2, 160], [M4 non-iso.], [0.266 [0.25]], [0.241 [0.25]], [0.064 [0.0625]], [1.013 [1]], [3.819 [4]],
  [C2, 80], [G, i = 16], [0.611 [1]], [15.98 [16]], [15.27 [16]], [0.050 [0.0625]], [0.038 [0.0625]],
  [C1, 100], [A], [0.211 [0.2]], [1.209 [1.16 corr.]], [1.010 [1]], [0.165 [0.172 corr.]], [0.198 [0.2]],
  [C1, 100], [M5], [0.209 [0.2]], [0.205 [0.2]], [0.041 [0.04]], [0.987 [1]], [4.851 [5]],
  [C1, 100], [A + M5], [0.042 [0.04]], [0.205 [0.2]], [0.041 [0.04]], [0.199 [0.2]], [0.977 [1]],
  [C1, 100], [G, i = 5], [0.824 [1]], [5.022 [5]], [4.906 [5]], [0.156 [0.2]], [0.159 [0.2]],
)
Intervals for every cell and the p = 1051 and 2131 runs are in `RESULTS.md`
and `plot_data.json`. The saturated m = 3 costs climb toward the prediction
as $r$ grows (M5: 4.38 → 4.77 → 4.85), as in §15.2.

= Limits

Toy sizes ($p <= 4111$, $r <= 3.3 dot 10^6$, $|F| <= 200$); the mitm cost
*model*, not solves; no algebraic oracle; no explicit Richelot image arm
(the image arm is $[2]J = G$ on C2). There is no claim about any deployed
curve and no sub-rho claim. Wall time was not used. The commands are in
`run_all.sh`.
