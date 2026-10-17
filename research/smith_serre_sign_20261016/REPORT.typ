// Results report. Every number is read from data/*.json; nothing here
// computes a result (rounding and the share of the gap are presentation).
#set page(paper: "a4", margin: (x: 1.6cm, y: 1.5cm), numbering: "1")
#set text(size: 9pt)
#set par(justify: true, leading: 0.55em)
#set heading(numbering: "1.")
#show heading: set text(size: 10.5pt)
#let s7 = json("/data/serre_p7.json")
#let s11 = json("/data/serre_p11.json")
#let f(x, n: 3) = str(calc.round(x, digits: n))
#let pc(x) = str(calc.round(x * 100, digits: 0)) + "%"
#let tbl(..args) = text(size: 7.8pt, table(inset: 3pt, stroke: 0.4pt + luma(170), ..args))
#align(center)[#text(size: 13pt, weight: "bold")[Smith's Prop. 4 sign follows $\#H(FF_p)$] \
`research/smith_serre_sign_20261016` · 2026-10-16 · protocol `PROTOCOL.md` frozen at `3be242141` before any $p = 11$ data]

= Question
Smith's construction fully succeeds on a genus-3 hyperelliptic curve $H$ when a rational tractable subgroup $S$
passes Prop. 2 (the trigonal map is rational) and Prop. 4 (the quotient $J_H\/S$ is the Jacobian of the plane
quartic $X$ rather than its quadratic twist: Serre's obstruction). Smith models both as fair coins. Two earlier
studies found that `full_success` clusters within $L$-classes beyond isomorphism-class granularity
(`smith_iso_granularity_20261015`) and that no 2-adic invariant of $J_H$ explains it (`smith_class_invariant_20261015`).
This study asks whether the Prop. 4 sign is biased by an $L$-class invariant.

= Verdicts (frozen wording; $p = 11$ is confirmatory, $p = 7$ found the law)
#tbl(columns: (10%, 45%, 45%),
  [*prime*], [*H1: sign law (monotone, max $|r - N_1\/(2p+2)| <= 0.10$)*], [*H2: clustering under $("stratum", N_1, "Prop. 2")$*],
  ..(s7, s11).map(s => ([$p = #s.p$#if s.p == 7 [ (expl.)]], [*#s.h1.verdict*: monotone #s.h1.monotone, max deviation #f(s.h1.max_abs_dev)],
    [*#s.h2.verdict*: closes #pc(s.h2.gap_closed_n1) of $U - C$ (stratum alone #pc(s.h2.gap_closed_stratum), random-label control
     #pc(s.h2.gap_closed_random)); null 95% [#f(s.h2.null_stratum_n1_trig.lo, n: 4), #f(s.h2.null_stratum_n1_trig.hi, n: 4)] vs $C = #f(s.h2.C, n: 4)$])).flatten(),
)

#figure(image("/figures/serre_sign.svg", width: 100%),
  caption: [(a) Prop. 4 pass rate $r(N_1)$ among non-degenerate Prop. 2 passes (subgroups where exactly one of $H$ and its
  twist passes Prop. 4), weighted by models, against $N_1\/(2p+2)$; hollow markers are bins under 1% of the Prop. 2 passes,
  outside the test. (b) Share of the gap $U - C$ closed by the median of each 999-permutation null. Source
  `figures/plot.py`, data `data/serre_p*.json`.])

= Findings
- *Measured (both primes).* The quadratic twist preserves Prop. 2 and flips Prop. 4, except on degenerate subgroups
  (Smith's construction fails at small $p$) where Prop. 4 fails on both sides (control C5′ passes; at most
  #pc(calc.max(..s11.h1.bins.map(b => b.degenerate_share))) of Prop. 2 passes per bin at $p = 11$). So for each Prop. 2 pass the
  quotient curve $X$ has either $\#X(FF_p) = N_1 = \#H(FF_p)$ (pass) or $2p + 2 - N_1$ (fail, $X$ in the twist class).
- *Measured, frozen rule (H1).* The pass rate rises with $N_1$ at both primes and, rescaled by $2p + 2$, the two primes
  lie on one curve close to $N_1\/(2p+2)$ (Figure 1a). At $p = 11$: #s11.h1.bins.filter(b => b.in_test).map(b =>
  [$N_1 = #b.N1$: #f(b.r) (model #f(b.model_N1_over_2p2))]).join("; "). The curve is steeper than linear in the tails.
- *Consequence.* $N_1 = p + 1 + a_1$ is an $L$-class invariant, so Prop. 4 outcomes inside a class share a bias. A class
  with few rational points almost always fails Prop. 4, and its twist class almost always passes. Combined with whether any
  subgroup passes Prop. 2, this absorbs #pc(s7.h2.gap_closed_n1) ($p = 7$) and #pc(s11.h2.gap_closed_n1) ($p = 11$) of
  the gap $U - C$, against #pc(s7.h2.gap_closed_random) / #pc(s11.h2.gap_closed_random) for a random class label of the same
  coarseness. The remaining #pc(1 - s11.h2.gap_closed_n1) at $p = 11$ is still unexplained (H2: partial).
- *Interpretation (labelled, not tested).* Between the two candidates for $X$, the one with more rational points is chosen
  roughly in proportion to its point count, as if Smith's rational data on $X$ (the trigonal structure) were attached to a
  rational point. If $r approx N_1\/(2p+2)$ persists, the bias is $|a_1|\/(2p+2) = O(1\/sqrt(p))$ and Smith's fair coin
  (Hypothesis 3) holds asymptotically. At small $p$ it does not: Prop. 4 is a biased coin whose bias is set by the class.
- *For the walk lead.* A walk inside an $L$-class cannot change $N_1$, so it cannot escape this bias. Classes with small
  $\#H(FF_p)$ are poor walk targets, and their twist classes are good ones.

= Validation
#tbl(columns: (auto, 1fr),
  [*control*], [*result*],
  [C1–C3 (2026-10-15 controls: census counts, invariance, partition)], [all pass at both primes],
  [C4: Prop. 2 / Prop. 4 subgroup counts constant on every isomorphism class], [pass ($p = 7$: #s7.controls.c4_subgroup_counts_inconsistent, $p = 11$: #s11.controls.c4_subgroup_counts_inconsistent inconsistent)],
  [C5 (exact twist complementarity; expected to fail on degenerate subgroups)], [#s7.controls.c5_twist_complementarity_violations / #s11.controls.c5_twist_complementarity_violations isomorphism classes with a degenerate subgroup],
  [C5′: curve-side + twist-side Prop. 4 passes ≤ Prop. 2 passes, per $N_1$ pair], [#s7.c5prime_pass / #s11.c5prime_pass],
  [random-label control (protocol)], [#pc(s11.h2.gap_closed_random - s11.h2.gap_closed_stratum) extra at $p = 11$ vs #pc(s11.h2.gap_closed_n1 - s11.h2.gap_closed_stratum) for $N_1$: under the 75% bar, so H2 is not attributed to coarse conditioning],
)

= Scope
Exact enumeration at $p = 7, 11$, degree-8 models. $p = 7$ was used to find the law; $p = 11$ was held out. The
$N_1\/(2p+2)$ form is suggested by the data, not derived. The residual quarter of the clustering is open. Nothing here
builds or prices an isogeny walk, or says anything at cryptographic sizes or in characteristic 2.

= Reproduce
```
cargo build --release --example smith_iso_granularity
B=target/release/examples/smith_iso_granularity; D=research/smith_serre_sign_20261016/data
$B --p 7 --out $D --subgroups; $B --p 11 --out $D --subgroups   # 17 s and 16 min on 2 vCPU
$B summarize-serre --dir $D 7 11
python research/smith_serre_sign_20261016/figures/plot.py         # layout only
```
Source: B. Smith, J. Cryptology 22 (2009), arXiv:0806.2995, Props. 2 and 4, Hypothesis 3.
