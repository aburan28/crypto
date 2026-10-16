// Results report. Every number is read from data/*.json; nothing here
// computes a result (rounding is presentation only).
#set page(paper: "a4", margin: (x: 1.6cm, y: 1.5cm), numbering: "1")
#set text(size: 9pt)
#set par(justify: true, leading: 0.55em)
#set heading(numbering: "1.")
#show heading: set text(size: 10.5pt)
#let v = json("/data/verdict.json")
#let R = (101, 1009, 10007).map(p => json("/data/scaling_p" + str(p) + ".json"))
#let f(x, n: 3) = str(calc.round(x, digits: n))
#let tbl(..args) = text(size: 7.8pt, table(inset: 3pt, stroke: 0.4pt + luma(170), ..args))
#align(center)[#text(size: 13pt, weight: "bold")[The Prop. 4 sign bias shrinks like $1\/sqrt(p)$] \
`research/smith_sign_scaling_20261016` · 2026-10-16 · protocol `PROTOCOL.md` frozen at `5fd6133d0` before any data]

= Question
At $p = 7, 11$ the Prop. 4 pass rate of Smith's construction follows $N_1\/(2p+2)$ with $N_1 = \#H(FF_p)$
(`smith_serre_sign_20261016`), i.e. $r = 1\/2 + beta a_1\/(2p+2)$ with $beta approx 1$. If $beta$ stays bounded, the
bias at a typical curve ($|a_1| tilde sqrt(p)$) is $O(1\/sqrt(p))$ and Smith's fair coin (Hypothesis 3) holds
asymptotically; if $beta$ grows like $sqrt(p)$, the bias persists at every size. This study measures $beta$ on
$10^6$ random curves at each of $p = 101, 1009, 10007$.

= Verdict (frozen wording)
#tbl(columns: (auto, 1fr),
  [*H3: scaling exponent $gamma$ of $beta_p$ in $p$ (WLS on $ln beta$ vs $ln p$)*],
  [$gamma = #f(v.gamma)$, 95% CI [#f(v.gamma_ci.at(0)), #f(v.gamma_ci.at(1))] → *#v.h3* (rule: upper end < 0.25; persistence predicts 0.5)],
)
#tbl(columns: 7,
  [*p*], [*squarefree curves*], [*non-degenerate Prop. 2 passes*], [*degenerate*], [*$alpha$ (control: $0.5 plus.minus 0.02$)*],
  [*$beta$ [95% CI, CR1]*], [*H1′: CI contains 1*],
  ..R.map(r => ([#r.p], [#r.squarefree], [#r.nondegenerate], [#r.degenerate], [#f(r.alpha, n: 4)],
    [#f(r.beta) [#f(r.beta_ci.at(0)), #f(r.beta_ci.at(1))]], [#r.beta_ci_contains_1])).flatten(),
)

#figure(image("/figures/sign_scaling.svg", width: 100%),
  caption: [(a) Prop. 4 pass rate of non-degenerate Prop. 2 passes against the normalised trace $t = a_1\/(2 sqrt(p))$
  (bins with ≥ 2,000 subgroups, 95% binomial intervals): the dependence flattens toward 1/2 as $p$ grows. (b) Fitted
  slope $beta_p$ on $a_1\/(2p+2)$ with 95% curve-clustered intervals, against $beta = 1$ and the persistence alternative
  $beta prop sqrt(p)$. Source `figures/plot.py`, data `data/scaling_p*.json`, `data/verdict.json`.])

= Findings
- *Measured, frozen rule.* $beta_p$ does not grow with $p$: $gamma = #f(v.gamma)$ [#f(v.gamma_ci.at(0)), #f(v.gamma_ci.at(1))],
  far from the 0.5 that persistence predicts. The Prop. 4 sign bias is $beta a_1\/(2p+2)$ with $beta$ bounded, so at a
  typical curve it is $O(1\/sqrt(p))$; Smith's fair-coin Hypothesis 3 is recovered as $p$ grows, and the small-$p$
  clustering of `full_success` within $L$-classes is a small-prime effect.
- *Measured.* $beta$ stays close to 1 from $p = 7$ to $p = 10007$ (#R.map(r => f(r.beta, n: 2)).join(", ") at
  $p = 101, 1009, 10007$): the pass rate is close to $\#H(FF_p)\/(2p+2)$ throughout. At $p = 1009$ the interval
  [#f(R.at(1).beta_ci.at(0)), #f(R.at(1).beta_ci.at(1))] just excludes 1; the three primes do not show a trend in that
  deviation.
- *Measured.* Degenerate Prop. 2 passes (Prop. 4 fails on both twists) fall from #R.at(0).degenerate at $p = 101$ to
  #R.at(2).degenerate at $p = 10007$ per about $5 dot 10^5$ passes, consistent with a small-$p$ construction failure.
- *Interpretation (labelled, not tested).* $r = \#H(FF_p)\/(2p+2) = \#X_+(FF_p)\/(\#X_+(FF_p) + \#X_-(FF_p))$, where
  $X_plus.minus$ are the two candidate quotient quartics: the outcome behaves as if a rational point were chosen
  uniformly on $X_+ union X_-$. No mechanism in Smith's construction is identified here that would explain this.

= Validation
All three primes valid: intercept within $0.5 plus.minus 0.02$ (#R.map(r => f(r.alpha, n: 4)).join(", ")); native point
count equal to `smith_prop24::h_point_count` on the first 1,000 squarefree curves (#R.map(r => str(r.controls.point_count_mismatches)).join(", ") mismatches);
no curve with more Prop. 4 passes than Prop. 2 passes. Standard errors are clustered by curve (CR1).

= Scope
Three primes up to $10007$, random degree-8 models. "Shrinks like $1\/sqrt(p)$" is the trend over this range, not a
theorem. The residual quarter of the within-class clustering at $p = 7, 11$ stays open; a post hoc check (protocol,
"Also recorded") found that the number of Prop. 2 passes adds nothing and $N_3$ at most a few points beyond a random
label. The mechanism behind the $\#X$-proportional rule is not tested.

= Reproduce
```
cargo build --release --example smith_sign_scaling
B=target/release/examples/smith_sign_scaling; D=research/smith_sign_scaling_20261016/data
for p in 101 1009 10007; do $B run --p $p --n 1000000 --out $D; done   # 50 s, 82 s, 225 s on 2 vCPU
$B verdict --dir $D 101 1009 10007
python research/smith_sign_scaling_20261016/figures/plot.py             # layout only
```
