#set page(paper: "a4", margin: 1.6cm)
#set text(size: 9.5pt)
#set par(justify: true)

= Pair-46 factor-base campaign (2026-10-15)

*Study:* `research/pair46_factor_bases_20261015/`. *Branch:* `claude/pair-46-elliptic-curves-71ny4x`. Pair 46 has $p = 5t^2+19t+19$ with $t = 34489701$, a 53-bit prime $N = p+1-t$, a CM($-19$) crater with $j=-96^3$, and a floor curve. There were 14 items, each with a preregistered protocol, a run, a results report and an independent red-team review. Every verdict below is the red team's corrected verdict. 12 of the 13 items that ran were upheld, item 14 was refuted (`exact_obstruction` becomes `empirical_negative`), and item 11 was blocked.

#figure(
  image("fig_campaign_summary.svg", width: 100%),
  caption: [Headline ratios against each item's matched reference, with 95% intervals as recorded in the item's RESULTS.md (no recomputation). Item 01 is a cost ratio, so values above 1 are worse. The rest are yield, coverage or density ratios. Data: `summary_data.json`. Source: `plot_summary.py`.],
)

== Summary

#table(
  columns: (auto, 1.6fr, 1.3fr, 3.2fr),
  inset: 4pt,
  stroke: 0.4pt + luma(170),
  table.header([*\#*], [*item*], [*verdict*], [*headline (measured / exact)*]),
  [01], [endo-chain], [negative], [ψ-chain S = 1.47 [1.28, 1.69] against ±BSGS 0.955 [0.86, 1.05] (ratio 1.54 [1.30, 1.82]). Exact: a known-head chain is BSGS with $c_E approx 1$.],
  [02], [deuring-lift], [exact obstruction], [The lift is 361a1 with rank 1, so lifting a target is the DLP. Twist bases give 0 Q-relations. λ1/GH difference +0.14 [−0.18, +0.26] against random.],
  [03], [floor-pullback], [exact obstruction], [φ is an isomorphism, so the yield is identical to the crater's. Toy ratio 0.9978 to 1.0023, |z| ≤ 1.40; the positive control had |z| ≥ 562.],
  [04], [family-sweep], [null], [M1 0.984 [0.964, 1.005], M3 1.0011 [0.998, 1.004], 20 to 30 bits. All positive controls had p ≤ 0.0005.],
  [05], [qr-window], [constant only], [Density gain 2.00004 [1.99967, 2.00042]. The relation ratio is 1 exactly (8000/8000).],
  [06], [pkm-rational-map], [negative], [Cost per system about $p^(2.0)$ (local toy slope), against $p^(0.52)$ for rho. About 4e6× rho at n1 = 5. Yield equals the flat set exactly.],
  [07], [torsion-division], [exact obstruction], [No $psi_m$ root on E. Preimage bases change coverage by exactly 0. Size/degree ρ ≤ 0.25.],
  [08], [learned], [exact obstruction (counting); empirical part inconclusive], [Learned selectors are within ±0.5% of random. $sum_R c(R) = |S|^m$. The advantage is at most about μ/(1−e^(−μ)).],
  [09], [cm-norm], [exact obstruction], [$F_X + F_X subset.eq F_(4X)$. Coverage is 0.0159 [0.0156, 0.0161]× random at B = 256.],
  [10], [bounded-degree], [exact obstruction], [|F|/p ≥ 0.167 for d ≤ 8. |F| ≤ √p needs d ≥ √p/5.],
  [11], [ff-xedni], [blocked], [Not run. B2 and function-field xedni remain open.],
  [12], [model-small-x], [exact obstruction], [No Montgomery, Edwards or Hessian model over $FF_p$ for pair 46. Toy tests within ±1% of random.],
  [13], [xy-windows], [null], [Matched-size y-windows give D/D(X) 0.986 to 0.995, CIs within [0.956, 1.026], at 4 to 64× the scan price.],
  [14], [dor-verdict], [*empirical negative* (refuted)], [Binary toy model. The abandon rule fires for generic bases (0.86 to 1.36 per doubling) but never for the standard window. E1 fails for the degree-2 Macaulay rank on 937 of 4638 targets. Does not transfer to pair 46.],
)

== Conclusion, as strong as the evidence

- *No item found a factor-base choice that gives index calculus an advantage on pair 46 or on its toy analogues.*
- *The exact obstructions are scoped* to the tested families and metrics. They do not establish that no factor base can help. For example, bases chosen to make solving cheaper (item 08) and compact high-genus descriptions (item 10) are not covered.
- *The empirical findings are toy-scale* (12 to 30 bits, m = 2). The nulls have a resolution of about ±2 to 3%.
- *The cost results 01 and 06 are constants or local slopes*, not asymptotic exponents.
- *Not tested:* m ≥ 3, except through exact identities; relation collection with linear algebra at pair-46 size; F4/msolve; and item 11.

== Positive-control status (item 05)

The positive control was detected. The predicted 2× density gain was measured at 2.00004 [1.99967, 2.00042], and the negative controls gave 1. The agreement of the main checks holds by construction. The independent validation is the E1 full-field counts and the A2/A3 halving controls, which keep about 1/4 of the m = 2 decompositions. Every other item that ran had its own positive control detected. The weak points are 08 (the random2 null failed at 16 bits) and 06 (the preregistered sugar-degree match failed; the amended lcm-degree match passed).

== Blocked and cut

- *11, blocked:* no artifacts.
- *Cut:*
  - 03: no pair-46-size isogeny.
  - 04: no floor members, and no M1/M2 at 53 bits.
  - 06: chained arms stopped at n1 = 5, and there was no F4.
  - 07: no Gröbner invariants.
  - 14: m = 3, sym-S4 and n = 131 were not ported.
- *Hosts:* timings come from a shared host at isolation level L0. No ecbench sessions were used.

Each item's own `report.pdf`, `RESULTS.md` and `REDTEAM.md` hold its full data, intervals, deviations and reproduce commands.
