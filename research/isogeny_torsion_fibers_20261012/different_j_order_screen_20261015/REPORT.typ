#set page(paper: "a4", margin: 18mm)
#set text(size: 9.5pt)
= Bounded larger-order different-j isogeny screen
*2026-10-15 · native Rust · Darwin arm64, macOS 26.6 · public synthetic curves.* The protocol and exact source/binary hashes were committed before execution. All 20 eligible primes with 61 ≤ p ≤ 271 and p = 1 (mod 3) were emitted. Seventeen failed the frozen subgroup criterion, and three passed. The first passing p was selected as specified.

#table(
  columns: (0.8fr, 1.35fr, 0.7fr, 0.9fr, 1fr, 1fr),
  stroke: 0.5pt + rgb("cbd9de"), inset: 5pt,
  [*p*], [*Both orders over F(p²)*], [*r*], [*Cofactor*], [*Source j*], [*Dest. j*],
  [*211*], [*44,655*], [*229*], [*195*], [*77*], [*160*],
  [223], [50,175], [223], [225], [183], [222],
  [241], [58,563], [241], [243], [33], [164],
)

The selected source is y²=x³−21x+37 with kernel generator (3,1); the explicit degree-three destination is y²=x³−81x−243. Its 44,655 points factor as 3·5·13·229. Independent field and point formulas checked *192/192* sampled map witnesses. A separately frozen direct F(p²) quadratic-character census counted both curves at each accepted p and matched all three trace-derived orders. Different j values exclude an isomorphism between the two selected models.

#image("GRAPH_V1.svg", width: 100%)

The graph is a candidate census, not a DLP cost plot. It shows the largest prime factor of each group order and the r ≥ 211 selection threshold. The selected order-229 subgroup has 228 distinct nonidentity points, enough for the planned 64-target batch. A subsequent separately frozen full pointwise map/fiber audit passed; see FIBER_REPORT.pdf. Relation collection, scalar recovery, direct BSGS/rho comparison, counted S, isolation and independent-host performance replay remain pending.

*Evidence:* RAW_V1.jsonl SHA-256 31238914df366380; REPLAY_V1.json ac2c13d5496755d5; CENSUS_V1.jsonl 50f2b79d9811c72c. The source-linked RESULTS.md records full hashes, scripts, empty stderr receipts, checker scope, and the next test. This screen does not change the canonical ECDLP frontier.
