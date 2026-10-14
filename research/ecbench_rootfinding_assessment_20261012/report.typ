#set page(paper: "us-letter", margin: (x: 0.72in, y: 0.72in))
#set text(size: 10pt)
#set par(justify: true)
#set heading(numbering: none)

#align(center)[
  #text(size: 17pt, weight: "bold")[Point Decomposition by Root Finding]
  #linebreak()
  #text(size: 12pt)[Source Assessment and Ecbench Measurement]
  #linebreak()
  12 October 2026
]

*Scope.* This record assesses the supplied 13-page manuscript and measures its direct-subtraction and pair-table strategies inside ecbench on one 20-bit prime curve and one 17-bit Koblitz curve. All 256 session executions verified and passed local replay. The plotted metric is counted group-addition equivalent work normalized by $sqrt(r)$; wall observations are L0. The manuscript's oracle ratios, composed larger-field totals, and these new counts have distinct scopes.

== Algorithm and source audit

For two summands, fixing a factor-base point makes the third summation-polynomial equation a quadratic with abscissa roots for $R + F_i$ and $R - F_i$. Point subtraction followed by a membership lookup therefore replaces a polynomial-system solve for that query. The prime implementation batches inversions; the binary implementation builds a pair-sum table and probes it for two to four summands. The table has quadratic setup and memory. Witness re-addition verifies soundness. Completeness requires the stated sign and factor-base convention; the source's binary code uses full point keys and both signs.

The manuscript reports 383–1535× faster binary two-summand queries against matrix F4 at $n = 9, 13$, excluding table setup. Its 3021–5802× NIST P-curve query ratios use a six-point base and a timed `buchberger` call, despite the program's printed F4 label. The directly timed 16-bit prime complete solve reports 100×; 20- and 24-bit totals are composed estimates. The manuscript is a draft whose novelty audit and stronger F4 cross-check remain open. The new measurements below use ecbench's generic arithmetic, not the paper's optimized prime kernel.

A separate committed native PDP root-finder study (`3745d06911c8`) checked the same mechanisms on polynomial-grid workloads: 707/707 prime brute-force agreements and binary corpus digest agreements, with growth at the generic pair-table or enumeration floors. That study did not run an ecbench complete solve. Its oracle ratios and the pipeline ratios below have different units and are kept separate.

== Fresh ecbench measurements

Eight public targets per curve, three measured rounds and one warmup, matched seeds and factor bases, and full-rank checked linear algebra produced 24 verified measurements per arm. The frozen source commit is `4ffb2c91600a`; sealed sessions are `ECBS1h40838be1a86e` and `ECBS1ha82cdeaf82db`. All 128 executions in each session passed full local replay with zero audit problems. A dagger marks unpriced native work.

#figure(
  image("mean_s.svg", width: 100%),
  caption: [Cold one-target counted cost. Each arm has 24 verified measured runs.]
)

#pagebreak()

On the prime curve, pair table / direct subtraction has paired counted-$S$ quotient 0.018861 (95% cluster interval 0.018547–0.019242); on the Koblitz curve it is 0.057047 (0.041159–0.092069). These are 24-pair local count comparisons. The Koblitz table spends 3537 GAE in setup, about 68% of its cold counted cost; its relation phase uses 1307 GAE versus 91030 for subtraction. The prime relation phases use 19659 versus 1142460 GAE, but prime native operations are unpriced. Matched rho has smaller counted $S$ on both curves. Rho itself has unpriced operations, so these cross-method ratios are bounded rather than complete-cost conclusions.

#figure(
  image("pipeline.svg", width: 100%),
  caption: [Phase map. The replacement acts at decomposition; setup is charged in a cold solve.]
)

== Phase accounting

#text(size: 8.5pt)[
  #table(
    columns: (1.4fr, 0.7fr, 0.7fr, 1fr, 0.7fr, 0.7fr, 0.8fr),
    inset: 3pt,
    align: (left, right, right, right, right, right, right),
    table.header([Arm], [Base], [Setup], [Relations], [LA], [Verify], [Trials]),
    [Prime subtraction], [0†], [0], [1,142,460], [0†], [849], [17,869],
    [Prime pair table], [0†], [1,056], [19,659], [0†], [849], [17,869],
    [Koblitz subtraction], [54], [0], [91,030], [14], [302], [321],
    [Koblitz pair table], [54], [3,537], [1,307], [14], [302], [322],
  )
]

Mean ecbench GAE by phase; † denotes native work present but unpriced. Relation trials are counts, not GAE. Setup is charged once for each cold solve.

Measured child-process peak RSS was 30,992 KiB for both prime IC arms; 31,328 KiB for Koblitz subtraction and 31,360 KiB for its pair table. These are whole-process maxima and do not isolate table storage at these base sizes.

== Interpretation and limits

The algebraic replacement removes the measured polynomial-oracle overhead. A precomputed pair table can strongly reduce counted relation cost for repeated queries, with setup charged in the one-target solve. These two small-curve results do not establish a larger-field scaling law, a qualified wall-time advantage, or transfer of the manuscript's NIST query ratios into a complete ECDLP cost. Independent-host replay, priced native operations, a matched production-size baseline, and the manuscript's open novelty and F4 checks remain necessary for stronger claims. The full report, exact specs, sealed records, replay receipts, and source-derived SVGs accompany this PDF.

== Reproducibility binding

The supplied PDF has SHA-256 `743d3d3e7e3aab6f980813c0944c1d2cbacc906e075e7283de0abc7b923b7c01`, identical to the PDF in the paper branch at commit `5bb2f20edb`. The ecbench release executable has SHA-256 `be956553a5a2cd1227e8f1add8dd24e11dffba5068c799234b55e44ffc3a1b61`. The all-record replay receipts and their hashes are listed in `REPORT.md`; the two specs and sealed sessions accompany this document. No canonical frontier graph was changed because the local wall level is L0, the rho baselines contain unpriced native operations, and an independent host has not replayed these sessions.
