#set page(paper: "a4", margin: 18mm)
#set text(size: 10pt)
#set par(leading: 0.65em)

= Torsion-coordinate relations over a quadratic extension

*Cold ECDLP cost on four frozen public prime-order subgroups, 2026-10-14*

The measured group is $E(F_(p^2))[r]$ for $E: y^2=x^3+1$ and subgroup orders $r=19,61,211,811$. The rational point $T=(0,1)$ has order three. The candidate coordinate $chi=(y+1)/(y-1)$ has a public order-three Mobius action under translation by $T$. The experiment prices a projected pair table, an exact $chi$ orbit index on that table, and a degree-three isogeny control against the same-target negation rho, A/A rho, and charged six-class rho.

The table below gives *cold current-thread CPU nanoseconds divided by a measured source-curve addition and by* $sqrt(r)$. It includes setup, collection, algebra, recovery, and scalar verification. This host-dependent unit is denoted $S_("CPU/GAE")$. The displayed ratio is the median of paired, same-round complete costs; it is not a phase ratio. Three rounds were measured at each size.

#table(
  columns: (auto, auto, auto, auto, auto, auto),
  table.header([Order $r$], [Projected], [$chi$ index], [Isogeny], [Negation rho], [$chi$/rho]),
  [19], [164.19], [195.05], [172.81], [19.93], [9.74],
  [61], [138.43], [156.62], [140.38], [12.74], [12.59],
  [211], [159.01], [190.65], [160.28], [10.72], [19.47],
  [811], [202.56], [284.96], [217.45], [7.30], [51.89],
)

*Result.* The $chi$ index, projected table, and isogeny control had the same relation outcomes, rank, and recovered scalar in all twelve paired size-by-round cases. The index added CPU work. The 84 measured arm rows include 75 independent-arithmetic scalar replays; nine raw-pair rows reached the frozen cap. The three-round spread is substantial, and the host was shared (L0). The graph on the next page shows all seven arms and the complete ranges.

The exact pair census found projected supports $18/18$, $58/60$, $196/210$, and $742/810$. Torsion shifts added two *raw* subgroup support points at $p=787$ on the frozen base, yet created no new projected support and did not include the frozen target. This precise local signal does not establish a faster DLP. The unimplemented $chi$ polynomial oracle remains the next discriminating experiment.

#pagebreak()

== Complete arm ranges

#image("GRAPH.svg", width: 100%)

Each dot is the median cold local $S_("CPU/GAE")$; whiskers span the three measured rounds. Hollow raw-pair markers show capped rows and have no scalar answer.

#pagebreak()

== What was charged and verified

The cold interval starts before line lifting and cofactor projection, then includes base construction, pair table, unsuccessful queries, exact relation checks, modular elimination, scalar recovery, and internal verification. Isogeny map work and six-class eigenvalue discovery are included in their arms. The second native field/curve implementation checked every recovered $[k]G=Q$ against the frozen fixture. The CPU conversion calibrations were 169.76625, 147.45333, 162.66167 and 156.87292 nanoseconds per source addition. The generic two-class count floor fixed in the protocol is $sqrt(pi/4)=0.886226925$ counted group operations per $sqrt(r)$; a CPU-converted equivalent is only a local comparison to that floor.

The base-field embedding degrees $"ord"_r(p)$ are 9, 10, 105 and 810; relative to $q=p^2$, the degrees are 9, 5, 105 and 405. This only screens a possible pairing transfer into a multiplicative extension field. A pairing, its transfer cost and a multiplicative-field DLP were not constructed or measured. Per-arm peak RSS, isolated L2 timing, independent-host replay, direct Semaev resultant, and the nonfactorable $chi$ polynomial system remain open. Setup reuse count was one target.

The first timed output was retained despite a sandboxed system-time wrapper failure. The second included JSON serialization in its method interval. The published third version stops the clock before serialization; the earlier bytes remain in the study. The finite counts and CPU ratios are observations on four exact models, not a general curve theorem. The conditional torsion-support Lean theorem and the two polynomial identity certificates are linked from the prior proof map; no new general Lean claim is made here.

*Source and replay:* `REPORT.md`, `RUN_SPEC.json`, `FIXTURES.jsonl`, `RUN_V3.jsonl`, `REPLAY_INDEPENDENT_V4.json`, `SUMMARY_V3.jsonl`, `SUPPORT.jsonl`, and `ARTIFACTS.sha256` in this study directory. The complete 28-row table, hashes, per-phase counters, proof scope, and next gate are in `REPORT.md`.
