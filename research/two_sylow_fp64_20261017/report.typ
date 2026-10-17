// Results report: two_sylow_search on fp_backend::Fp64, before and after the
// branch-free selections. Tables read the run records; compile with the
// study directory as root.
#set page(paper: "a4", margin: 1.7cm, numbering: "1")
#set text(size: 9.5pt)
#set heading(numbering: "1.")
#show table: set text(size: 7.8pt)

#let R = "runs/"
#let F = "reports/run-20261017T0018Z/"
#let r(x, d) = str(calc.round(float(x), digits: d))
#let csum = json(R + "compare/summary.json")
#let cmp = csv(R + "compare/compare.csv", row-type: dictionary)
#let c(id) = csum.claims.find(x => x.id == id)
#let arm(a) = cmp.find(x => x.arm == a)
#let tsy = cmp.filter(x => x.arm.starts-with("tsy."))
#let r1s = tsy.map(x => float(x.r1))
#let r2s = tsy.map(x => float(x.r2))
#let w1 = cmp.filter(x => x.arm.starts-with("kernel.w1.") and int(x.arm.split(".").at(2)) >= 4096)
#let pc = read(R + "protocol_commit.txt").trim()

#align(center)[
  #text(size: 14pt, weight: "bold")[`two_sylow_search` on `fp_backend::Fp64`] \
  #text(size: 10pt)[Before and after the branch-free selections: results report, 2026-10-17]
]

*Study:* `research/two_sylow_fp64_20261017/`; protocol `PROTOCOL.md` frozen at
#raw(pc.slice(0, 9)) before the first run. *Run report:*
`reports/run-20261017T0018Z/run_report.pdf`. *Follows:*
`research/fp64_branchfree_20261016/` (PR \#2034), whose "Next" item 1 asked
whether a latency-bound user of `fp_backend::Fp64` relies on the
predictable `sub` chains that the change slowed.

*Scope.* Wall time of one production module's `Fp64` instantiation. No
counted cost changes; not an index-calculus result; no gap-to-rho row.

= Findings

+ *No detected gain, no slowdown beyond the bound.* Across the six arms
  (two families × three search arms) the candidate/baseline ratio is
  #r(calc.min(..r1s), 2)–#r(calc.max(..r1s), 2) in the first pair of runs
  and #r(calc.min(..r2s), 2)–#r(calc.max(..r2s), 2) in the second (S1
  #c("S1").verdict, worst #r(c("S1").value, 3) against 1.10). Five arms show
  no detected effect; `tsy.z2xz2.fast-chi` is classed slower
  (#r(arm("tsy.z2xz2.fast-chi").r1, 3), #r(arm("tsy.z2xz2.fast-chi").r2, 3)).
+ *The second pair reads high for unchanged code too.* Every `tsy` arm's
  second ratio is above its first, and the candidate's second run (B2) is
  the noisiest of the four (Figure 1). In the same pair, the unchanged
  `FpWide<1>` kernel at B ≥ 4096 reads #r(calc.min(..w1.map(x => float(x.r2))), 3)–#r(calc.max(..w1.map(x => float(x.r2))), 3)
  (C1 #c("C1").verdict, max deviation #r(c("C1").value, 3)). Four of the
  20 memory-read arms, whose code is identical in both binaries, were also
  classed slower. Part of the second pair's 2–8% is host drift, not the
  change.
+ *The comparison sees what it should.* The positive control, the `fp64`
  rho kernel, is #r(c("K1").value, 3)× or less of baseline at every
  B ≥ 4096 (K1 #c("K1").verdict), replicating PR \#2034. Gate G: both
  binaries computed the same search (identical depth sums, point counts and
  per-arm checksums over 2,000 draws per family, in all four runs).
+ *Why no gain.* `Fp64` has no specialised `HalvingField` methods. A draw is
  dominated by `pow` (square roots and `(p − 3)/4` powers, dependent
  `mul`/`sqr` chains) and `jacobi_u64`, neither of which the change
  touches. The dependent `mul` and `sqr` latency arms measured in the same
  runs are unchanged (`fp64.mul.latency` #r(arm("fp64.mul.latency").r1, 3),
  #r(arm("fp64.mul.latency").r2, 3); `fp64.sqr.latency`
  #r(arm("fp64.sqr.latency").r1, 3), #r(arm("fp64.sqr.latency").r2, 3)).
  The dependent `sub` chain is again slower (`fp64.add.latency`
  #r(arm("fp64.add.latency").r1, 2), #r(arm("fp64.add.latency").r2, 2)), as
  in PR \#2034. The search's handful of `add`/`sub` per level is not enough
  for that to move a draw beyond the noise.

*Reading.* PR \#2034's change is neutral for `two_sylow_search` on
`fp_backend::Fp64` within this host's resolution: no arm slowed by 10%, and
the first pair (the less noisy one) is within ±2.3%. Production runs at P-192
use `FpP192`, which the change does not touch.

*Limits.* Host L0 without a PMU; one prime (2#super[63] − 25, the largest
prime ≡ 3 mod 4 below 2#super[63]); record depth 6; baseline-ISA builds only.

= Claims

#table(
  columns: 5,
  align: (left, left, right, right, left),
  [id], [statement], [value], [threshold], [verdict],
  ..csum.claims.map(x => ([#x.id], [#x.statement], [#r(x.value, 4)], [#x.threshold], [#x.verdict])).flatten(),
  [G], [gate lines identical in all four runs], [—], [—], [#read(R + "compare/gate.txt").trim().split(":").at(0).replace("G ", "")],
)

#figure(image(F + "fig_tsy.svg", width: 100%), caption: [Every sample of the six
`two_sylow_search` arms, µs per draw, in the four runs.])

#figure(image(F + "fig_ratios.svg", width: 100%), caption: [Candidate/baseline
ratios of the `tsy` arms and of the kernel controls (positive: `fp64`, negative:
`w1`) at B ≥ 4096.])
