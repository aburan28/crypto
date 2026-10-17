// Results report: branch-free selections in production Fp64 (fp_backend and
// prime_fast). Tables read the run records at build time; compile with the
// study directory as root.
#set page(paper: "a4", margin: 1.7cm, numbering: "1")
#set text(size: 9.5pt)
#set heading(numbering: "1.")
#show table: set text(size: 7.6pt)

#let R = "runs/"
#let F = "reports/run-20261016T2319Z/"
#let r(x, d) = str(calc.round(float(x), digits: d))
#let csum = json(R + "compare/summary.json")
#let pair = json(R + "cm-pair/claims.json")
#let cmp = csv(R + "compare/compare.csv", row-type: dictionary)
#let dis = csv(R + "disasm/counts.csv", row-type: dictionary)
#let pc = read(R + "protocol_commit.txt").trim()
#let claims = csum.claims + pair
#let c(id) = claims.find(x => x.id == id)
#let arm(a) = cmp.find(x => x.arm == a)
#let classes = ("candidate faster", "candidate slower", "no detected effect")
#let count(g, k) = {
  let hit = csum.classes.find(x => x.group == g and x.class == k)
  if hit == none { 0 } else { hit.arms }
}

#align(center)[
  #text(size: 14pt, weight: "bold")[Branch-free selections in production `Fp64`] \
  #text(size: 10pt)[`fp_backend::Fp64` and `prime_fast::Fp64` with
  `select_unpredictable`: results report, 2026-10-16]
]

*Study:* `research/fp64_branchfree_20261016/`; protocol `PROTOCOL.md` frozen at
#raw(pc.slice(0, 9)) before the first run. *Run report:*
`reports/run-20261016T2319Z/run_report.pdf`. *Origin:* `hwcost` round 2
(`research/hwcost_round2_20261016/`, addendum 1). *Code:*
`src/cryptanalysis/fp_backend.rs`, `src/cryptanalysis/prime_fast.rs`.

*Scope.* An engineering change to two field implementations, measured on
`hwcost`'s rho kernel and on the production prime-field rho
(`prime_fast::rho_parallel`). It changes no algorithm and no counted cost
(every result is bit for bit the previous one), so it is not an
index-calculus improvement, adds no gap-to-rho row and ranks no ECDLP
methods (AGENTS.md success criterion, §5a).

#let k(f) = cmp.filter(x => x.arm.starts-with("kernel." + f + ".") and int(x.arm.split(".").at(2)) >= 4096).map(x => (float(x.r1), float(x.r2))).flatten()
#let rng(xs) = r(calc.min(..xs), 2) + "–" + r(calc.max(..xs), 2)
#let pf = arm("pfrho.44")
#let rw = csv(R + "disasm/rho_worker_counts.csv", row-type: dictionary).filter(x => not x.function.contains("closure"))
#let vals = ("base-1", "cand-1", "base-2", "cand-2").map(n => json(R + n + "/validation.json"))
#let v2(i, f) = vals.at(i).claims.find(x => x.id == "V2." + f).value

= Findings

+ *`fp_backend::Fp64`.* The rho kernel (`hwcost`, field 2#super[63] − 25) takes
  #rng(k("fp64"))× of the baseline's time at every B ≥ 4096 in both pairs of
  runs (P1 #c("P1").verdict, worst #r(c("P1").value, 3) against 0.85). The
  candidate's kernel step compiles to `Fp64Cm`'s machine code (D2: 303 → 285
  instructions, 32 → 20 conditional jumps, 4 → 12 `cmov`). It runs within
  #r((c("P3").value - 1) * 100, 1)% of `Fp64Cm` in the same runs (P3
  #c("P3").verdict).
+ *`prime_fast::Fp64`.* The same kernel through `prime_fast`'s arithmetic
  (`pf64`, 2#super[62] − 57) takes #rng(k("pf64"))× at every B ≥ 4096 (P1pf
  #c("P1pf").verdict, worst #r(c("P1pf").value, 3)).
+ *The production rho.* `prime_fast::rho_parallel` on a 44-bit prime-order
  curve, one thread, 7 known-answer solves per run: the median goes from
  #r(pf.a1, 1) to #r(pf.b1, 1) ns per step in the first pair and from
  #r(pf.a2, 1) to #r(pf.b2, 1) in the second (ratios #r(pf.r1, 3) and
  #r(pf.r2, 3), drift #r(float(pf.drift) * 100, 1)%; P4 #c("P4").verdict
  against 0.95). That is #r(1 / float(pf.r2), 2)–#r(1 / float(pf.r1), 2)×
  as many steps per second. Its worker loop `prime_fast::rho_worker` has
  #rw.at(0).jcc → #rw.at(1).jcc conditional jumps and #rw.at(0).cmov →
  #rw.at(1).cmov `cmov`.
+ *Guards and controls.* Dependent multiply and squaring chains move by at
  most #r((calc.max(c("P2a").value, c("P2b").value, c("P2c").value, c("P2d").value) - 1) * 100, 1)%
  (P2a–d PASS against +10%). The unchanged `Fp64Cm` and `FpWide<1>` kernels
  stay within #r(c("C1a").value * 100, 1)% and #r(c("C1b").value * 100, 1)%
  (C1 PASS against ±10%). Memory bandwidth (20 arms) and clock (2 arms)
  show no detected effect. Inversion arms show no detected effect either
  (the binary extended Euclid is unchanged; only its final product is not).
+ *One cost, stated.* The `add.latency` arms time a dependent chain
  `a ← a − b` with a fixed `b`, whose wrap-around repeats in a pattern a
  predictor learns. On that chain the select costs latency where the branch
  was free:
  - `fp64`: #r(arm("fp64.add.latency").a1, 2)/#r(arm("fp64.add.latency").a2, 2)
    → #r(arm("fp64.add.latency").b1, 2)/#r(arm("fp64.add.latency").b2, 2) ns
    (ratios #r(arm("fp64.add.latency").r1, 2), #r(arm("fp64.add.latency").r2, 2)),
    classed no detected effect only because the baseline drifted
    #r(float(arm("fp64.add.latency").drift) * 100, 0)% between its runs;
  - `pf64`: ratios #r(arm("pf64.add.latency").r1, 2) and
    #r(arm("pf64.add.latency").r2, 2) (candidate slower).

  The same arm also moved for code that did not change (`fp64cm`
  #r(arm("fp64cm.add.latency").r1, 2), `w1` #r(arm("w1.add.latency").r1, 2)).
  This 1-ns chain is therefore sensitive to code placement between
  binaries, and part of the move may not be the change. Rho never runs such
  a chain: its subtractions take independent, effectively random operands.

*The cost model agrees.* The a priori op-count model (V2, no fitted
parameter) missed the baseline `fp64` kernel by #r(v2(0, "fp64") * 100, 0)% and
#r(v2(2, "fp64") * 100, 0)%, and `pf64` by #r(v2(0, "pf64") * 100, 0)% and
#r(v2(2, "pf64") * 100, 0)% (round 1's `fp64` anomaly, which `pf64` shares). On the candidate it predicts
`fp64` within #r(v2(1, "fp64") * 100, 0)% and #r(v2(3, "fp64") * 100, 0)%, and
`pf64` within #r(v2(1, "pf64") * 100, 0)% and #r(v2(3, "pf64") * 100, 0)%
(descriptive; `validation.json` per run). Counted operations now explain both
production kernels.

*Reading.* Both production fields now do what round 2's `Fp64Cm` did, with
the same results bit for bit, and the gain carries into the production rho:
about 1.4× more steps per second for `prime_fast::rho_parallel` on this host.
Every caller inherits it, including the F4 engine's row reductions,
`two_sylow_search`, `Fp52`'s scalar operations and `prime_fast`'s
index-calculus relation search, none of which was timed here.

Counted costs (GAE, field operations) do not change, so no index-calculus
comparison in counted units moves. Wall-time comparisons against the
`prime_fast` rho now face a faster rho.

*Limits.*
- Host L0: one shared 4-vCPU cloud VM with no PMU, so branch misses are
  inferred, not counted.
- Baseline-ISA builds only. A `-C target-cpu=native` build already turned
  some of these branches into `cmov` (round 2).
- `pfrho.44` is one curve size at one thread, with its set-up spread over a
  random number of steps per solve.
- `step_pf64` was never called, so it was linked in neither binary. It has
  been removed since, and `rho_worker` stands in for it in D2.

= Claims

#table(
  columns: 5,
  align: (left, left, right, right, left),
  [id], [statement], [value], [threshold], [verdict],
  ..claims.map(x => ([#x.id], [#x.statement], [#r(x.value, 4)], [#x.threshold], [#x.verdict])).flatten()
)

Arm classes, candidate against baseline (an effect is called only when both
ratios clear max(2·drift, 2%) on the same side):

#table(
  columns: 4,
  align: (left, right, right, right),
  [arm group], [candidate faster], [candidate slower], [no detected effect],
  ..csum.classes.map(x => x.group).dedup().map(g => ([#g], ..classes.map(k => [#count(g, k)]))).flatten()
)

#figure(image(F + "fig_kernel.svg", width: 100%), caption: [Rho kernel time per field,
baseline (blue) and candidate (red) binaries, runs A1 B1 A2 B2.])

#figure(image(F + "fig_ratios.svg", width: 100%), caption: [Candidate / baseline ratio
for every arm.])

#figure(image(F + "fig_pfrho.svg", width: 65%), caption: [The production rho
`prime_fast::rho_parallel` (44-bit curve, one thread): every sample.])

#figure(image(F + "fig_disasm.svg", width: 75%), caption: [Machine code of
`hwcost`'s non-inlined kernel step and of `prime_fast::rho_worker`, the
production rho's loop.])
