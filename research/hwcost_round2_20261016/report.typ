// Results report, hwcost round 2. Tables read the run records at build time;
// compile with the study directory as root.
#set page(paper: "a4", margin: 1.7cm, numbering: "1")
#set text(size: 9.5pt)
#set heading(numbering: "1.")
#show table: set text(size: 7.6pt)

#let R = "runs/"
#let F = "reports/run-20261016T1654Z/"
#let r(x, d) = str(calc.round(float(x), digits: d))
#let isa = ("isa-base-1", "isa-native-1", "isa-base-2", "isa-native-2")
#let val = isa.map(n => json(R + n + "/validation.json"))
#let cmv = ("cm-1", "cm-2").map(n => json(R + n + "/validation.json"))
#let tr = ("isa-native-1", "isa-native-2").map(n => json(R + n + "/transfer/validation.json"))
#let cmp = csv(R + "compare/compare.csv", row-type: dictionary)
#let csum = json(R + "compare/summary.json")
#let cmcl = json(R + "cm-pair/claims.json")
#let cmp2 = csv(R + "cm-pair/pairs.csv", row-type: dictionary)
#let f3 = csv(R + "f3/counts.csv", row-type: dictionary)
#let gcl = json(R + "gpu/claims.json")
#let gsum = read(R + "gpu/summary.txt").trim().split("\n")
#let ic = ("koblitz", "prime").map(d => (
  name: d,
  pairs: csv(R + "ic-" + d + "/pairs.csv", row-type: dictionary),
  fits: json(R + "ic-" + d + "/fits.json"),
  proj: csv(R + "ic-" + d + "/projections.csv", row-type: dictionary),
  gae: csv(R + "ic-" + d + "/projections_gae.csv", row-type: dictionary),
  sum: json(R + "ic-" + d + "/summary.json"),
))
#let cl(v, id) = v.claims.find(c => c.id == id)
#let gc(id) = gcl.find(c => c.id == id)
#let v(c) = [#r(c.value, 3) #c.verdict]
#let fields = ("fp64", "fp64bf", "w1", "w2", "w3", "w4")
// Display only: a positive number as m.mm × 10^e.
#let sci(x) = {
  let e = calc.floor(calc.log(float(x)))
  [#r(float(x) / calc.pow(10.0, e), 2)#sym.times;10#super[#e]]
}

#align(center)[
  #text(size: 14pt, weight: "bold")[hwcost round 2] \
  #text(size: 10pt)[The ISA effect, the `fp64` branch lead, a GPU profile, and index
  calculus in time and dollars: results report, 2026-10-16]
]

*Study:* `research/hwcost_round2_20261016/`; protocol `PROTOCOL.md` frozen at
`597502b4d` before every run and analysis, addendum 1 frozen at `4e5994fc9`
before its runs. *Run report:* `reports/run-20261016T1654Z/run_report.pdf`
(every run, its host, gates and raw values). *Round 1:*
`research/hwcost_model_20261016/`. *Code:* `src/cryptanalysis/hwcost/`,
`src/bin/hwcost.rs` (Rust; Python only lays out graphs).

*Scope.* This round measures and prices; it changes no index-calculus
pipeline and solves no discrete logarithm. It is *not* an index-calculus
improvement and adds no gap-to-rho row (AGENTS.md success criterion). Part D
compares IC with rho from frozen records and ranks two arms only where the
wall-time order and the order at every §5a weight (0, 0.1, 1 GAE per unpriced
native operation) agree (AGENTS.md §5a).

= Findings

+ *ISA (part A).* A `-C target-cpu=native` build runs the `fp64` and
  `fp64bf` rho kernels 1.6–2.6× faster at every batch size (native/baseline
  0.38–0.63 in both pairs of runs; 0.38–0.47 at B ≥ 256); the multi-limb `FpWide` kernels show no
  detected effect at 37 of 40 batch sizes. The model calibrated on the native
  build predicts it (I2 passes for every field), but a residual fitted on the
  baseline build does not transfer to the native one for `fp64`/`fp64bf`
  (I3 fails, error 1.3–1.5) while it does for `w1`–`w4` (≤ 0.10).
+ *Branches (part B, addendum 1).* The pre-registered test F1 was invalid:
  LLVM compiled `Fp64Bf`'s masks back into the same branchy code as `Fp64`
  (F3: identical 303 instructions, 32 conditional jumps). With the selections
  written as `select_unpredictable` (`Fp64Cm`: 20 jumps, 12 `cmov`), the
  kernel takes *#r(cmcl.at(0).value, 2)× or less* of `fp64`'s time at every
  B ≥ 4096 in both runs (F1′ PASS), and the a priori model, with no fitted
  parameter, predicts it within #r(cl(cmv.at(0), "V2.fp64cm").value * 100, 0)% and
  #r(cl(cmv.at(1), "V2.fp64cm").value * 100, 0)% (F2′ PASS) where it misses `fp64` by
  #r(cl(cmv.at(0), "V2.fp64").value * 100, 0)%. Round 1's `fp64` anomaly is the
  data-dependent branches in `Fp64::sub`/`redc`: they cost half or more of
  the kernel's time on random operands. Causal evidence by substitution; branch
  misses themselves could not be counted (no PMU on this host).
+ *GPU (part C).* With SMs as cores and the busiest pipe binding, the a priori
  profile predicts all 46 recorded ECC2K-130 builds within a median
  #r(gc("G1").value * 100, 1)% (G1), never above the measured rate (no build
  with η > 1); one efficiency η per GPU brings it to #r(gc("G2").value * 100, 1)% (G2);
  carried from the RTX PRO 6000 to the B200, #r(gc("G3").value * 100, 1)% (G3). The
  cost arithmetic reproduces the paper's ECC2K-130 figure (\$123,749 vs \$124k).
+ *Index calculus (part D).* On every measured curve where §5a admits a
  rank, IC costs more wall time than the matched rho; on no curve does IC rank
  below rho. The prime band *understates* IC's wall cost from log2 r = 32 up
  (measured ratio above the band's whole range, 1.1–2.8×); the Koblitz band
  holds the measured ratio on its 3 largest curves. Both projections, the
  pre-registered wall-time fit and the post hoc GAE fit, rank IC > rho at
  every target from $2^(48)$ to $2^(128)$ and at every weight; in the post hoc
  projection IC costs 29 to 3.7 × $10^(4)$ times rho's time (Koblitz) and
  290 to 3.8 × $10^(7)$ times (prime).

= Part A: the ISA effect

Four full `measure` runs alternated baseline (A) and native (B) builds,
A1 B1 A2 B2, and each arm was classified only when both ratios $r_k = B_k\/A_k$
cleared $max(2 d, 2%)$ on the same side, $d$ the A1/A2 drift (claim I1).

#let classes = ("candidate faster", "candidate slower", "no detected effect")
#let count(g, k) = {
  let hit = csum.classes.find(c => c.group == g and c.class == k)
  if hit == none { 0 } else { hit.arms }
}
#table(
  columns: 4,
  align: (left, right, right, right),
  [arm group], [native faster], [native slower], [no detected effect],
  ..csum.classes.map(c => c.group).dedup().map(g => ([#g], ..classes.map(k => [#count(g, k)]))).flatten()
)

The kernel arms called "faster" are all 20 `fp64` and `fp64bf` arms; among
the 40 `FpWide` kernel arms, `w3` at B = 16 and `w4` at B = 64 are faster by
4–7% and `w1` at B = 65536 slower by 3–5%; every other one shows no detected
effect. Of the 48 field-operation arms, 35 show no detected effect, 10 are
slower on the native build (by 2% to 59%, the largest `fp64bf.add` latency)
and 3 faster (by 2–12%): isolated operations got no faster, so the kernel
speedup comes from the kernel's code generation (fewer branches, `mulx`,
vector moves; table in part B). The memory-read arms are faster on the native
build in 9 of 20 cases (ratios down to 0.43), so the read bandwidth this
harness measures, and with it round 1's roofline constants, depends on the
build.

#table(
  columns: 7,
  align: (left, right, right, right, right, right, right),
  [field], [V1 isa-base-1], [V1 isa-base-2], [I2 (V1) isa-native-1], [I2 (V1) isa-native-2], [I3 isa-native-1], [I3 isa-native-2],
  ..fields.map(f => ([#f], v(cl(val.at(0), "V1." + f)), v(cl(val.at(2), "V1." + f)), v(cl(val.at(1), "V1." + f)), v(cl(val.at(3), "V1." + f)), v(cl(tr.at(0), "V1." + f)), v(cl(tr.at(1), "V1." + f)))).flatten()
)

*Replication of round 1* on the baseline runs: C0 (clock) and V1 (calibrated
model, test batches) pass for every field; V2 (a priori) fails for `fp64`
(#r(cl(val.at(0), "V2.fp64").value, 2), #r(cl(val.at(2), "V2.fp64").value, 2)) and `fp64bf`
as in round 1, and passes for `w1`–`w4`. On the native build V2 passes for
`fp64` (#r(cl(val.at(1), "V2.fp64").value, 3), #r(cl(val.at(3), "V2.fp64").value, 3)): once
the branches are gone, the op counts explain the kernel. I3's failure for
`fp64` says the same from the other side: the baseline residual
(≈ 27 ns per iteration) is the branch cost, which the native build does not pay.

#figure(image(F + "fig_isa_kernel.svg", width: 100%), caption: [Rho kernel time per
field, runs A1 B1 A2 B2 (medians of 7). Native (red) is lower only for `fp64` and
`fp64bf`.])

= Part B and addendum 1: the `fp64` branch lead

Round 1 found the `fp64` kernel 1.6× slower than `w1`'s although its
operations are cheaper in isolation, and named the data-dependent branches in
`Fp64::sub` and `Fp64::redc` as an unverified lead. This host has no PMU, so
the lead was tested by substitution: the same field, the same values, the
branches replaced.

*F1 was not a test.* `Fp64Bf` wrote the two selections as masks. F1 failed
(worst ratio #r(csum.claims.at(0).value, 3)), but F3 shows why: in the baseline
binary the two kernels are the same instruction sequence. Addendum 1 replaced
the masks with `std::hint::select_unpredictable` (`Fp64Cm`), checked the
machine code *before* freezing, and re-tested.

#table(
  columns: 7,
  align: (left, left, right, right, right, right, right),
  [build], [kernel], [instructions], [cond. jumps], [cmov], [vector-register], [mulx],
  ..f3.map(x => ([#x.label], [#raw(x.function)], [#x.instructions], [#x.jcc], [#x.cmov], [#x.vector], [#x.mulx])).flatten()
)

#table(
  columns: 4,
  align: (left, left, right, left),
  [claim], [rule], [value], [verdict],
  [F1], [`fp64bf` ≤ 0.85 × `fp64`, every B ≥ 4096, A1 and A2], [#r(csum.claims.at(0).value, 4)], [#csum.claims.at(0).verdict (invalid: same code)],
  [F2], [V2 on `fp64bf` ≤ 0.30, A1 and A2], [#r(cl(val.at(0), "V2.fp64bf").value, 3), #r(cl(val.at(2), "V2.fp64bf").value, 3)], [FAIL],
  [F1′], [`fp64cm` ≤ 0.85 × `fp64`, every B ≥ 4096, cm-1 and cm-2], [#r(cmcl.at(0).value, 4)], [#cmcl.at(0).verdict],
  [F2′], [V2 on `fp64cm` ≤ 0.30, cm-1 and cm-2], [#r(cl(cmv.at(0), "V2.fp64cm").value, 3), #r(cl(cmv.at(1), "V2.fp64cm").value, 3)], [#cl(cmv.at(0), "V2.fp64cm").verdict, #cl(cmv.at(1), "V2.fp64cm").verdict],
)

`fp64cm / fp64` lies between
#r(calc.min(..cmp2.filter(x => int(x.batch) >= 4096).map(x => float(x.ratio))), 3) and
#r(calc.max(..cmp2.filter(x => int(x.batch) >= 4096).map(x => float(x.ratio))), 3) at the 12
cells with B ≥ 4096 (all batch sizes: #r(calc.min(..cmp2.map(x => float(x.ratio))), 3)–#r(calc.max(..cmp2.map(x => float(x.ratio))), 3)).
Model C's fitted residual per iteration drops from
#r(cmv.at(0).residual_ns.fp64, 1) ns (`fp64`) to #r(cmv.at(0).residual_ns.fp64cm, 1) ns (`fp64cm`)
in cm-1 (cm-2: #r(cmv.at(1).residual_ns.fp64, 1) and #r(cmv.at(1).residual_ns.fp64cm, 1) ns).

*Reading.* Only the selections differ between `fp64` and `fp64cm`, and the
results are bit for bit equal (unit test and run gate). Removing 12 of 32
conditional jumps halves the kernel's time and makes it what its op counts
predict; on random operands each of those branches is a coin flip. The native
build's 2.1–2.6× `fp64` speedup at B ≥ 256 (part A) is therefore mostly the same effect
reached by LLVM through `cmov` (22 jumps, 14 `cmov`), with `mulx` and vector
moves on top; how the two shares split is not measured here.

*Scope and limits.* Measured on `hwcost::kernel`'s rho step (an r-adding walk
with batched inversion), not on the repository's production rho or IC code;
host L0 (shared cloud VM). Branch misses are inferred from the substitution,
not counted.

#figure(image(F + "fig_fp64_cm.svg", width: 100%), caption: [Runs cm-1 and cm-2: kernel
time, the `fp64cm / fp64` ratio (part A's `fp64bf / fp64` dotted), measured against the
a priori model, and the F3 instruction counts.])

= Part C: a GPU profile

Each of the 46 recorded ECC2K-130 builds (21 RTX PRO 6000, 25 B200;
`ecc2k130/benchmarks/roofline/calibration.json`) has dynamic lane-instruction
counts per update. The profile treats SMs as cores and an update's pipes as
overlapping, so SM-clocks per update are the busiest pipe's count divided by
its measured ceiling (`Phase::pipes_overlap`). Every constant has a source;
issue width 128 and three B200 ceilings are flagged as assumed.

#table(
  columns: 5,
  align: (left, left, right, right, left),
  [id], [statement], [value], [threshold], [verdict],
  ..gcl.map(c => ([#c.id], [#c.statement], [#r(c.value, 4)], [#if c.threshold == none [—] else [#c.threshold]], [#c.verdict])).flatten()
)

The ceiling model is never beaten (no η > 1), as a ceiling should be. The
one-parameter-per-GPU model is good to about 1.5%; the shared-integer-datapath
variant fits better a priori (12%) and across GPUs (11%). Cost check:
#raw(gsum.at(-2)). The paper states \$124k.

#figure(image(F + "fig_gpu.svg", width: 100%), caption: [Predicted against measured
update rate for the 46 builds, separate (filled) and shared-int (open) pipe models.])

= Part D: index calculus in time and dollars

Frozen Phase 1 records of `research/ic_frontier_campaign_20261015` (32
verified runs per arm and curve, host Xeon 2.30 GHz, L0): *Koblitz*
`ic-compact` vs `rho-frob`, *prime* `ic-pair` vs `rho-neg`. The band prices
each arm's work $P(w)$ at $w ∈ {0, 0.1, 1}$ GAE per unpriced native operation;
$w^*$ is the weight at which the band's ratio equals the measured wall ratio.

#for d in ic [
  #parbreak() *#d.name*

  #table(
    columns: 9,
    align: (left, right, right, right, right, right, right, right, left),
    [curve], [log2 r], [IC ms], [rho ms], [wall ratio], [band w=0], [band w=1], [w\*], [§5a rank],
    ..d.pairs.map(x => ([#x.session], [#r(x.log2_r, 1)], [#r(x.ic_wall_ms, 2)], [#r(x.rho_wall_ms, 2)], [#r(x.wall_ratio, 3)],
      [#r(x.band_ratio_w0, 2)], [#r(x.band_ratio_w1, 2)], [#r(x.w_star, 3)], [#x.rank])).flatten()
  )
  #d.sum.map(c => [#c.id = #r(c.value, 3)]).join("; ").
  #parbreak()
]

*Reading.* On small curves rho's wall time sits on a per-solve floor of
about 3–4 ms (setup, not iterations), so IC's wall time is lower on the
smallest ones (Koblitz log2 r ≤ 22, prime 20) while every band weight ranks
IC above rho; §5a leaves those curves unranked. From
log2 r ≈ 28 (Koblitz) and 24 (prime) every curve ranks IC > rho. On Koblitz
the three largest curves fall inside the band ($w^*$ 0.005–0.17): the band's
weights describe IC's wall cost there. On prime the measured ratio leaves the
band from above at log2 r 32, 36 and 40 (by 1.1×, 1.6×, 2.8×): even $w = 1$
charges IC's unpriced work too little, and $w^*$ is meaningless
(the prime band's ratio barely moves with $w$).

*Projections (extrapolated, not measured).* The pre-registered projection fits
$log_2 y = a + alpha log_2 r$ to wall time over the tier. Rho's floor flattens
its wall exponent (α ≈ 0.26 Koblitz, 0.34 prime, against 0.46–0.47 for its
GAE count), so that projection understates rho at large r by orders of
magnitude; the post hoc projection (named in addendum 1) uses the fitted
$P(w)$ exponents at rho's measured ns per GAE instead.

#table(
  columns: 7,
  align: (left, right, right, right, right, right, left),
  [domain], [log2 r], [pre-reg. IC core-yr], [pre-reg. rho core-yr], [post hoc IC core-yr (w=0…1)], [post hoc rho core-yr], [rank],
  ..ic.map(d => d.proj.map(p => {
    let g = d.gae.filter(x => x.log2_r == p.log2_r)
    let ics = g.map(x => float(x.ic_core_years))
    let rhos = g.map(x => float(x.rho_core_years))
    ([#d.name], [#p.log2_r], [#sci(p.ic_core_years)], [#sci(p.rho_core_years)],
     [#sci(calc.min(..ics)) – #sci(calc.max(..ics))],
     [#sci(calc.min(..rhos)) – #sci(calc.max(..rhos))],
     [#p.rank])
  })).flatten().flatten()
)

Dollars per solve at one vCPU of an illustrative c7i.xlarge (\$0.1785/h ÷ 4)
are in `projections.csv` and `projections_gae.csv`. At $2^(128)$ the post hoc
projection prices a Koblitz rho solve at about 1.4–2.4 × $10^(3)$ core-years and
IC at 4.5 × $10^(5)$ to 5.3 × $10^(7)$; on prime, rho 1.3–1.6 × $10^(4)$ and IC
1.4–6.1 × $10^(11)$. The ranks are the same under both projections and every
weight: IC > rho at every target.

#figure(image(F + "fig_ic_koblitz.svg", width: 100%), caption: [Koblitz: measured and
band-predicted wall time (mean of 32 per arm and curve), implied weight, and both
projections.])
#figure(image(F + "fig_ic_prime.svg", width: 100%), caption: [Prime: as above.])

= Claims table

#table(
  columns: 4,
  align: (left, left, right, left),
  [id], [rule], [value], [verdict],
  [I1], [per-arm classification (descriptive)], [—], [counts above],
  [I2], [V1 ≤ 0.30, native build, every field, B1 and B2], [max #r(calc.max(..fields.map(f => cl(val.at(1), "V1." + f).value), ..fields.map(f => cl(val.at(3), "V1." + f).value)), 3)], [PASS],
  [I3], [baseline residual on native, ≤ 0.30, every field], [`fp64` #r(cl(tr.at(0), "V1.fp64").value, 2), `fp64bf` #r(cl(tr.at(0), "V1.fp64bf").value, 2)], [FAIL (`w1`–`w4` pass)],
  [F0], [`fp64bf`, `fp64cm` walks equal `Fp64`'s], [—], [PASS (all runs)],
  [F1], [`fp64bf` ≤ 0.85 × `fp64`], [#r(csum.claims.at(0).value, 3)], [FAIL; invalid test (F3)],
  [F2], [V2 on `fp64bf` ≤ 0.30], [#r(cl(val.at(0), "V2.fp64bf").value, 3)], [FAIL],
  [F3], [disassembly (descriptive)], [—], [`fp64` = `fp64bf` code],
  [F1′], [`fp64cm` ≤ 0.85 × `fp64`], [#r(cmcl.at(0).value, 3)], [#cmcl.at(0).verdict],
  [F2′], [V2 on `fp64cm` ≤ 0.30], [#r(calc.max(cl(cmv.at(0), "V2.fp64cm").value, cl(cmv.at(1), "V2.fp64cm").value), 3)], [PASS],
  [G1], [a priori, 46 builds, ≤ 0.30], [#r(gc("G1").value, 3)], [#gc("G1").verdict],
  [G2], [η per GPU, ≤ 0.30 (and ≤ 0.10)], [#r(gc("G2").value, 4)], [#gc("G2").verdict (#gc("G2.10").verdict at 10%)],
  [G3], [η carried across GPUs, ≤ 0.30], [#r(gc("G3").value, 3)], [#gc("G3").verdict],
  [D1], [in band; $w^*$ (descriptive)], [Koblitz 3/8, prime 0/6], [—],
  [D2], [projections and §5a rank (descriptive)], [—], [IC > rho at every target],
)

= Limits

- Host L0: a shared 4-vCPU cloud VM without pinning, frequency control or a PMU.
  Same-build drift is in every margin; branch misses are inferred, not counted.
- Parts A and B measure `hwcost`'s own rho step, not the production rho or IC.
- Part C validates against recorded GPU rates; no GPU ran in this study.
- Part D's band and wall times come from another host (2.30 GHz) and
  campaign; both projections extrapolate exponents fitted on 4–5 curves
  over 13–18 bits of log2 r to 128 bits.
- The per-curve §5a `rank` column and `projections_gae.csv` were added after
  the analysis; they change no pre-registered value.

= Next steps

+ Port the `select_unpredictable` selections to the production
  `fp_backend::Fp64` in its own PR, with a frozen baseline and the rho
  benchmarks as the gate: a candidate 2× on the baseline-ISA `Fp64` kernel.
+ Re-measure parts A and B on a host with a PMU to count branch misses and
  split the native speedup between `cmov` and `mulx`/vector code.
+ Price IC's memory-bound phases (sieve, linear algebra) with the roofline,
  so the band's unpriced native work gets a measured weight per phase.
