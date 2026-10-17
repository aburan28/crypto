// Results report: hwcost, a calibrated analytical time-and-cost model,
// validated on the Pollard rho iteration. Tables read the run records.
#set page(paper: "a4", margin: 1.8cm, numbering: "1")
#set text(size: 9.5pt)
#set heading(numbering: "1.")
#show table: set text(size: 8pt)

#let run = "runs/final/"
#let rep = "reports/run-20261016T0459Z/"
#let m = json(run + "measurements.json")
#let v = json(run + "validation.json")
#let hw = json(run + "profile.json")
#let rho = csv(run + "rho.csv", row-type: dictionary)
#let rhon = csv(run + "rho_negation.csv", row-type: dictionary)
#let commit = read(run + "commit.txt").trim()
#let started = read(run + "started_utc.txt").trim()
#let r(x, d) = str(calc.round(x, digits: d))
#let pct(x) = r(100 * x, 1) + "%"
#let fields = ("fp64", "w1", "w2", "w3", "w4")
#let usd(x) = if x < 0.01 [< 0.01] else if x < 1e6 [#r(x, 2)] else [2^#r(calc.log(x, base: 2), 1)]

#align(center)[
  #text(size: 14pt, weight: "bold")[hwcost: a calibrated analytical model of compute time and cost] \
  #text(size: 10pt)[validated on the Pollard rho iteration — results report, 2026-10-16]
]

*Study:* `research/hwcost_model_20261016/` (protocol `PROTOCOL.md`, frozen at
#raw(commit.slice(0, 9)) before the run). *Run:* `runs/final`, started #started
UTC; descriptive run report `reports/run-20261016T0459Z/run_report.pdf`. *Code:* `src/cryptanalysis/hwcost/`,
CLI `src/bin/hwcost.rs`.

= What the tool does

A hardware profile and a workload go in; a time and cost breakdown comes out.

- *Compute*: iterations × Σ (operations per iteration × cycles per operation) / clock.
  Cycles per operation come from the profile: measured on the hot loop by
  `hwcost measure`, or entered from instruction tables, llvm-mca or uiCA with
  the source named. An operation on a serial chain is priced at its latency
  (`op@latency`), the rest at reciprocal throughput.
- *Memory*: bytes moved per iteration / bandwidth of the first level that holds
  the per-core working set; a phase takes the larger of compute and memory
  time (the roofline, applied to a CPU).
- *Parallelism*: Amdahl with the *measured* scaling efficiency, interpolated
  between measured core counts and flagged when extrapolated.
- *Transfers*: messages × latency + bytes / bandwidth per data movement;
  transfers that overlap compute take a max with it, others add.
- *Cost*: wall hours × machines × \$/h, plus egress and storage.

Every constant carries its source; every prediction lists the unmeasured
constants it rests on. An operation the workload names but the profile does
not price is an error, never a zero.

= Calibration and validation on the rho iteration

The rho kernel (`hwcost::kernel`) advances `B` affine walks of an r-adding walk
together with one Montgomery batch inversion. Its operation counts per iteration
were counted, not asserted: 3 multiplications with independent work around
them, 2 on the batch inversion's prefix-product chain, 1 squaring, 6
additions/subtractions and 1/B inversion (counted step on the timed walks,
every field and batch size). Four predictions were compared with the measured
median ns per iteration (protocol, "Models"): *A* a priori with no fitted
parameter; *C* = A plus one per-iteration residual per field fitted on
B ≤ #v.fit_max_batch and tested on larger B; *S* for `w4` (P-256 field) with
every operation cost fitted against limb count on 1–3 limbs and extrapolated.

#table(
  columns: (auto, 1fr, auto, auto, auto),
  align: (left, left, right, right, center),
  [claim], [statement], [value], [threshold], [verdict],
  ..v.claims.map(c => ([#c.id], [#c.statement], [#r(c.value, 4)], [≤ #c.threshold], [*#c.verdict*])).flatten()
)

#figure(image(rep + "fig_validation.svg", width: 100%), caption: [Measured (median, IQR, n = #m.params.reps) and predicted ns per rho iteration against batch size; grey is the fit set. Run `runs/final`, #started.])

#figure(image(rep + "fig_errors.svg", width: 82%), caption: [Prediction error per field and batch size; band ±30%.])

== What the results say

- *Calibrated (C): every field passes*, with median absolute error between
  #pct(calc.min(..fields.map(f => v.claims.find(c => c.id == "V1." + f).value))) and
  #pct(calc.max(..fields.map(f => v.claims.find(c => c.id == "V1." + f).value))) on the
  test batches, 4× to 4096× larger than the largest batch the residual was
  fitted on. The worst single test row is
  #pct(calc.max(..v.rows.filter(x => x.set == "test").map(x => calc.abs(x.err_calibrated)))).
  One fitted constant per field is enough for this kernel on this host, in line
  with the 10–30% expectation for a calibrated compute-bound model, at the low
  end.
- *A priori (A): 4 of 5 fields pass, `fp64` fails* (median |error|
  #pct(v.claims.find(c => c.id == "V2.fp64").value)). On the `FpWide` fields A's
  error lies between
  #pct(calc.min(..v.rows.filter(x => x.field != "fp64").map(x => x.err_a_priori))) and
  #pct(calc.max(..v.rows.filter(x => x.field != "fp64").map(x => x.err_a_priori)))
  (it under-predicts `w2`–`w4` at every batch size); on `fp64` it lies between
  #pct(calc.min(..v.rows.filter(x => x.field == "fp64").map(x => x.err_a_priori))) and
  #pct(calc.max(..v.rows.filter(x => x.field == "fp64").map(x => x.err_a_priori))):
  the microbenchmarked operation costs explain only about half of the measured
  time there. The fitted residuals
  (measured − A, ns per iteration):
  #fields.map(f => raw(f) + " " + r(v.residual_ns.at(f), 1)).join(", ").
  The `fp64` backend's operations are cheaper than `w1`'s in isolation (mul
  #r(m.samples.filter(s => s.arm == "fp64.mul.throughput").map(s => s.value).sorted().at(3), 2)
  vs #r(m.samples.filter(s => s.arm == "w1.mul.throughput").map(s => s.value).sorted().at(3), 2) ns throughput),
  yet its kernel is about 1.6× slower (≈ 45 vs ≈ 28 ns per iteration). Its
  operations are inlined (no out-of-line `Fp64::mul`/`sub` symbol), so call
  overhead is excluded. *Unverified lead:* `Fp64::sub` and `Fp64::redc` take
  data-dependent `if` branches that a microbenchmark's single chain may
  predict but random walk data may not; no `perf` counters were available to
  check. This is exactly the gap the calibration step exists to close.
- *Size extrapolation (S)*: fitting every operation cost against limb count on
  1–3 limbs and the residual linearly, then extrapolating to 4 limbs, predicts
  the P-256-field kernel with median |error|
  #pct(v.claims.find(c => c.id == "V3.w4").value) (over-predicting by 1–12%).
  The fitted costs at n = 4 come out above the measured ones
  (fig. 4 of the run report: e.g. mul throughput extrapolated
  #r(v.size_fits.find(s => s.op == "mul.throughput").extrapolated_ns, 1) ns vs measured
  #r(v.size_fits.find(s => s.op == "mul.throughput").measured_ns, 1) ns).
- *Memory roofline*: the model assigns every arm to compute (the kernel draws
  at most ≈ 3.5 GB/s against ≥ 8.8 GB/s single-thread DRAM read bandwidth),
  and the measurement agrees to first order: once the inversion is amortized,
  growing the working set from L2 to beyond the last-level cache moves the
  per-iteration time by 2–10% depending on the field (B = 2^16 to 2^22). That
  rise is a cost the `max()` roofline does not represent, and C's error drifts
  negative there (down to
  #pct(calc.min(..v.rows.map(x => x.err_calibrated)))).
- *Clock*: the reg-reg `add` and `imul` chains agree within
  #pct(v.claims.find(c => c.id == "C0").value) at #r(v.clock_ghz, 3) GHz (the
  nominal clock is 2.1 GHz). A chain of `add reg, imm` reads about 12 GHz on
  this CPU: dependent immediate adds are folded at rename, so such a chain is
  not a clock measurement.
- *Scaling*: efficiency at #m.params.threads.last() threads
  #fields.map(f => raw(f) + " " + r(v.scaling.at(f).at(str(m.params.threads.last())), 3)).join(", ")
  (L0 host; the 4 vCPUs may share physical cores). The profile carries the
  median over fields.

#figure(image(rep + "fig_roofline.svg", width: 70%), caption: [Bandwidth available (sequential read) vs bandwidth the kernel drew.])

= Application: parallel rho on prime-order curves (descriptive, no claim)

`hwcost rho` turns the kernel's counted operations into a van Oorschot–Wiener
workload. It uses √(πn/2) expected iterations (√(πn/4) with the negation map,
fruitless cycles unpriced) plus the distinguished-point tail. Distinguished-point
bits are chosen so the tail stays ≤ 1%. Start points are about 1.5 log₂ n group
operations per walk. Distinguished points are uploaded over each machine's link
and stored for the run. The workload is predicted on the fitted profile at
B = 1024 walks per core, 4 cores per machine, with the field the smallest
measured one holding the group. The price is illustrative (AWS us-east-1
on-demand c7i.xlarge, \$0.1785/h); the network and disk constants are assumed
(the output lists them).

#table(
  columns: 8,
  align: (right, left, right, right, right, right, right, right),
  [log₂ n], [field], [iterations (log₂)], [ns/iter], [1 machine], [1,000 machines], [cost, 1 machine (\$)], [same, negation map (\$)],
  ..rho.filter(e => e.machines == "1").map(e => {
    let k = rho.find(x => x.group_log2 == e.group_log2 and x.machines == "1000")
    let n = rhon.find(x => x.group_log2 == e.group_log2 and x.machines == "1")
    (
      [#e.group_log2], [#e.field], [#r(float(e.iterations_log2), 2)], [#r(float(e.ns_per_iter_compute), 1)],
      [#e.wall_human], [#k.wall_human],
      [#usd(float(e.total_usd))],
      [#usd(float(n.total_usd))],
    )
  }).flatten()
)

#figure(image(rep + "fig_rho.svg", width: 100%), caption: [Predicted wall time and cost of parallel rho against group size on this profile (4 cores per machine). Extrapolated from run `runs/final`; price and network illustrative/assumed.])

Reading the table: the time per iteration grows with the limb count of the
field, so cost tracks √n × (cost of one field's arithmetic). At n = 2^48 on
1,000 machines the start-point setup dominates (4 million walks at about 72
group operations each, against 2^24.6 walk iterations), so 1,000 machines are
only #{
  let a = rho.find(x => x.group_log2 == "48" and x.machines == "1")
  let b = rho.find(x => x.group_log2 == "48" and x.machines == "1000")
  r(float(a.wall_s) / float(b.wall_s), 0)
}× faster than one. At 2^112
(the size of the solved secp112r1 instance) the model gives about 3 weeks on
1,000 four-core machines; at 2^256 (P-256) it gives about 2^69 years on the same
fleet. Both are extrapolations of a 3-minute calibration of one 4-vCPU host and
one kernel; they inherit every limit below.

= Scope and limits

- *One host, isolation L0* (shared 4-vCPU KVM guest, no pinning, turbo
  uncontrolled), *one kernel*, *one build* (x86-64 baseline: no BMI2/ADX/AVX2
  code generation even though the CPU has them). An ISA-specific build is a
  different profile and needs its own `hwcost measure` run.
- The rho kernel does not restart walks at distinguished points and solves no
  logarithm; the application's DP tail, start points, network and storage
  terms are modelled, not measured.
- GPU profiles are supported by the schema (SMs as cores, instruction
  classes as operations) but none was measured here.
- Not an index-calculus result, so no gap-to-rho row; no two ECDLP methods
  are ranked, so the AGENTS.md §5a unpriced-work band does not apply.
- Correction during the study: the application first capped distinguished-point
  bits at 64, which for groups above 2^192 left up to 2^64 points stored and
  made storage dominate the cost. The cap was removed before the reported
  `rho*.csv` were produced (commit after #raw(commit.slice(0, 9)); kernel and
  measurement code unchanged).

= Next steps

+ Run `hwcost measure` on a `-C target-cpu=native` build (MULX/ADX, AVX-512 IFMA
  `Fp52`) and compare profiles: the ISA effect on these field operations.
+ Check the `fp64` lead with hardware counters (branch misses) on a host with
  `perf`; if confirmed, a branch-free `Fp64::sub`/`redc` is an engineering
  change to measure on the rho kernel.
+ Add a GPU profile from the measured pipe ceilings in
  `benchmarks/hardware-limits/` and the ECC2K-130 per-pipe instruction counts
  (`ecc2k130/ROOFLINE.md`), and validate it against the recorded 15.1 B it/s.
+ Price an index-calculus workload phase by phase (relation search, linear
  algebra, descent) from ecbench's counted units, so IC and rho are compared
  in time and dollars as well as in GAE.
