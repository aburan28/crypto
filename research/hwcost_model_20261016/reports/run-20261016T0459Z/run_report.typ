// Run report (descriptive): hwcost calibration run 20261016T0459Z.
// Every number in the tables is read from the run directory at build time.
#set page(paper: "a4", margin: 1.8cm, numbering: "1")
#set text(size: 9.5pt)
#set heading(numbering: "1.")
#show table: set text(size: 8pt)

#let run = "../../runs/final/"
#let m = json(run + "measurements.json")
#let v = json(run + "validation.json")
#let med = csv(run + "medians.csv", row-type: dictionary)
#let rows = csv(run + "validation.csv", row-type: dictionary)
#let commit = read(run + "commit.txt").trim()
#let started = read(run + "started_utc.txt").trim()
#let finished = read(run + "finished_utc.txt").trim()
#let r(x, d) = str(calc.round(x, digits: d))
#let pct(x) = r(100 * x, 1) + "%"
#let mv(arm) = {
  let hit = med.find(e => e.arm == arm)
  if hit == none { none } else { float(hit.median) }
}
#let fields = ("fp64", "w1", "w2", "w3", "w4")

= hwcost calibration run `runs/final` (started #started UTC)

*Kind:* benchmark / calibration run, descriptive report (observations only; the
interpretation is in the study report `../../report.pdf`). \
*Study:* `research/hwcost_model_20261016/`, protocol `PROTOCOL.md` frozen at
commit #raw(commit) before this run. \
*Command:* `./target/release/hwcost measure --out research/hwcost_model_20261016/runs/final`,
then `hwcost fit --run … --fit-max-batch` #v.fit_max_batch. Binary built from
commit #raw(commit.slice(0, 9)) (`cargo build --release --bin hwcost`). Started
#started, finished #finished (UTC); exit status 0. \
*Runs requested / completed:* every arm below, #m.params.reps repetitions each,
all completed; no arm failed or timed out.

== Host and isolation

#let st = m.host.stable
- CPU: #st.cpu_model (#st.cpu_vendor, signature #st.cpu_signature), #st.logical_cpus
  logical CPUs, #st.numa_nodes NUMA node(s); virtualization: #st.virtualization.pairs().map(((k, x)) => k + " = " + x).join(", ").
- Caches (sysfs, as the guest sees them):
  #m.caches.filter(c => c.kind != "Instruction").map(c => [L#c.level #c.kind #r(c.size_bytes / 1024, 0) KiB (CPUs #c.shared_cpu_list)]).join(", ").
- Binary compiled with target features: #if m.compiled_features.len() == 0 [none of those probed (SSE4.2, AVX, AVX2, FMA, BMI1/2, ADX, AVX-512F/IFMA, PCLMULQDQ): the default x86-64 baseline] else [#m.compiled_features.join(", ")]. The CPU itself supports BMI2, ADX, AVX2 and AVX-512 (host capsule features); this build does not use them.
- Isolation level *L0*: shared cloud container, no core pinning, no frequency control.
  Arms were interleaved and their order rotated per repetition.
- Host capsule class id: #raw(m.host.env_class_id).

== Correctness gates (before any timing)

For each field the batched walk (one batch inversion) was run 4 steps on 33
walks next to an unbatched reference (one inversion per walk); both agreed and
every point satisfied $y^2 = x^3 + 7$ (stderr: `runs/final/stderr.txt`). Every
timed kernel arm asserted zero degenerate steps.

== Measured values

All values are medians over #m.params.reps repetitions (IQR where shown),
each repetition about #m.params.target_ms ms.

*Clock (M0).* Dependent reg-reg `add` chain: #r(mv("clock.add_chain"), 3) GHz;
dependent `imul` chain (3-cycle latency): #r(mv("clock.imul_chain"), 3) GHz.
Cycles below are ns × #r(v.clock_ghz, 3) GHz.

*Field operations (M1)*, ns per operation (cycles in parentheses):

#table(
  columns: 9,
  align: (left,) + (right,) * 8,
  [field], [mul thr], [mul lat], [sqr thr], [sqr lat], [add thr], [add lat], [inv thr], [inv lat],
  ..fields.map(f => (
    [#f],
    ..("mul", "sqr", "add", "inv").map(op => ("throughput", "latency").map(mode => {
      let x = mv(f + "." + op + "." + mode)
      [#r(x, if op == "inv" { 0 } else { 2 }) (#r(x * v.clock_ghz, 0))]
    })).flatten(),
  )).flatten()
)

*Memory read bandwidth (M2)*, GB/s, sequential read of the given size per thread:

#table(
  columns: 1 + m.params.mem_sizes.len(),
  align: (left,) + (right,) * m.params.mem_sizes.len(),
  [threads],
  ..m.params.mem_sizes.map(s => if s >= 1048576 [#r(s / 1048576, 0) MiB] else [#r(s / 1024, 0) KiB]),
  ..(1, m.params.threads.last()).map(t => (
    [#t],
    ..m.params.mem_sizes.map(s => {
      let x = mv("mem." + str(t) + "t." + str(s))
      if x == none [–] else [#r(x, 1)]
    }),
  )).flatten()
)

*Rho kernel (M3)*, ns per iteration on one thread, median and IQR, with the
exact operation counts per iteration from the counted step (M5):

#table(
  columns: 9,
  align: (left, right, right, right, right, right, right, right, left),
  [field], [B], [median ns], [IQR], [mul], [mul (chain)], [sqr / add], [inv], [model memory level],
  ..rows.map(e => {
    let c = m.counts.find(k => k.field == e.field and str(k.batch) == e.batch)
    (
      [#e.field], [#e.batch], [#r(float(e.measured_ns), 1)],
      [#r(float(e.measured_q1), 1)–#r(float(e.measured_q3), 1)],
      [#r(c.mul, 3)], [#r(c.mul_dep, 3)], [#r(c.sqr, 0) / #r(c.add, 0)], [1/#r(1 / c.inv, 0)], [#e.memory_level],
    )
  }).flatten()
)

*Thread scaling (M4)*, B = #m.params.scaling_batch walks per thread. The
efficiency is ns(1 thread) / (T × aggregate ns per iteration at T threads):

#table(
  columns: 1 + m.params.threads.len(),
  align: (left,) + (right,) * m.params.threads.len(),
  [field], ..m.params.threads.map(t => [T = #t]),
  ..fields.map(f => ([#f], ..m.params.threads.map(t => [#r(v.scaling.at(f).at(str(t)), 3)]))).flatten()
)

== Graphs

#figure(image("fig_validation.svg", width: 100%), caption: [Measured ns per iteration (median, IQR; n = #m.params.reps per point) against batch size, with the protocol's model predictions A, C and S drawn over it. Grey: fit set. Run `runs/final`, #started.])

#figure(image("fig_errors.svg", width: 85%), caption: [Prediction error `predicted / measured − 1` of models A, C and S per field and batch size; the band is the protocol's ±30% threshold.])

#figure(image("fig_roofline.svg", width: 75%), caption: [Measured sequential-read bandwidth (1 and #m.params.threads.last() threads, aggregate) against size per thread, and the bandwidth each kernel arm drew (12 element-sized moves per iteration divided by its measured time).])

#figure(image("fig_opcost.svg", width: 100%), caption: [Measured cycles per field operation against limb count for `FpWide<n>`, with the size fits on n = 1..3 and their extrapolation to n = 4 (stars); `fp64` (dedicated backend) for reference.])

== Raw records

`runs/final/`: `measurements.json` (parameters, host capsule, caches, counts,
every sample), `samples.csv` (one row per arm and repetition), `medians.csv`,
`stderr.txt`, `commit.txt`, `started_utc.txt`, `finished_utc.txt`; fit outputs
`profile.json`, `validation.json`, `validation.csv`, `fit.txt`. Graph source:
`../../plot.py` (layout of the recorded values only).
