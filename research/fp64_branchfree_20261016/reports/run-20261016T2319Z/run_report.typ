// Run report (descriptive): Fp64 branch-free study, runs started 20261016T2319Z.
// Every number in the tables is read from the run records at build time;
// compile with the study directory as root.
#set page(paper: "a4", margin: 1.6cm, numbering: "1")
#set text(size: 9.3pt)
#set heading(numbering: "1.")
#show table: set text(size: 7.4pt)

#let R = "../../runs/"
#let r(x, d) = str(calc.round(float(x), digits: d))
#let runs = ("base-1", "cand-1", "base-2", "cand-2")
#let meas = runs.map(n => json(R + n + "/measurements.json"))
#let val = runs.map(n => json(R + n + "/validation.json"))
#let txt(n, f) = read(R + n + "/" + f).trim()
#let cmp = csv(R + "compare/compare.csv", row-type: dictionary)
#let csum = json(R + "compare/summary.json")
#let pair = json(R + "cm-pair/claims.json")
#let pairs = csv(R + "cm-pair/pairs.csv", row-type: dictionary)
#let dis = csv(R + "disasm/counts.csv", row-type: dictionary)
#let rw = csv(R + "disasm/rho_worker_counts.csv", row-type: dictionary)
#let samples = runs.map(n => csv(R + n + "/samples.csv", row-type: dictionary).filter(s => s.group == "pfrho"))
#let cl(v, id) = v.claims.find(c => c.id == id)
#let verdict(c) = if c == none [—] else [#r(c.value, 4) #c.verdict]
#let st = meas.at(0).host.stable
#let pc = read(R + "protocol_commit.txt").trim()

= Fp64 branch-free study: run records (runs started 20261016T2319Z)

*Kind:* benchmark runs, descriptive report (observations only; the
interpretation is in the study report `../../report.pdf`). \
*Study:* `research/fp64_branchfree_20261016/`, protocol `PROTOCOL.md` frozen at
the protocol commit #raw(pc.slice(0, 9)) (`runs/protocol_commit.txt`) before
the first run. \
*Commands:* the protocol's "Commands" section, verbatim. Each run is
`hwcost measure --fields fp64,fp64cm,pf64,w1 --prime-fast-rho-bits 44`
(#meas.at(0).params.reps repetitions of about #meas.at(0).params.target_ms ms;
batches #meas.at(0).params.batches.map(str).join(", ")). Every run exited 0
and completed every arm; no arm failed or timed out.

== Runs and binaries

#table(
  columns: 4,
  align: (left, left, left, left),
  [run], [started / finished UTC], [binary SHA-256 (first 16)], [source],
  ..runs.map(n => {
    let b = txt(n, "binary.txt").split("\n")
    ([#n], [#txt(n, "started_utc.txt") / #txt(n, "finished_utc.txt")],
     [#raw(b.at(1).replace("sha256: ", "").slice(0, 16))], [#b.at(2).replace("source: ", "")])
  }).flatten()
)

Both binaries: `cargo build --release --bin hwcost` (default x86-64 target;
compiled target features: #if meas.at(0).compiled_features.len() == 0 [none
probed] else [#meas.at(0).compiled_features.join(", ")]).

== Host and isolation

- CPU: #st.cpu_model (#st.cpu_vendor, signature #st.cpu_signature), #st.logical_cpus
  logical CPUs; virtualization: #st.virtualization.pairs().map(((k, x)) => k + " = " + x).join(", "); kernel #st.kernel_release; #st.rustc_vv.split("\n").at(0).
- Isolation *L0*: shared cloud container, no pinning, no frequency control;
  arms interleaved, order rotated per repetition. No other build or test ran
  during the four runs.
- *No hardware performance counters* (guest kernel: no PMU driver for this
  CPU): branch misses are missing data, not zero.

== Correctness gates

Every run's `stderr.txt` opens with one line per field: the batched walk
equals the unbatched reference and stays on the curve, and the timed
arithmetic's walk equals `fp_backend::Fp64`'s at the same prime point for
point (for `pf64`, `prime_fast::Fp64` and `fp_backend::Fp64` share the
Montgomery form). All passed in all four runs. Every `pfrho` sample is a
known-answer solve whose logarithm was asserted. Before the runs, the
targeted unit tests (92, including both new bit-identity tests) passed on the
protocol commit.

== Claims (computed by `hwcost compare` and `hwcost pair`)

#table(
  columns: 5,
  align: (left, left, right, right, left),
  [id], [statement], [value], [threshold], [verdict],
  ..(csum.claims + pair).map(c => ([#c.id], [#c.statement], [#r(c.value, 4)], [#c.threshold], [#c.verdict])).flatten()
)

Counts of arm classes (candidate vs baseline, effect called only beyond
max(2·drift, 2%) on the same side in both pairs):

#table(
  columns: 3,
  align: (left, left, right),
  [group], [class], [arms],
  ..csum.classes.map(c => ([#c.group], [#c.class], [#c.arms])).flatten()
)

== Rho kernel, B ≥ 4096 (ns per iteration, medians of 7)

#table(
  columns: 10,
  align: (left,) + (right,) * 8 + (left,),
  [arm], [A1], [B1], [A2], [B2], [r1], [r2], [drift], [margin], [class],
  ..cmp.filter(x => x.arm.starts-with("kernel.") and int(x.arm.split(".").at(2)) >= 4096).map(x => (
    [#x.arm.replace("kernel.", "")], [#r(x.a1, 1)], [#r(x.b1, 1)], [#r(x.a2, 1)], [#r(x.b2, 1)],
    [#r(x.r1, 3)], [#r(x.r2, 3)], [#r(x.drift, 3)], [#r(x.margin, 3)], [#x.class],
  )).flatten()
)

#figure(image("fig_kernel.svg", width: 100%), caption: [Rho kernel time per field
against batch size in all four runs (medians of 7 repetitions).])

== Field operations and the production rho

#table(
  columns: 10,
  align: (left,) + (right,) * 8 + (left,),
  [arm], [A1], [B1], [A2], [B2], [r1], [r2], [drift], [margin], [class],
  ..cmp.filter(x => x.arm.split(".").at(0) in ("fp64", "pf64", "pfrho")).map(x => (
    [#x.arm], [#r(x.a1, 2)], [#r(x.b1, 2)], [#r(x.a2, 2)], [#r(x.b2, 2)],
    [#r(x.r1, 3)], [#r(x.r2, 3)], [#r(x.drift, 3)], [#r(x.margin, 3)], [#x.class],
  )).flatten()
)

Units: field operations ns per operation; `pfrho.44` ns per walk step of
`prime_fast::rho_parallel` (one thread, one known-answer solve per sample).
Every `pfrho.44` sample (ns per step):

#table(
  columns: 2,
  align: (left, left),
  [run], [samples (repetition order)],
  ..runs.enumerate().map(((i, n)) => ([#n], [#samples.at(i).map(s => r(s.value, 2)).join(", ")])).flatten()
)

#figure(image("fig_ratios.svg", width: 100%), caption: [Candidate / baseline ratio
for every arm, $r_1$ and $r_2$, with the drift margin.])

#figure(image("fig_pfrho.svg", width: 70%), caption: [`pfrho.44`: every sample of
the production rho in each run.])

== `fp64` against `fp64cm` in the candidate runs (P3)

#table(
  columns: 5,
  align: (left, right, right, right, right),
  [run], [B], [fp64cm ns], [fp64 ns], [fp64 / fp64cm],
  ..pairs.filter(x => int(x.batch) >= 4096).map(x => ([#x.run], [#x.batch], [#r(x.fp64cm_ns, 2)], [#r(x.fp64_ns, 2)], [#r(x.ratio, 4)])).flatten()
)

== Disassembly (D2)

#table(
  columns: 6,
  align: (left, left, right, right, right, right),
  [binary], [function], [instructions], [cond. jumps], [cmov], [mulx],
  ..dis.map(x => ([#x.label], [#raw(x.function)], [#x.instructions], [#x.jcc], [#x.cmov], [#x.mulx])).flatten()
)

`step_pf64` is linked in neither binary (`disasm.sh` reports it absent: no
caller), so `prime_fast::rho_worker`, the production rho's per-thread loop,
is counted instead (`disasm_rho_worker.sh`, whole symbol by its `nm -S` size):

#table(
  columns: 6,
  align: (left, left, right, right, right, right),
  [binary], [function], [instructions], [cond. jumps], [cmov], [mulx],
  ..rw.map(x => ([#x.label], [#raw(x.function)], [#x.instructions], [#x.jcc], [#x.cmov], [#x.mulx])).flatten()
)

#figure(image("fig_disasm.svg", width: 80%), caption: [Instruction counts of the
non-inlined kernel step and of `prime_fast::rho_worker` in each binary
(`runs/disasm`).])

== Validation of the cost model per run (descriptive)

#let ids = val.at(0).claims.map(c => c.id)
#table(
  columns: 5,
  align: (left, right, right, right, right),
  [claim], ..runs.map(n => [#n]),
  ..ids.map(id => ([#id], ..range(4).map(i => verdict(cl(val.at(i), id))))).flatten()
)

== Raw records

`runs/{base,cand}-{1,2}/` (`measurements.json`, `samples.csv`, `medians.csv`,
`validation.{json,csv}`, `profile.json`, `fit.txt`, `stderr.txt`,
`binary.txt`, timestamps), `runs/compare/`, `runs/cm-pair/`, `runs/disasm/`,
`runs/protocol_commit.txt`. Graph sources: `../../plot.py` (layout only).
