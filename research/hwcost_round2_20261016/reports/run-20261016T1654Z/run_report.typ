// Run report (descriptive): hwcost round 2, runs started 20261016T1654Z.
// Every number in the tables is read from the run records at build time.
#set page(paper: "a4", margin: 1.6cm, numbering: "1")
#set text(size: 9.3pt)
#set heading(numbering: "1.")
#show table: set text(size: 7.4pt)

#let R = "../../runs/"
#let r(x, d) = str(calc.round(float(x), digits: d))
#let isa = ("isa-base-1", "isa-native-1", "isa-base-2", "isa-native-2")
#let cm = ("cm-1", "cm-2")
#let runs = isa + cm
#let meas = runs.map(n => json(R + n + "/measurements.json"))
#let val = runs.map(n => json(R + n + "/validation.json"))
#let txt(n, f) = read(R + n + "/" + f).trim()
#let tr = ("isa-native-1", "isa-native-2").map(n => json(R + n + "/transfer/validation.json"))
#let cmp = csv(R + "compare/compare.csv", row-type: dictionary)
#let csum = json(R + "compare/summary.json")
#let fpp = csv(R + "compare/fp64_pairs.csv", row-type: dictionary)
#let cmp2 = csv(R + "cm-pair/pairs.csv", row-type: dictionary)
#let cmcl = json(R + "cm-pair/claims.json")
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
#let verdict(c) = if c == none [—] else [#r(c.value, 4) #c.verdict]
#let st = meas.at(0).host.stable

= hwcost round 2: run records (runs started 20261016T1654Z)

*Kind:* benchmark and analysis runs, descriptive report (observations only;
interpretation is in the study report `../../report.pdf`). \
*Study:* `research/hwcost_round2_20261016/`, protocol `PROTOCOL.md` frozen at
commit `597502b4d` before parts A–D ran; addendum 1 frozen at `4e5994fc9`
before runs cm-1 and cm-2. \
*Commands:* the protocol's "Commands" sections, verbatim. Every run below
exited 0 and completed every requested arm (#meas.at(0).params.reps
repetitions of about #meas.at(0).params.target_ms ms each, batches
#meas.at(0).params.batches.map(str).join(", ")); no arm failed or timed out.

== Runs

#table(
  columns: 6,
  align: (left, left, left, left, left, left),
  [run], [binary], [fields], [started UTC], [finished UTC], [checkout HEAD (`commit.txt`)],
  ..runs.enumerate().map(((i, n)) => (
    [#n],
    [#if n.contains("native") [`target-native` (`-C target-cpu=native`)] else [`target` (default x86-64)]],
    [#meas.at(i).params.fields.join(", ")],
    [#txt(n, "started_utc.txt")],
    [#txt(n, "finished_utc.txt")],
    [#raw(txt(n, "commit.txt").slice(0, 9))],
  )).flatten()
)

`commit.txt` is the checkout's HEAD when `measure` started. Both part-A
binaries were built at 16:54 UTC from the protocol commit's tree; `b427150af`
(committed during the runs) changes only claim printing in `compare`,
`gpu-check` and `ic-price` and `fit.rs`'s claim helpers: the timed code
(`kernel.rs`, `measure.rs`, `mod.rs`, `rho.rs`) is identical at both commits
(`git diff --stat 597502b4d b427150af`). cm-1 and cm-2 ran a baseline binary
built at `4e5994fc9`.

Compiled target features: baseline runs
#if meas.at(0).compiled_features.len() == 0 [none probed (x86-64 baseline)];
native runs #meas.at(1).compiled_features.join(", ").

== Host and isolation

- CPU: #st.cpu_model (#st.cpu_vendor, signature #st.cpu_signature), #st.logical_cpus
  logical CPUs, #st.numa_nodes NUMA node; virtualization: #st.virtualization.pairs().map(((k, x)) => k + " = " + x).join(", "); kernel #st.kernel_release; #st.rustc_vv.split("\n").at(0).
- Caches (sysfs): #meas.at(0).caches.filter(c => c.kind != "Instruction").map(c => [L#c.level #c.kind #r(c.size_bytes / 1024, 0) KiB]).join(", ").
- Isolation *L0*: shared cloud container, no pinning, no frequency control; arms interleaved, order rotated per repetition. Host class #raw(meas.at(0).host.env_class_id).
- *No hardware performance counters:* the guest kernel reports `Performance
  Events: unsupported CPU family 6 model 207 no PMU driver, software events
  only.` Branch misses, retired instructions and cycles could not be counted.
  This is missing data, not a measured zero.

== Correctness gates (before any timing)

Every run's `stderr.txt` opens with one line per field: the batched walk equals
the unbatched reference and stays on the curve, and the timed arithmetic
equals the field's (for `fp64bf` and `fp64cm`: the walk is `Fp64`'s point for
point, claim F0). All gates passed in all six runs; the unit test
`branch_free_fp64_is_bit_for_bit_fp64` compares 10,000 random operations of
each variant with `Fp64`.

== Part A: baseline vs native build

*Validation claims per run* (median |relative error|; threshold 0.30; C0
threshold 0.05). V1 = calibrated model on test batches B ≥ 4096, V2 = a priori
model, V3 = size-extrapolated model.

#let ids = val.at(0).claims.map(c => c.id)
#table(
  columns: 5,
  align: (left, right, right, right, right),
  [claim], ..isa.map(n => [#n]),
  ..ids.map(id => ([#id], ..range(4).map(i => verdict(cl(val.at(i), id))))).flatten()
)

*Calibration transfer (I3)*: native op costs plus isa-base-1's fitted
residuals (`hwcost fit --residual-from`), median |error| on test batches:

#table(
  columns: 3,
  align: (left, right, right),
  [claim], [isa-native-1], [isa-native-2],
  ..tr.at(0).claims.filter(c => c.id.starts-with("V1")).map(c => ([#c.id], verdict(c), verdict(cl(tr.at(1), c.id)))).flatten()
)

*Per-arm comparison (I1)* (`runs/compare/compare.csv`, #(cmp.len()) arms):
$r_k = B_k \/ A_k$, drift $d = |A_2\/A_1 - 1|$, an effect called only when both
ratios clear $max(2d, 0.02)$ on the same side. Counts per arm group:

#table(
  columns: 3,
  align: (left, left, right),
  [group], [class], [arms],
  ..csum.classes.map(c => ([#c.group], [#c.class], [#c.arms])).flatten()
)

Kernel arms at four batch sizes (ns per rho iteration, medians of 7):

#let kb = ("64", "4096", "262144", "4194304")
#table(
  columns: 10,
  align: (left,) + (right,) * 8 + (left,),
  [arm], [A1], [B1], [A2], [B2], [r1], [r2], [drift], [margin], [class],
  ..cmp.filter(x => x.arm.starts-with("kernel.") and kb.contains(x.arm.split(".").at(2))).map(x => (
    [#x.arm.replace("kernel.", "")], [#r(x.a1, 1)], [#r(x.b1, 1)], [#r(x.a2, 1)], [#r(x.b2, 1)],
    [#r(x.r1, 3)], [#r(x.r2, 3)], [#r(x.drift, 3)], [#r(x.margin, 3)], [#x.class],
  )).flatten()
)

#figure(image("fig_isa_kernel.svg", width: 100%), caption: [Rho kernel time per field
against batch size, runs isa-base-1, isa-native-1, isa-base-2, isa-native-2 (medians
of 7 repetitions each).])

#figure(image("fig_isa_ratios.svg", width: 100%), caption: [Native / baseline ratio per
arm, $r_1$ and $r_2$, with the drift margin, all #(cmp.len()) arms of `compare.csv`.])

== Part B: `fp64` vs `fp64bf` (baseline runs)

`fp64bf / fp64` kernel ratio at B ≥ 4096 in isa-base-1 and isa-base-2 (F1 cells):

#table(
  columns: 5,
  align: (left, right, right, right, right),
  [run], [B], [fp64 ns], [fp64bf ns], [fp64bf / fp64],
  ..fpp.filter(x => int(x.batch) >= 4096 and x.run.contains("base")).map(x => ([#x.run], [#x.batch], [#r(x.fp64_ns, 2)], [#r(x.fp64bf_ns, 2)], [#r(x.ratio, 4)])).flatten()
)

#{
  let c = csum.claims.at(0)
  [*#c.id* (#c.statement): #r(c.value, 4), threshold #c.threshold, #c.verdict.]
}

*F3, disassembly* (`runs/f3/`, `f3_disasm.sh`; `objdump -d` of each
`#[inline(never)]` kernel function; binaries in `runs/f3/binaries.txt`):

#table(
  columns: 7,
  align: (left, left, right, right, right, right, right),
  [build], [function], [instructions], [cond. jumps], [cmov], [vector-register], [mulx],
  ..f3.map(x => ([#x.label], [#raw(x.function)], [#x.instructions], [#x.jcc], [#x.cmov], [#x.vector], [#x.mulx])).flatten()
)

#figure(image("fig_fp64_branches.svg", width: 85%), caption: [`fp64` and `fp64bf`
kernels in all four part-A runs, and their ratio; grey band = F1 cells.])

== Addendum 1: `fp64` vs `fp64cm` (runs cm-1, cm-2)

#table(
  columns: 5,
  align: (left, right, right, right, right),
  [run], [B], [fp64 ns], [fp64cm ns], [fp64cm / fp64],
  ..cmp2.map(x => ([#x.run], [#x.batch], [#r(x.fp64_ns, 2)], [#r(x.fp64cm_ns, 2)], [#r(x.ratio, 4)])).flatten()
)

#for c in cmcl [*#c.id* (#c.statement): #r(c.value, 4), threshold #c.threshold, #c.verdict. ]

#table(
  columns: 3,
  align: (left, right, right),
  [claim], [cm-1], [cm-2],
  ..val.at(4).claims.map(c => ([#c.id], verdict(c), verdict(cl(val.at(5), c.id)))).flatten()
)

Fitted residual per iteration (model C): cm-1 fp64 #r(val.at(4).residual_ns.fp64, 1) ns,
fp64cm #r(val.at(4).residual_ns.fp64cm, 1) ns; cm-2 fp64 #r(val.at(5).residual_ns.fp64, 1) ns,
fp64cm #r(val.at(5).residual_ns.fp64cm, 1) ns.

#figure(image("fig_fp64_cm.svg", width: 100%), caption: [Runs cm-1 and cm-2: kernel time,
`fp64cm / fp64` ratio (part A's `fp64bf / fp64` dotted for reference), measured against
the a priori model, and the F3 instruction counts.])

== Part C: GPU profile against 46 ECC2K-130 builds

Input: `ecc2k130/benchmarks/roofline/calibration.json` (recorded rates; no
GPU ran in this study). Output `runs/gpu/` (per-build rows in `gpu_check.csv`).

#table(
  columns: 5,
  align: (left, left, right, right, left),
  [id], [statement], [value], [threshold], [verdict],
  ..gcl.map(c => ([#c.id], [#c.statement], [#r(c.value, 4)], [#if c.threshold == none [—] else [#c.threshold]], [#c.verdict])).flatten()
)

Cost-arithmetic check (`summary.txt`): #gsum.slice(-2).map(l => raw(l)).join(linebreak())

#figure(image("fig_gpu.svg", width: 100%), caption: [Predicted against measured
update rate, 46 builds, both pipe models; the three panels are G1, G2, G3.])

== Part D: IC and rho in wall time (frozen campaign records)

Input: `research/ic_frontier_campaign_20261015` Phase 1 sessions and §5a band
files (32 counted, verified runs per arm and curve; host Xeon 2.30 GHz, 4 vCPU,
L0). No solve ran in this study. Dollars at \$0.1785/h ÷ 4 per vCPU-hour
(illustrative). The `rank` column applies the §5a rule (wall ratio and all three
band ratios on one side of 1) per measured curve; it was added to `ic-price`
after the analysis and changes no other value.

#for d in ic [
  *#d.name*: #d.pairs.at(0).ic_arm vs #d.pairs.at(0).rho_arm.

  #table(
    columns: 11,
    align: (left,) + (right,) * 9 + (left,),
    [curve], [log2 r], [IC ms], [rho ms], [wall IC/rho], [band w=0], [w=0.1], [w=1], [in band], [w\*], [rank],
    ..d.pairs.map(x => ([#x.session], [#r(x.log2_r, 1)], [#r(x.ic_wall_ms, 2)], [#r(x.rho_wall_ms, 2)], [#r(x.wall_ratio, 3)],
      [#r(x.band_ratio_w0, 2)], [#r(x.band_ratio_w01, 2)], [#r(x.band_ratio_w1, 2)], [#x.in_band], [#r(x.w_star, 3)], [#x.rank])).flatten()
  )

  #d.sum.map(c => [#c.id: #r(c.value, 4) (#c.statement)]).join(linebreak())

  Fitted exponents $log_2 y = a + alpha log_2 r$ over the tier: #d.fits.map(f => [#f.arm #f.quantity α = #r(f.alpha, 3) (#f.points pts)]).join("; ").

  Pre-registered projection (wall-time exponents; extrapolated, not measured):

  #table(
    columns: 7,
    align: (right,) * 6 + (left,),
    [log2 r], [IC core-yr], [rho core-yr], [wall ratio], [band w=0], [w=1], [rank],
    ..d.proj.map(x => ([#x.log2_r], [#x.ic_core_years], [#x.rho_core_years], [#r(x.wall_ratio, 1)], [#r(x.band_ratio_w0, 1)], [#r(x.band_ratio_w1, 1)], [#x.rank])).flatten()
  )

  Post hoc projection (`projections_gae.csv`, addendum 1): fitted $P(w)$ exponents at the matched rho's ns per GAE at its largest measured size, log2 r = 128:

  #table(
    columns: 6,
    align: (right,) * 6,
    [w], [ns/GAE], [IC core-yr], [rho core-yr], [IC \$], [IC/rho],
    ..d.gae.filter(x => x.log2_r == "128").map(x => ([#x.w], [#x.ns_per_gae], [#x.ic_core_years], [#x.rho_core_years], [#x.ic_usd], [#r(x.ratio, 1)])).flatten()
  )

  #figure(image("fig_ic_" + d.name + ".svg", width: 100%), caption: [#d.name: measured and
  band-predicted wall time per solve (mean of 32), implied weight $w^*$, and the two
  projections (extrapolations, not measurements).])
]

== Deviations and missing data

- No PMU: branch-miss counts are missing; part B and addendum 1 test the branch
  hypothesis by substitution instead.
- F3 found `step_fp64` and `step_fp64bf` identical in the baseline binary, so
  F1's two arms ran the same instructions; addendum 1 (frozen before cm-1 and
  cm-2) added `fp64cm`.
- Part D's rho wall time carries a per-solve floor of about 3–4 ms on the smaller curves
  (the `rho ms` column); the post hoc GAE projection was added for that reason
  and is labelled post hoc.
- `run,batch,…` headers of `cm-pair/pairs.csv` name the compared fields
  (`fp64cm_ns`); the values were unchanged by that fix.

== Raw records

`runs/isa-*/` and `runs/cm-*/`: `measurements.json`, `samples.csv`,
`medians.csv`, `validation.{json,csv}`, `profile.json`, `fit.txt`,
`stderr.txt`, timestamps, `commit.txt`; `runs/isa-native-*/transfer/`;
`runs/compare/`; `runs/cm-pair/`; `runs/f3/`; `runs/gpu/`; `runs/ic-koblitz/`,
`runs/ic-prime/`. Graph sources: `../../plot.py` (layout only).
