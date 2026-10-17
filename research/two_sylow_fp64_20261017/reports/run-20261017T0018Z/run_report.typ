// Run report (descriptive): two_sylow_search on fp_backend::Fp64, runs
// started 20261017T0018Z. Tables read the run records; compile with the study
// directory as root.
#set page(paper: "a4", margin: 1.6cm, numbering: "1")
#set text(size: 9.3pt)
#set heading(numbering: "1.")
#show table: set text(size: 7.6pt)

#let R = "../../runs/"
#let r(x, d) = str(calc.round(float(x), digits: d))
#let runs = ("base-1", "cand-1", "base-2", "cand-2")
#let meas = runs.map(n => json(R + n + "/measurements.json"))
#let txt(n, f) = read(R + n + "/" + f).trim()
#let cmp = csv(R + "compare/compare.csv", row-type: dictionary)
#let csum = json(R + "compare/summary.json")
#let samples = runs.map(n => csv(R + n + "/samples.csv", row-type: dictionary).filter(s => s.group == "tsy"))
#let st = meas.at(0).host.stable
#let pc = read(R + "protocol_commit.txt").trim()

= two_sylow_search on fp_backend::Fp64: run records (started 20261017T0018Z)

*Kind:* benchmark runs, descriptive (interpretation: `../../report.pdf`). \
*Study:* `research/two_sylow_fp64_20261017/`, protocol `PROTOCOL.md`, frozen at
#raw(pc.slice(0, 9)) before the first run. \
*Command:* `hwcost measure --fields fp64,w1 --two-sylow` per run
(#meas.at(0).params.reps repetitions, about #meas.at(0).params.target_ms ms
per arm); every run exited 0 and completed every arm.

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

Host: #st.cpu_model (signature #st.cpu_signature), #st.logical_cpus vCPUs,
#st.virtualization.pairs().map(((k, x)) => k + " = " + x).join(", "); isolation
L0; no PMU; no other build or test ran during the runs.

== Gate (G)

Each run's `gate.txt` (the `two_sylow` lines of its `stderr.txt`):

#raw(txt("base-1", "gate.txt"), block: true)

#raw(read(R + "compare/gate.txt").trim(), block: true)

== Claims (`hwcost compare`)

#table(
  columns: 5,
  align: (left, left, right, right, left),
  [id], [statement], [value], [threshold], [verdict],
  ..csum.claims.map(c => ([#c.id], [#c.statement], [#r(c.value, 4)], [#c.threshold], [#c.verdict])).flatten()
)

== two_sylow_search arms (ns per draw, medians of 7)

#table(
  columns: 10,
  align: (left,) + (right,) * 8 + (left,),
  [arm], [A1], [B1], [A2], [B2], [r1], [r2], [drift], [margin], [class],
  ..cmp.filter(x => x.arm.starts-with("tsy.")).map(x => (
    [#x.arm], [#r(x.a1, 0)], [#r(x.b1, 0)], [#r(x.a2, 0)], [#r(x.b2, 0)],
    [#r(x.r1, 3)], [#r(x.r2, 3)], [#r(x.drift, 3)], [#r(x.margin, 3)], [#x.class],
  )).flatten()
)

Every sample (ns per draw, repetition order):

#table(
  columns: 5,
  align: (left, left, left, left, left),
  [arm], ..runs.map(n => [#n]),
  ..cmp.filter(x => x.arm.starts-with("tsy.")).map(x => ([#x.arm], ..range(4).map(i => [#samples.at(i).filter(s => s.arm == x.arm).map(s => r(s.value, 0)).join(", ")]))).flatten()
)

#figure(image("fig_tsy.svg", width: 100%), caption: [Every `tsy` sample per run.])

== Controls and other arms

Arm classes (effect called only beyond max(2·drift, 2%) on the same side in
both pairs):

#table(
  columns: 3,
  align: (left, left, right),
  [group], [class], [arms],
  ..csum.classes.map(c => ([#c.group], [#c.class], [#c.arms])).flatten()
)

#figure(image("fig_ratios.svg", width: 100%), caption: [Ratios of the `tsy` arms and
of the kernel controls at B ≥ 4096.])

== Raw records

`runs/{base,cand}-{1,2}/` (`measurements.json`, `samples.csv`, `medians.csv`,
`stderr.txt`, `gate.txt`, `binary.txt`, timestamps), `runs/compare/`,
`runs/protocol_commit.txt`. Graph source: `../../plot.py` (layout only).
