#set page(paper: "a4", margin: 1.6cm)
#set text(size: 9.2pt)
#set table(stroke: 0.4pt)
#show table: set text(size: 8pt)

#let S = json("summary_edwards.json")
#let sci(x, d: 2) = {
  if x == none { return "—" }
  let v = float(x)
  if v == 0 { return "0" }
  let e = calc.floor(calc.log(calc.abs(v), base: 10))
  let m = v / calc.pow(10.0, e)
  [#calc.round(m, digits: d)·10#super[#e]]
}

= Run report: $k = 5$ Edwards route with F4 trace replay, 2026-10-15T17:23Z onwards

Descriptive record of the runs; observations only. The interpretation is in `../../REPORT_EDWARDS.pdf` and `../../README.md`. This report was rendered from `summary_edwards.json` (a copy of the study's file, placed beside it).

*Spec.* `PROTOCOL_EDWARDS.md`, frozen at commit `3a8599b3f`; orchestration `run_edwards.sh`, then `run_edwards_resume.sh`; binary `examples/joux_vitse_edwards.rs`. Deviations are listed in `data/edwards/DEVIATION.md`:
+ the overlapped first launch was discarded;
+ the remaining arms were reordered behind the capped `p = 41` arm;
+ a container restart killed two arms (partial records kept);
+ a resumable driver was added, checked for equivalence on `p = 31`, seed 1, split at trial 12,000 (`data/edwards/resume_equivalence/CHECK.txt`).

*Host.* Claude Code cloud container, 4 vCPU Intel Xeon \@ 2.10 GHz, 15 GB, rustc 1.97.0. Isolation: none. Arms overlapped on the host at times.

*Measured.* The unit is $FF_p$ multiplications. These are counts, not wall time. Checks failed: *#S.checks_failed*.

#table(columns: 7, align: (left,) + (right,) * 6,
  [record], [$|F|$], [state], [trials], [relations], [IC online], [rho $FF_p$ mult. (ops)],
  ..S.dlp.map(d => (
    d.file, str(int(d.base_points)),
    if d.complete { if d.verified { "complete, verified" } else { "stopped (cap)" } } else { "incomplete" },
    if d.complete { str(int(d.trials)) } else { "—" },
    if d.complete { str(d.relations) } else { str(int(d.relations_recorded)) },
    if d.complete { sci(d.online_muls) } else { "—" },
    if d.complete { sci(d.rho_fp_muls) + " (" + str(d.rho_group_ops) + ")" } else { str(d.rho_group_ops) + " ops" },
  )).flatten()
)

#table(columns: 5, align: (left,) + (right,) * 4,
  [cost record], [residuals], [committed solve, mean], [replay, mean], [checks failed],
  ..S.cost.map(c => (c.file, str(c.residuals), sci(c.full_mean), sci(c.replay_mean), str(c.checks_failed))).flatten()
)

#figure(image("edwards_replay.svg", width: 78%), caption: [Cost of one test per residual, plotted from the cost records by `../../figures/plot_edwards.py`.])
