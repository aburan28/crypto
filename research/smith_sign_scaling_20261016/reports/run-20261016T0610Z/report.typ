// Run report: values read from the recorded files.
#set page(paper: "a4", margin: 1.6cm)
#set text(size: 9pt)
#let R = (101, 1009, 10007).map(p => json("/data/scaling_p" + str(p) + ".json"))
= Run report: Prop. 4 sign scaling, p = 101, 1009, 10007
Protocol frozen at `5fd6133d0`. Host: EC2, Amazon Linux 2023, 2 vCPU, 7 GB, no isolation (L0); timings descriptive.
Command per prime: `smith_sign_scaling run --p P --n 1000000 --out data`; then `verdict --dir data 101 1009 10007`.
#table(columns: 7, inset: 3pt, stroke: 0.4pt + luma(170),
  [*p*], [*seed*], [*git head*], [*drawn / squarefree*], [*β ± SE*], [*valid*], [*seconds, threads*],
  ..R.map(r => ([#r.p], [#r.seed], [#r.git_head.slice(0, 9)], [#r.n_drawn / #r.squarefree],
    [#calc.round(r.beta, digits: 3) ± #calc.round(r.beta_se, digits: 3)], [#r.valid],
    [#calc.round(r.seconds, digits: 0), #r.threads])).flatten(),
)
#figure(image("/figures/sign_scaling.svg", width: 100%), caption: [Recorded values; see `../../REPORT.pdf`.])
No runs failed or timed out. Raw records: `data/scaling_p*.json`, `data/verdict.json`, `data/run.log`.
