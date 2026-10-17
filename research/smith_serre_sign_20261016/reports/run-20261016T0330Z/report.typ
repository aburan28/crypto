// Run report: values read from the recorded files.
#set page(paper: "a4", margin: 1.6cm)
#set text(size: 9pt)
#let c7 = json("/data/controls_p7.json")
#let c11 = json("/data/controls_p11.json")
= Run report: Prop. 4 sign census, p = 7 and p = 11
Study `research/smith_serre_sign_20261016`, protocol frozen at `3be242141`. Host: EC2, Amazon Linux 2023, 2 vCPU,
7 GB, no isolation (L0); timings descriptive. Commands: `smith_iso_granularity --p 7|11 --out data --subgroups`, then
`smith_iso_granularity summarize-serre --dir data 7 11`.
#table(columns: 3, inset: 3pt, stroke: 0.4pt + luma(170),
  [*recorded*], [*p = 7*], [*p = 11*],
  [git head], [#c7.git_head.slice(0, 9)], [#c11.git_head.slice(0, 9)],
  [octics / iso classes / L-classes], [#c7.octics / #c7.iso_classes / #c7.l_classes], [#c11.octics / #c11.iso_classes / #c11.l_classes],
  [C1 / C2 / C3 / C4], [#c7.c1_pass / #c7.c2_pass / #c7.c3_pass / #c7.c4_pass], [#c11.c1_pass / #c11.c2_pass / #c11.c3_pass / #c11.c4_pass],
  [C5 violations (degenerate, expected)], [#c7.c5_twist_complementarity_violations], [#c11.c5_twist_complementarity_violations],
  [seconds total, threads], [#calc.round(c7.seconds_total, digits: 1), #c7.threads], [#calc.round(c11.seconds_total, digits: 1), #c11.threads],
)
#figure(image("/figures/serre_sign.svg", width: 100%), caption: [Recorded values; see `../../REPORT.pdf`.])
No runs failed or timed out. Raw records: `data/iso_p*.csv.gz`, `data/iso_sub_p*.csv.gz`, `data/controls_p*.json`,
`data/serre_p*.json`, `data/serre_verdict.json`.
