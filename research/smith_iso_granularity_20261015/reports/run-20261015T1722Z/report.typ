// Run report: values read from the recorded controls files.
#set page(paper: "a4", margin: 1.6cm)
#set text(size: 9pt)
#let c7 = json("/data/controls_p7.json")
#let c11 = json("/data/controls_p11.json")
#let s7 = json("/data/summary_p7.json")
#let s11 = json("/data/summary_p11.json")
= Run report: isomorphism-class census, p = 7 and p = 11
Study `research/smith_iso_granularity_20261015`, protocol frozen at `05f3d9203`. Host: EC2, Amazon Linux 2023,
2 vCPU, 7 GB, no isolation (L0); timings are descriptive only. Commands:
`smith_iso_granularity --p 7|11 --out data`, then `smith_iso_granularity summarize --dir data 7 11`.
#table(columns: 3, inset: 3pt, stroke: 0.4pt + luma(170),
  [*recorded*], [*p = 7*], [*p = 11*],
  [git head], [#c7.git_head.slice(0, 9)], [#c11.git_head.slice(0, 9)],
  [octics analysed], [#c7.octics], [#c11.octics],
  [models / iso classes / L-classes], [#c7.models / #c7.iso_classes / #c7.l_classes], [#c11.models / #c11.iso_classes / #c11.l_classes],
  [C1 / C2 / C3], [#c7.c1_pass / #c7.c2_pass / #c7.c3_pass], [#c11.c1_pass / #c11.c2_pass / #c11.c3_pass],
  [seconds (analysis / total), threads], [#calc.round(c7.seconds_analysis, digits: 1) / #calc.round(c7.seconds_total, digits: 1), #c7.threads],
  [#calc.round(c11.seconds_analysis, digits: 1) / #calc.round(c11.seconds_total, digits: 1), #c11.threads],
  [observed C; iso-null 2.5–97.5%], [#calc.round(s7.C, digits: 4); #calc.round(s7.null_primary.lo, digits: 4)–#calc.round(s7.null_primary.hi, digits: 4)],
  [#calc.round(s11.C, digits: 4); #calc.round(s11.null_primary.lo, digits: 4)–#calc.round(s11.null_primary.hi, digits: 4)],
)
#figure(image("/figures/iso_granularity.svg", width: 100%), caption: [Recorded values of this run; see `../../REPORT.pdf`.])
No runs failed or timed out. Raw records: `data/iso_p{7,11}.csv.gz`, `data/controls_p*.json`, `data/summary_p*.json`.
