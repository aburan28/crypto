// Run report: values read from the recorded controls and summary files.
#set page(paper: "a4", margin: 1.6cm)
#set text(size: 9pt)
#let c7 = json("/data/controls_p7.json")
#let c11 = json("/data/controls_p11.json")
#let s7 = json("/data/summary_p7.json")
#let s11 = json("/data/summary_p11.json")
#let r(x) = str(calc.round(x, digits: 4))
= Run report: class invariants, p = 7 and p = 11
Study `research/smith_class_invariant_20261015`. Code at #raw(c7.git_head.slice(0, 9)) (both primes). Host: EC2,
Amazon Linux 2023, 2 vCPU, 7 GB, no isolation (L0); timings are descriptive only. Commands:
`smith_class_invariant compute --p 7|11 --out data`, then `smith_class_invariant summarize --dir data 7 11`.
Started 2026-10-15T23:32Z. Inputs: frozen `research/smith_iso_granularity_20261015/data/iso_p{7,11}.csv.gz`.
#table(columns: 3, inset: 3pt, stroke: 0.4pt + luma(170),
  [*recorded*], [*p = 7*], [*p = 11*],
  [git head], [#c7.git_head.slice(0, 9)], [#c11.git_head.slice(0, 9)],
  [isomorphism classes / tractable (computed)], [#c7.iso_classes / #c7.tractable_classes], [#c11.iso_classes / #c11.tractable_classes],
  [X2 uncertified classes], [#c7.counts.k3_uncertified], [#c11.counts.k3_uncertified],
  [K2 / K3 / K4 / K4′ / K5 / K5c],
  [#c7.k2_pass_j4 #c7.k2_pass_sylow / #c7.k3_pass / #c7.k4_pass / #c7.k4c_pass / #c7.k5_pass / #c7.k5c_pass],
  [#c11.k2_pass_j4 #c11.k2_pass_sylow / #c11.k3_pass / #c11.k4_pass / #c11.k4c_pass / #c11.k5_pass / #c11.k5c_pass],
  [K4 mismatching classes (frozen descent)], [#c7.counts.k4_mismatch], [#c11.counts.k4_mismatch],
  [seconds, threads], [#calc.round(c7.seconds, digits: 0), #c7.threads], [#calc.round(c11.seconds, digits: 0), #c11.threads],
  [observed C; unconditioned null 2.5–97.5%], [#r(s7.C); #r(s7.null_unconditioned.lo)–#r(s7.null_unconditioned.hi)],
  [#r(s11.C); #r(s11.null_unconditioned.lo)–#r(s11.null_unconditioned.hi)],
)
#figure(image("/figures/class_invariant_nulls.svg", width: 100%), caption: [Recorded values of this run; see `../../REPORT.pdf`.])
Both compute runs completed; none failed or timed out. K4 and K5 fail for the frozen descent (X3, X3+), as recorded
in Amendment 1; these candidates are reported invalid. Raw records: `data/inv_p{7,11}.csv.gz`,
`data/controls_p*.json`, `data/summary_p*.json`, `data/run_p*.log`.
