
#set page(paper: "a4", margin: 15mm)
#set text(font: "DejaVu Sans", size: 9pt)
#set par(justify: false, leading: 0.65em)
#set heading(numbering: none)
#let data = json("graph-data.json")
#let arms = data.arms
#let counts = table(
  columns: (1.3fr, 1fr, 1fr, 1fr, 1fr),
  inset: 4pt,
  align: (left, right, right, right, right),
  table.header([Arm], [Baseline requests], [Other requests], [Baseline bytes], [Other bytes]),
  ..arms.map(a => (
    [#a.arm], [#a.baseline_allocations], [#a.candidate_allocations],
    [#a.baseline_allocated_bytes], [#a.candidate_allocated_bytes],
  )).flatten()
)

= Five alternating allocation-audit replications
UTC run id: `run-20261011T1716Z`. Native run and comparison: completed, exit status 0.

Five complete baseline/candidate pairs completed in alternating order. Each pair contains 72 baseline and 72 candidate records with three repetitions per workload. Native comparison passed every pair; byte-for-byte cmp reproduced both primary JSONL files in every round. The graph displays the recorded values for one pair, identical across all five. No timings, pooled wall-time values or confidence intervals were computed.

#figure(image("allocation-audit.svg", width: 100%), caption: [#data.caption. Recorded counts and cumulative requested bytes; 24 calls per arm.])

#counts

Run start: records/replication-start.txt. Same revisions and binary hashes as the candidate run. Exact records and comparisons: records/paired/. All five baseline and five candidate files reproduced the primary records exactly; native comparison exited 0 for all pairs.

No failures, timeouts or allocation regressions occurred in this replication run. Alternation is an allocation reproducibility control. CPU isolation was not established, so no runtime inference follows. There are 24 calls per arm per variant per pair; the graph does not sum across pairs.

*Host and scope.* Linux x86-64, AMD EPYC 9V74, 5 logical CPUs,
18,882,699,264 bytes RAM; Rust 1.99.0 (b940084d7), release, one Rayon worker.
Exclusive CPU isolation was not established. Wall time and peak memory were
not measured. This is a standalone Boolean matrix allocation diagnostic;
full ECDLP cost, matched rho ratios, m=83 scaling and general-field gains remain
unestablished.

*Inputs and accounting.* Eight planted quadratic systems: (n,m,d) =
(8,8,3), (8,8,4), (10,10,4), (10,10,5), seeds
`0x2545f4914f6cdd1d` and `0x9e3779b97f4a7c15`. Warmup precedes three repetitions.
Each bar covers 24 calls. Fixture generation, hashing and verification are
outside the audit; packing, elimination and output construction are included.
Bytes are cumulative requests, not peak memory. Raw repeated values are
retained; no statistical interval was estimated.

*Evidence.* Repository `aburan28/crypto`; study
`research/f4_f5_row_build_20261011`. Native driver:
`examples/f4_f5_row_build_bench.rs`. Full commands and paths: this directory's
`REPORT.md`; graph values: `graph-data.json`; rendering source:
`../../render_reports.py`; raw records and hashes: `../../records/`.
