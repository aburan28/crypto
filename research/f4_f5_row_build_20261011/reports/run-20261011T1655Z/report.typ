
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

= Baseline Boolean F4/F5 allocation audit
UTC run id: `run-20261011T1655Z`. Native run and comparison: completed, exit status 0.

The repaired, unoptimized baseline completed all 72 requested records. Exact F4 and F5 reduced rows agreed on every input; repeat correctness and counters matched. There was no baseline/candidate comparison in this initial run.

#figure(image("allocation-audit.svg", width: 100%), caption: [#data.caption. Recorded counts and cumulative requested bytes; 24 calls per arm.])

#counts

Baseline production revision: e957e48c2b8cf1f84a7fc44a1d673103964a789d. Saved executable and source hashes: records/baseline-sha256.txt; raw values: records/baseline.jsonl. Baseline requests/bytes are shown; the second table arm is the later identical A/A control.

Prerequisite builds failed on duplicated source fragments at the original main revision c95a14de9bfa2711643a9318f3728d5d1109367b. The repair commits retain intact existing constructors, helpers, returns and assertions. Failed and successful build logs are retained as `records/f4-f5-baseline*.log`. These failures produced no valid benchmark values. The native comparison completed only after repairs. The formatted driver was subsequently independently calibrated against this baseline.

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
