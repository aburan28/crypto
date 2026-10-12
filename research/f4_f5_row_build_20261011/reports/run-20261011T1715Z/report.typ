
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

= Candidate Boolean F4/F5 allocation audit
UTC run id: `run-20261011T1715Z`. Native run and comparison: completed, exit status 0.

All 72 requested candidate records completed and matched baseline outputs and non-allocation counters exactly. Native comparison admitted the frozen allocation thresholds. Observed allocation request reductions: F4 0.25\%, F5 7.70\%, criterion 35.12\%. Cumulative requested-byte reductions: F4 42.47\%, F5 45.47\%, criterion 50.30\%. Percentages are native derived values, 100 x (1 - candidate/baseline), over the 24 calls per arm. Every workload met or improved both allocation measures.

#figure(image("allocation-audit.svg", width: 100%), caption: [#data.caption. Recorded counts and cumulative requested bytes; 24 calls per arm.])

#counts

Candidate kernel revision: 2e3dffda681a4771be511653180da8354de3bce3. Baseline production: e957e48c2b8cf1f84a7fc44a1d673103964a789d, with the same formatted diagnostic. Source/binary hashes: records/candidate-sha256.txt. Raw records: records/candidate.jsonl; matched comparison: records/comparison.json; per-system values: records/per-system-counts.csv.

The first library-test build failed on two pre-existing fixtures: an undefined parser fixture and a removed historical sorting utility. Both were repaired in test-only code while keeping all assertions. Failed and successful logs: `records/f4-f5-candidate-test-build*.log`. All 73 focused F4/F5 and Macaulay tests passed before this audit. No candidate audit records failed.

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
