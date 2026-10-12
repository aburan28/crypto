# Candidate Boolean F4/F5 allocation audit

UTC run id: run-20261011T1715Z. Completed native audit; exit status 0.

All 72 requested candidate records completed and matched baseline outputs and non-allocation counters exactly. Native comparison admitted the frozen allocation thresholds. Observed allocation request reductions: F4 0.25%, F5 7.70%, criterion 35.12%. Cumulative requested-byte reductions: F4 42.47%, F5 45.47%, criterion 50.30%. Percentages are native derived values, 100 x (1 - candidate/baseline), over the 24 calls per arm. Every workload met or improved both allocation measures.

`target/release/examples/f4_f5_row_build_bench 3 > records/candidate.jsonl`\n\n`target/release/examples/f4_f5_row_build_bench compare records/baseline-aa.jsonl records/candidate.jsonl > records/comparison.json`

Candidate kernel revision: 2e3dffda681a4771be511653180da8354de3bce3. Baseline production: e957e48c2b8cf1f84a7fc44a1d673103964a789d, with the same formatted diagnostic. Source/binary hashes: records/candidate-sha256.txt. Raw records: records/candidate.jsonl; matched comparison: records/comparison.json; per-system values: records/per-system-counts.csv.

The first library-test build failed on two pre-existing fixtures: an undefined parser fixture and a removed historical sorting utility. Both were repaired in test-only code while keeping all assertions. Failed and successful logs: records/f4-f5-candidate-test-build*.log. All 73 focused F4/F5 and Macaulay tests passed before this audit. No candidate audit records failed.

Host: Linux x86-64, AMD EPYC 9V74, 5 logical CPUs, 18,882,699,264 bytes
RAM; Rust 1.99.0 (b940084d7), release build, one Rayon worker. No exclusive
CPU reservation was earned. No wall-time or peak-memory result is recorded.

Workloads: eight planted quadratic Boolean systems, (n,m,degree) =
(8,8,3), (8,8,4), (10,10,4), (10,10,5), each with seeds
0x2545f4914f6cdd1d and 0x9e3779b97f4a7c15. Three repetitions after a warmup.
Each arm has 24 complete calls. Allocation generation and output verification
are excluded from the counted interval. Requested bytes include every
allocation/reallocation request and are cumulative, not peak memory.

[Protocol](../../PROTOCOL.md), [host and raw records](../../records/),
[native driver](../../../../examples/f4_f5_row_build_bench.rs).
Graph source: [render_reports.py](../../render_reports.py) and
[graph-data.json](graph-data.json). [SVG](allocation-audit.svg);
[PDF](report.pdf); editable [Typst](report.typ).
