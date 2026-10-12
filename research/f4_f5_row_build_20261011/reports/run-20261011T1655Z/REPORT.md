# Baseline Boolean F4/F5 allocation audit

UTC run id: run-20261011T1655Z. Completed native audit; exit status 0.

The repaired, unoptimized baseline completed all 72 requested records. Exact F4 and F5 reduced rows agreed on every input; repeat correctness and counters matched. There was no baseline/candidate comparison in this initial run.

`/workspace/scratch/f4-f5-row-build-baseline 3 > records/baseline.jsonl`

Baseline production revision: e957e48c2b8cf1f84a7fc44a1d673103964a789d. Saved executable and source hashes: records/baseline-sha256.txt; raw values: records/baseline.jsonl. Baseline requests/bytes are shown; the second table arm is the later identical A/A control.

Prerequisite builds failed on duplicated source fragments at the original main revision c95a14de9bfa2711643a9318f3728d5d1109367b. The repair commits retain intact existing constructors, helpers, returns and assertions. Failed and successful build logs are retained as records/f4-f5-baseline*.log. These failures produced no valid benchmark values. The native comparison completed only after repairs. The formatted driver was subsequently independently calibrated against this baseline.

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
