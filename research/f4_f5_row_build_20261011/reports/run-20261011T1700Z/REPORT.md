# Baseline A/A driver calibration

UTC run id: run-20261011T1700Z. Completed native audit; exit status 0.

All 72 requested A/A records completed. The formatted final driver, linked to the saved unoptimized baseline library, reproduced every baseline output, metadata field, allocation count and byte count exactly. The comparison's optimization acceptance flag is false because identical baseline copies cannot satisfy improvement thresholds; calibration passed.

`/workspace/scratch/f4-f5-row-build-baseline-formatted 3 > records/baseline-aa.jsonl`\n\n`/workspace/scratch/f4-f5-row-build-baseline-formatted compare records/baseline.jsonl records/baseline-aa.jsonl > records/baseline-aa-comparison.json`

Production baseline e957e48c2b8cf1f84a7fc44a1d673103964a789d; formatted diagnostic in ee344374. Original and formatted binaries were retained under /workspace/scratch. Hashes: records/baseline-aa-sha256.txt. Raw values and control: records/baseline.jsonl, records/baseline-aa.jsonl, records/baseline-aa-comparison.json.

The driver was compiled natively using rustc --edition=2021, -C opt-level=3, the saved baseline rlib and the original release serde_json, rayon and blake3 dependencies. Full reconstruction command is in the study report. Formatting did not alter the allocation accounting interval. No unsuccessful A/A records occurred.

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
