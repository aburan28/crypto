# Five alternating allocation-audit replications

UTC run id: run-20261011T1716Z. Completed native audit; exit status 0.

Five complete baseline/candidate pairs completed in alternating order. Each pair contains 72 baseline and 72 candidate records with three repetitions per workload. Native comparison passed every pair; byte-for-byte cmp reproduced both primary JSONL files in every round. The graph displays the recorded values for one pair, identical across all five. No timings, pooled wall-time values or confidence intervals were computed.

For `audit_round` from 1 through 5, run the saved formatted baseline with `3`, then the candidate with `3`, then native `compare`. Retain each baseline-${audit_round}.jsonl, candidate-${audit_round}.jsonl and comparison-${audit_round}.json under records/paired/. Compare every raw JSONL with its primary counterpart using `cmp`.

Run start: records/replication-start.txt. Same revisions and binary hashes as the candidate run. Exact records and comparisons: records/paired/. All five baseline and five candidate files reproduced the primary records exactly; native comparison exited 0 for all pairs.

No failures, timeouts or allocation regressions occurred in this replication run. Alternation is an allocation reproducibility control. CPU isolation was not established, so no runtime inference follows. There are 24 calls per arm per variant per pair; the graph does not sum across pairs.

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
