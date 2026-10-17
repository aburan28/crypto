# Third paired run: row-width mismatch

Started 2026-10-13 01:59 UTC (2026-10-12 18:59 PDT). Frozen
[protocol](../../PROTOCOL.md), commit `08276ed63638`; the source was the
v2 runner with the unmeasured lookup value removed while its header still
named that column. Darwin arm64, shared host, L0.

Command: `target/release/examples/isogeny_fiber_pair_20261012 > raw-v3.tsv 2> run-v3.log`.
The program logged 48/48 verified arms and both map checks, but the
[raw TSV](../../raw-v3.tsv) has 47 header fields and 48 data fields. Its
SHA-256 is `77b7fa9fe07e71b74c7df462577252be49b94e31d6dd09a910e80bea448f8d7a`.
The [map log](../../run-v3.log), [replay](../../replay-v3.tsv), and
[replay log](../../replay-v3.log) are preserved; replay bytes match.

This is an output-schema failure, so no graph or numeric comparison is
derived from it. The [flow diagram](flow.svg) and [PDF](report.pdf) show
the completed calculation and rejected record. The v4 runner asserts
header and row width before emitting each record.
