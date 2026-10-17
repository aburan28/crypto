# f4_fp: blocked Faugère–Lachartre pair reduction (2026-10-16)

A callgrind profile of `f4_fp` on `quad_n10_p65521` (1 thread) put 63 % of
instructions in `add_lazy_avx2`, mostly the Faugère–Lachartre step that
computes `D − C·B'` for the pair rows one row at a time, re-reading each
needed `B'` row from memory for every pair row. The change
(`src/cryptanalysis/f4_fp.rs`, `reduce_pairs_blocked`) takes 64 pair rows
at a time, groups their updates by reducer and walks the free columns in
tiles of 256, so a tile of a `B'` row is read once for every row of the
block while the block's accumulators stay in cache. The updates and their
exact sum are unchanged; the reduced bases are identical.
`F4_FP_FL_BLOCKED=0` restores the row-by-row product (same-binary control).

Run report: [`reports/run-20261016T1720Z/report.pdf`](reports/run-20261016T1720Z/report.pdf)
(source `report.typ`, records in `raw/`).

Result, 1 thread, median of 3 interleaved rounds on a shared 4 vCPU Xeon
(L0): blocked/control 0.78× on `quad_n10` and 0.89× on `quad_n11` (per-round
ranges disjoint), 0.88× on `quad_n9` (ranges overlap), within noise on the
other 15 cases. Against msolve 0.6.5 (`-t 1`, process wall time) in the same
rounds: 0.78× / 0.75× on `quad_n10` / `quad_n11`, 1.11–1.49× on `quad_n8`,
`quad_n9`, `cubic_n5`, `cubic_n6`. Not established: any effect at 4
threads (exploratory runs showed none), or on other hosts.

Reproduce: build `cargo build --release --example gb_cas_export`, export the
systems with `GB_CAS_LARGE=1 gb_cas_export fp <dir> 1`, then
`reports/run-20261016T1720Z/raw/run_ab.sh` (expects `gce_blk` = that binary
and `inL/` = `<dir>` beside it). Plot: `python plot.py raw/ab_full.tsv <out>`.
