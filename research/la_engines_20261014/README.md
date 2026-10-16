# IC linear-algebra engines: ours against CADO-NFS, and a port of its filter

2026-10-14/15. Results report: [`REPORT.pdf`](REPORT.pdf) (source `REPORT.typ`);
run report: [`reports/run-20261014T1728Z/report.pdf`](reports/run-20261014T1728Z/report.pdf).
Protocol: [`PROTOCOL.md`](PROTOCOL.md), frozen before data, and amendments
[1](PROTOCOL_AMENDMENT_1.md)–[5](PROTOCOL_AMENDMENT_5.md), each committed before its own data.

**Question.** Is the last step of our index calculus, the linear algebra, as
fast as it can be, and does it limit end-to-end solves?

**Findings** (host: 4 vCPU Xeon @ 2.80 GHz, shared, isolation L0; synthetic
m = 3 relation systems over a 62-bit prime; times are `solve_s`, relations in,
verified logarithms out; ratios within one session, per-seed range in brackets):

- End to end the linear algebra is not the bottleneck: 0.11–1.9% of real
  `ic logs` precomputes (27 runs); the largest real IC matrix in the
  repository (N73, 600 columns) took 5.1 h to collect and milliseconds to solve.
- Our block Wiedemann beats CADO-NFS `bwc` on the identical matrix at every
  size measured (1.86× at 131,072 columns, more at smaller sizes).
- CADO-NFS's filter left a 5–8× smaller matrix than ours; this study ports it
  (`FilterStrategy::Cado`: Ω₃₁ clique purge, Markowitz/spanning-tree merge) and
  speeds it up. At D = 262,144: 86.57 s → 8.35 s (**10.40×** [9.99, 11.06]),
  and **1.93×** [1.74, 1.96] faster than CADO-NFS's pipeline on the same
  matrices (7.83 s against 14.73 s).
- Below about 4,096 columns the new filter is slower (0.40× at 1,024; both in
  milliseconds). `ecbench`'s default `IncrementalGauss` grows 6–9× per doubling
  (19.1 s at 8,192 columns).

**Code.** `src/cryptanalysis/koblitz_sparse_la.rs` (`FilterStrategy`,
`FilterOptions::cado()`, faster duplicate removal and purge),
`src/bin/ic/experiment.rs` (`ic logs --filter-strategy cado|legacy`, default
`cado`), `examples/la_engine_bench.rs` (the instrument).

**Layout.** `raw/` run records (`*.jsonl`) and CADO logs; `scripts/` thin
orchestration (`grid*.sh`, `run_cado*.sh`, `e2e.sh`, `rerun_contended.sh`) and
`plots.py` (layout only); `figures/` SVG + 150 dpi PNG; `summary.json` medians
and ranges behind every figure.

**Reproduce.** Build CADO-NFS `692ecb7e` with
`-DBWC_GFP_ARITHMETIC_BACKENDS="p1;p3"` and the targets `bwc_full_gfp mf_scan2
mf_bal lingen_p1 lingen_p3 purge merge-dl replay-dl`; `cargo build --release
--example la_engine_bench --bin ic`; run the `scripts/grid*.sh` in order with
`TMPDIR` on a short path (CADO's `bwc.pl` mishandles very long ones).
