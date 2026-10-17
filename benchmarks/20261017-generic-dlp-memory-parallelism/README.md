# 2026-10-17: BSGS table prefetch and inversion-free rho start-up

These are copies for browsing. The source of truth is
[`research/generic_dlp_memory_parallelism_20261017`](../../research/generic_dlp_memory_parallelism_20261017/README.md),
which holds:

- the protocol;
- the frozen inputs;
- the raw records and `summary.json` (written natively by
  `examples/generic_dlp_analyze.rs`);
- the plotting source `scripts/plot.py`;
- the report sources.

The PDFs and SVGs here are byte-identical copies. Each PNG is the 150-dpi
render made from the same plotting source as its SVG.

- [`REPORT.pdf`](REPORT.pdf) is the results report. Its source is
  `research/generic_dlp_memory_parallelism_20261017/REPORT.typ`.
- [`run-report.pdf`](run-report.pdf) is the run report. Its source is
  `research/generic_dlp_memory_parallelism_20261017/reports/run-20261017T0201Z/report.typ`.
- The figures, each as SVG plus a 150-dpi PNG:
  - `fig_bsgs`: BSGS one-target medians, and A/B against the A/A′ noise
    band;
  - `fig_rho`: the rho solve and its start-up, with A/B at 1 and 3 threads;
  - `fig_external`: both native arms against PARI/GP 2.15.4 `elllog` on the
    same curves and targets.

**Class:** engineering. Every counted output is identical between the arms,
and `S` is unchanged.

**Host:** an Intel Xeon (Sapphire/Emerald Rapids class) at 2.1 GHz, a 4-vCPU
KVM guest, isolation L0. 26 of 438 timed runs were contended and are left
out.
