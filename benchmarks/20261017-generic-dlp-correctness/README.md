# 2026-10-17: generic-DLP correctness fixes (collab 2-cycle escape, Pohlig–Hellman)

These are copies for browsing. The source of truth is
[`research/generic_dlp_correctness_20261017`](../../research/generic_dlp_correctness_20261017/README.md),
which holds:

- the protocol;
- the frozen inputs;
- the raw counts and pinned-test transcripts;
- `summary.json` (written natively by
  `examples/generic_dlp_correctness_summary.rs`);
- the plotting source `scripts/plot.py`;
- the report source.

The PDF and SVG here are byte-identical copies. The PNG is the 150-dpi
render made from the same plotting source as the SVG.

- [`REPORT.pdf`](REPORT.pdf) is the results report. Its source is
  `research/generic_dlp_correctness_20261017/REPORT.typ`.
- `fig_correctness` (SVG plus PNG) has two panels:
  - collab dead trails per curve and mode, before and after;
  - Pohlig–Hellman `total_steps` per solve, before against after.

**What is measured:** deterministic counts only; nothing is timed.

**Class:** a correctness fix (collab) and an algorithmic change to counted
cost (Pohlig–Hellman). `S` and every IC/rho comparison are unchanged.

**Host:** counts do not depend on the host. Both arms were built with rustc
1.97.0, `cargo build --release`.
