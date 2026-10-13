# Koblitz higher-arity/rho benchmark copy

This directory is the benchmark-facing copy of the 2026-10-13 n=41 session.
The narrative, frozen protocol, PDF and probe stdout live in
[`research/koblitz_higher_arity_rho_20261013`](../../research/koblitz_higher_arity_rho_20261013/README.md).

- `summary.json`: machine-readable result and limitations.
- `records.jsonl`, `session.json`, `host.json`: complete ecbench measured session.
- `audit-replay-all.json`: 15/15 deterministic measured replays reproduced.
- `comparisons/`: saved paired operation and descriptive-wall comparisons.
- `*.isolation.jsonl`: same-cell B2 and B3 isolated probe conditions.

Counts are primary. The session earned L0 because PSI telemetry was unavailable;
all wall-time ratios are descriptive. All IC totals are lower bounds where the
records name unpriced native work.
