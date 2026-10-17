# Koblitz higher-arity/rho benchmark copy

This directory is the benchmark-facing copy of the 2026-10-13 n=41 session.
The narrative, frozen protocol, PDF and probe stdout live in
[`research/koblitz_higher_arity_rho_20261013`](../../research/koblitz_higher_arity_rho_20261013/README.md).

- `summary.json`: machine-readable result and limitations.
- `records.jsonl`, `session.json`, `host.json`: complete ecbench measured session.
- `audit-replay-all.json`: 15/15 deterministic measured replays reproduced.
- `comparisons/`: saved paired operation and descriptive-wall comparisons.
- `*.isolation.jsonl`: same-cell B2 and B3 isolated probe conditions.
- `E2E_REPORT.pdf`, `n41-e2e-timing.svg`, `n41-e2e-timing.png`: five-round
  cold end-to-end elapsed-time follow-up and phase attribution.

Counts are primary. The session earned L0 because PSI telemetry was unavailable;
all wall-time ratios are descriptive. All IC totals are lower bounds where the
records name unpriced native work.

The follow-up session `ECBS1hb77d7527b579` measured five complete cold solves
per arm. B3 averaged 3.137 s versus B2's 22.536 s and strong rho's 13.005 ms;
relation extraction consumed 98.981% of B3's solve. The exact report and raw
records remain in the linked study directory.
