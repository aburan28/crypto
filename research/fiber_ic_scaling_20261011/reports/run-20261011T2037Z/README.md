# Run report: `scaling_h4` (2026-10-11T20:37:11Z to 2026-10-11T20:39:14Z UTC)

Descriptive record of one `ecbench` session.  Interpretation is in the
study [README](../../README.md).  PDF: [`report.pdf`](report.pdf).

| item | value |
|:--|:--|
| protocol | [`PROTOCOL.md`](../../PROTOCOL.md) |
| spec | [`specs/scaling_h4.json`](../../specs/scaling_h4.json) |
| command | `ecbench run --wait`, then `ecbench verify --replay 12` ([`run_sessions.sh`](../../run_sessions.sh)) |
| session | `ECBS1h4a5009a8f716`, git commit `0314b6a7a`, `git_dirty` false |
| host | Claude Code cloud container, Xeon @ 2.80 GHz; isolation earned per record: L2 on 231, L1 on 105 (the spec asked for L2) |
| runs | 336 requested, 336 completed, 336 verified |
| audit | [`sessions/scaling_h4.audit.json`](../../sessions/scaling_h4.audit.json): ok, 12 replays reproduced, 0 problems |
| records | [`sessions/scaling_h4/records.jsonl`](../../sessions/scaling_h4/records.jsonl) |
| per-curve means | [`TABLES.md`](../../TABLES.md) section `scaling_h4` (from `examples/fiber_scaling_tables.rs`) |

The measured values are each arm's `S` (`ecbench.gae` cost over the
floor), per curve, over 6 runs (3 targets × 2 rounds).  The per-curve means
are in `TABLES.md`.  The graph [`s_by_arm.svg`](s_by_arm.svg) plots them
against `log₂ #E`, with min–max bars; its source is [`plot.py`](plot.py).
No wall times are reported.
