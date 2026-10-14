# Run report: `scaling_h72_small_b` (2026-10-11T20:44:14Z to 2026-10-11T20:44:48Z UTC)

Descriptive record of one `ecbench` session.  Interpretation is in the
study [README](../../README.md).  PDF: [`report.pdf`](report.pdf).

| item | value |
|:--|:--|
| protocol | [`PROTOCOL.md`](../../PROTOCOL.md) |
| spec | [`specs/scaling_h72_small_b.json`](../../specs/scaling_h72_small_b.json) |
| command | `ecbench run --wait`, then `ecbench verify --replay 12` ([`run_sessions.sh`](../../run_sessions.sh)) |
| session | `ECBS1h54835df66b11`, git commit `5c8a7694f`, `git_dirty` false |
| host | Claude Code cloud container, Xeon @ 2.80 GHz; isolation earned per record: L2 on 88, L1 on 56 (the spec asked for L2) |
| runs | 144 requested, 144 completed, 144 verified |
| audit | [`sessions/scaling_h72_small_b.audit.json`](../../sessions/scaling_h72_small_b.audit.json): ok, 12 replays reproduced, 0 problems |
| records | [`sessions/scaling_h72_small_b/records.jsonl`](../../sessions/scaling_h72_small_b/records.jsonl) |
| per-curve means | [`TABLES.md`](../../TABLES.md) section `scaling_h72_small_b` (from `examples/fiber_scaling_tables.rs`) |

The measured values are each arm's `S` (`ecbench.gae` cost over the
floor), per curve, over 6 runs (3 targets × 2 rounds).  The per-curve means
are in `TABLES.md`.  The graph [`s_by_arm.svg`](s_by_arm.svg) plots them
against `log₂ #E`, with min–max bars; its source is [`plot.py`](plot.py).
No wall times are reported.
