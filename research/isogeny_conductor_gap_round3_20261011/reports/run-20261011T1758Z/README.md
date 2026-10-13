# Run report: E5 session `e5_combined_n17_t1` (2026-10-11T17:58:53Z to 2026-10-11T18:14:18Z UTC)

Descriptive record of one `ecbench` session.  Interpretation is in the
round's [README](../../README.md).  PDF: [`report.pdf`](report.pdf).

| item | value |
|:--|:--|
| lead | E5 of [`PROTOCOL.md`](../../PROTOCOL.md) |
| spec | [`specs/e5_combined_n17_t1.json`](../../specs/e5_combined_n17_t1.json) |
| command | `ecbench run --wait --spec specs/e5_combined_n17_t1.json --out sessions/e5_combined_n17_t1`, then `ecbench verify --replay 12` (`run_sessions.sh`) |
| session | `ECBS1he4c2c379c257`, git commit `2fa5132e8`, `git_dirty` false |
| host | Claude Code cloud container, Xeon @ 2.80 GHz; isolation earned per record: L1 on 6, L2 on 2 (the spec asked for L2) |
| workload | `icv1-f2m17-tm101-00378d4e` (a = 1 Koblitz curve over F_2^17, h = 2), target index 0 = point (0xba59, 0x8230) |
| arms | `rho.negation`; `ic.pipeline` with `binary-subspace:dimension=8` and `descent-algebraic:m=2` (blind), `descent-algebraic-fiber:fiber=lifts`, `descent-algebraic-fiber:fiber=combined`, all on `f4-f2` |
| runs | 8 requested, 8 completed, 8 verified (no warm-up, 2 measured rounds per arm) |
| audit | [`sessions/e5_combined_n17_t1.audit.json`](../../sessions/e5_combined_n17_t1.audit.json): ok, 8 of 8 verified, 8 replays reproduced, 0 problems |
| raw records | [`sessions/e5_combined_n17_t1/records.jsonl`](../../sessions/e5_combined_n17_t1/records.jsonl) |

**Measured, per run** (each round has its own algorithm seed, shared by the
arms in that round; solver operations are the F4 solver's own count, which
`S` leaves unpriced):

| arm | round | online wall (s) | F4 solver operations | solving degree, mean (max) | peak RSS |
|:--|--:|--:|--:|--:|--:|
| rho-neg | 0, 1 | 0.0, 0.0 | – | – | 14 MB |
| f4-blind | 0, 1 | 2.0, 1.3 | 164,172,136; 84,877,444 | 2.99, 2.99 (3) | 15 MB |
| f4-lifts | 0, 1 | 1.5, 1.5 | 114,545,148; 116,831,490 | 2.96, 2.95 (3) | 15 MB |
| f4-combined | 0, 1 | 453.5, 464.5 | 38,387,092,945; 34,853,985,256 | 3.06, 3.03 (4) | 0.65–0.99 GB |

**Derived:** combined ÷ blind solver operations, paired by round: 234× (round 0) and 411× (round 1).

Graph: [`e5_arms.svg`](e5_arms.svg) (source [`plot.py`](plot.py), recorded values only).
