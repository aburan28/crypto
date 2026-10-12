# Run report: E5 session `e5_combined_n17_t0` (2026-10-11T17:27:31Z to 2026-10-11T17:43:15Z UTC)

Descriptive record of one `ecbench` session.  Interpretation is in the
round's [README](../../README.md).  PDF: [`report.pdf`](report.pdf).

| item | value |
|:--|:--|
| lead | E5 of [`PROTOCOL.md`](../../PROTOCOL.md) |
| spec | [`specs/e5_combined_n17_t0.json`](../../specs/e5_combined_n17_t0.json) |
| command | `ecbench run --wait --spec specs/e5_combined_n17_t0.json --out sessions/e5_combined_n17_t0`, then `ecbench verify --replay 12` (`run_sessions.sh`) |
| session | `ECBS1h05f404b1f89b`, git commit `1a63ec419`, `git_dirty` false |
| host | Claude Code cloud container, Xeon @ 2.80 GHz; isolation earned per record: L1 on 6, L2 on 2 (the spec asked for L2) |
| workload | `icv1-f2m17-tm101-00378d4e` (a = 1 Koblitz curve over F_2^17, h = 2), target index 0 = point (0x984b, 0x1c882) |
| arms | `rho.negation`; `ic.pipeline` with `binary-subspace:dimension=8` and `descent-algebraic:m=2` (blind), `descent-algebraic-fiber:fiber=lifts`, `descent-algebraic-fiber:fiber=combined`, all on `f4-f2` |
| runs | 8 requested, 8 completed, 8 verified (no warm-up, 2 measured rounds per arm) |
| audit | [`sessions/e5_combined_n17_t0.audit.json`](../../sessions/e5_combined_n17_t0.audit.json): ok, 8 of 8 verified, 8 replays reproduced, 0 problems |
| raw records | [`sessions/e5_combined_n17_t0/records.jsonl`](../../sessions/e5_combined_n17_t0/records.jsonl) |

**Measured, per run** (each round has its own algorithm seed, shared by the
arms in that round; solver operations are the F4 solver's own count, which
`S` leaves unpriced):

| arm | round | online wall (s) | F4 solver operations | solving degree, mean (max) | peak RSS |
|:--|--:|--:|--:|--:|--:|
| rho-neg | 0, 1 | 0.0, 0.0 | – | – | 14 MB |
| f4-blind | 0, 1 | 0.9, 1.8 | 65,072,207; 133,408,080 | 2.98, 2.97 (3) | 15 MB |
| f4-lifts | 0, 1 | 1.7, 1.3 | 106,694,797; 104,216,038 | 2.96, 2.97 (3) | 15 MB |
| f4-combined | 0, 1 | 324.5, 613.2 | 32,012,560,285; 70,613,284,746 | 3.04, 3.04 (4) | 0.64–0.66 GB |

**Derived:** combined ÷ blind solver operations, paired by round: 492× (round 0) and 529× (round 1).

Graph: [`e5_arms.svg`](e5_arms.svg) (source [`plot.py`](plot.py), recorded values only).
