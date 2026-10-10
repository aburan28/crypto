# First frozen n19 sparse S3 growth attempt

**Decision: CENSORED at the independent verifier's 600 s child wall cap.** The first and only release of freeze SHA-256 `15f38ea4eb18244234fa86c497ef521805ce2f9877a806936796d3d7d2642fe3` was Actions run [36526326477](https://github.com/aburan28/crypto/actions/runs/36526326477), PR #795 head `c3422a179522db64c1f49ff396431fe31da2bcd5`. Its unique label event ID was `32052162480`. The exact 12-file runner archive is in [evidence/final](evidence/final), with receipt SHA-256 `e5de2330e8fa424e095117f700e3a328763b943453bb2e1c1da75cb23c4fa839` and manifest SHA-256 `8e8dba393698372760b448ac041dff50bf011851e04cc1f4f2f85156abc40fd9`. `ci_replay.py --evidence evidence/final` rehashes and admits the archive as `CENSORED`; the labeled job's archive-audit and upload steps succeeded. No second run of this freeze is allowed.

The exporter completed in 68.30 s outer wall, 68.13 s reported child wall, and 141,070,336 B direct-child peak RSS. It emitted a 3,235,040-byte DIMACS file with 97,771 variables, 179,668 clauses, and 58,825 primary paths. These are **producer observations**, not independently certified CNF semantics. All five reported prefix groups stayed within the preregistered state, pair, path, DIMACS and RSS caps:

| Prefix | Admitted states | Cumulative transition input pairs |
| --- | ---: | ---: |
| S2 | 25 | 16 |
| S3 | 167 | 116 |
| S4 | 1,085 | 784 |
| S5 | 6,976 | 5,124 |
| S6 | 40,614 | 33,028 |

The independent verifier ran for 600.15 s outer wall, 599.92 s CPU, and 102,875,136 B direct-child peak RSS. It wrote `CENSORED: independent verifier wall cap`, exiting 1 before target/path/negative-control admission. Its traceback is inside the inherited `check_onehot_block` check, which first compares the exact Sinz clause block and then tests a satisfying auxiliary assignment for *every* possible chosen primary literal when the group has more than five members. That witness sweep is quadratic in a group's size despite the block already having been checked against the exact expected clauses. The archive does not identify a failed clause, wrong group law, or invalid target. It also cannot establish that the exporter semantics pass.

The next falsifiable step is a **new freeze**, not a retry: keep these factor/target inputs and the archived CNF, replace the large-group witness sweep with a separately proved linear structural check of the Sinz block, retain exhaustive truth-table controls for small groups, add a large-group corruption control, and independently replay the archived path/point-target semantics under a new verifier wall cap chosen before observation. Report a verifier speed and semantic result separately from exporter performance. Only after a complete independent pass should this n19 layout enter a frozen SAT-solver comparison. No n131 scaling, PDP success, ECDLP, or matched-rho speed claim follows from this censored attempt.
