# IC basis and setup study

Status: training and fresh cold holdouts are independently verified. The final
cohort stage is RUNNING under the native CPU3 reservation. Full-study statistics,
comprehensive PDF/evidence publication and review gates remain pending.
[Draft study PR1899](https://github.com/aburan28/crypto/pull/1899).

| Authorized requirement | Actual status and evidence |
|---|---|
| Relation/rank/log timing, probes, rank gain and memory | Native streaming instrumentation; all8completed training rank/log transcripts independently verified. [Profiles](reports/run-20261015T1516Z/native_profiles.json) and [report](reports/run-20261015T1516Z/REPORT.md). |
| Structurally different bases including dense against K600 | Five candidates frozen before timing. Native audit verified exact labels, prefix/window recipes and disjoint C2 support. All10IC training processes attempted;8completed,2C4allocation failures retained. C3 selected. |
| Independently verified rank/log transcript |8native C++ full-row/rank/backsolve/representative/witness/answer audits PASS; external audit resources separately metered. [Receipts](records/training-audit/). All16 cold IC rank/log audits also passed; cohort audits pending. |
| Complete matched costs and independent preparations | Training18processes attempted/16complete; three A/A work pairs verified. C3 selection committed before8freshholdouts. Remaining37processes/100plannedanswers frozen;26-process cold stage completed and independently verified;11-process cohort stage running. |
| Reports, graphs, PDFs and evidence publication | Inspected4pagetrainingPDF+SVG/150dpiPNG published. Inspected4page cold-stage PDF published too. Final addendum, updated comprehensive PDF, complete-study evidenceZIP/canonical diagnostics pending. |

[Prospective protocol](PROTOCOL.md), [frozen training selection](selection/),
[training PDF](reports/run-20261015T1516Z/n83_basis_training_report.pdf),
[holdout continuation](HOLDOUT_CONTINUE.md).

Training mean actual process costs: C0 1250.669494426s, C1 1157.2533628935s,
C2 1205.620140093s, selected C3 622.6626231075s. C4's93.4687809205s failed-
attempt mean is ineligible, not a fast solution. No holdout IC/rho winner yet.
Two known training failures mean118 is planned, not an actual verified count.

Selection commit19664f02fdcb5f3e5674f5a26e989a117be2f0bc preceded fresh
holdout generation outside1913signedorbits. Remaining37matrix, controller,
identities and inspectedtrainingreport published1e93e4ca1 before cold launch.
CPU3 reserve and first cold IC worker observed active; cold.exit.txt absent.
Private current handoff /workspace/n83_basis_setup/CONTINUE.md. Completionwatch
remains enabled through all remaining stages and deliverables. No heavy local
work while any measured reservation remains active.

Training source freeze e07da4efdc092d6d7d1d7b660b29225a665d0775 and all frozen
predecessors remain unchanged. Dense layout and base selection are distinct.
The draft first600 candidate equalled K600 and was corrected before timing to
the next600disjoint orbits. Its draft dump and all preflight failures remain.
Original retained K600 construction and OS-cache flushing remain unpriced;
no1685targetrun, calibratedS, L2, independent-host, asymptotic or securityclaim.

Validation: unchanged core release library4089pass/0fail/117ignored;
frozenworkers19IC/4rho; trainingcontroller6 and holdoutcontroller6 tests pass;
independent transcript1positive/14negative; native selection tests exclude
fastest censored arm and reject missing audit/overwrite. Conductor unavailable;
sparse/loopback/test-fixture failures retained and corrected before timing.

Cold observed means: C3 611.139215194875s, C0 1223.7821487055s,
rho2.10772421875s; C3/rho289.952171995877 conditional95%CI
[245.967289930960,375.396780301831]. Basiscomparison C3/C0
0.499385626634 [0.484205609743,0.514186089767]. All26coldprocesses/
26canonicalscalars/16fullrank-log audits+oneAAworkpair passed. Native20k paired
bootstrap; no cold IC-over-rho win. [Cold PDF](reports/run-20261015T1842Z/n83_basis_cold_holdout_report.pdf).
Current routing: [cohort continuation](COHORT_CONTINUE.md); historical earlier
running notes above are superseded by this result and the live private handoff.

The final native cohort worker was observed active after launch about23:11UTC,
2026-10-15, using the already-frozen37matrix. Cold-stage data/native statistics/
verification/report published158c334d7 beforelaunch; completionwatch confirmed
enabled through all final deliverables. No heavy local work during reservation.
