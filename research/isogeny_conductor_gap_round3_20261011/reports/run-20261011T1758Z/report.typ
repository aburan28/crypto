#set page(paper: "a4", margin: 1.8cm)
#set text(size: 10pt)
= Run report: E5 session e5_combined_n17_t1
_2026-10-11T17:58:53Z to 2026-10-11T18:14:18Z UTC. research/isogeny_conductor_gap_round3_20261011, lead E5._

*Spec* specs/e5_combined_n17_t1.json; *session* ECBS1he4c2c379c257 at commit 2fa5132e8 (git_dirty false); run with `ecbench run --wait`, audited with `ecbench verify --replay 12`: ok, 8 of 8 verified, 8 replays reproduced, 0 problems. *Host:* cloud container, Xeon @ 2.80 GHz; isolation earned L1 on 6 records and L2 on 2 (the spec asked for L2).

*Workload:* icv1-f2m17-tm101-00378d4e (a = 1 Koblitz curve over F#sub[2#super[17]], h = 2), target index 0. *Arms:* rho.negation; ic.pipeline on binary-subspace:dimension=8 with the blind F4 descent, the all-lifts descent and the combined descent (f4-f2). 8 runs requested, 8 verified; no warm-up, 2 measured rounds.

*Derived:* combined ÷ blind F4 solver operations, paired by round: 234× (round 0) and 411× (round 1). Solver operations are the solver's own count; the session's S leaves them unpriced (lower bound).

#figure(image("e5_arms.svg", width: 95%), caption: [Session ECBS1he4c2c379c257: online wall per target and F4 solver operations per run, 2 measured rounds per arm (log scales). Records sessions/e5_combined_n17_t1/records.jsonl.])
