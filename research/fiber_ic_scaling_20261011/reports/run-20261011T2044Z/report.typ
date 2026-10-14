#set page(paper: "a4", margin: 1.8cm)
#set text(size: 10pt)
= Run report: scaling_h72_small_b
_2026-10-11T20:44:14Z to 2026-10-11T20:44:48Z UTC. research/fiber_ic_scaling_20261011._

*Spec* specs/scaling_h72_small_b.json. *Session* ECBS1h54835df66b11 at commit 5c8a7694f (git_dirty false), run with `ecbench run --wait` and audited with `ecbench verify --replay 12`: ok, 144 of 144 verified, 12 replays reproduced, 0 problems. *Host:* cloud container, Xeon @ 2.80 GHz; isolation earned L2 on 88, L1 on 56 (spec asked for L2).

*Measured:* each arm's S (ecbench.gae cost over the floor) per curve, 6 runs each (3 targets × 2 rounds). Per-curve means: TABLES.md, section scaling_h72_small_b.

#figure(image("s_by_arm.svg", width: 95%), caption: [Session ECBS1h54835df66b11: mean S per arm against log#sub[2] \#E (bars min–max over 6 runs). Records sessions/scaling_h72_small_b/records.jsonl.])
