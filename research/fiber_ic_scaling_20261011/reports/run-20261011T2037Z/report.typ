#set page(paper: "a4", margin: 1.8cm)
#set text(size: 10pt)
= Run report: scaling_h4
_2026-10-11T20:37:11Z to 2026-10-11T20:39:14Z UTC. research/fiber_ic_scaling_20261011._

*Spec* specs/scaling_h4.json. *Session* ECBS1h4a5009a8f716 at commit 0314b6a7a (git_dirty false), run with `ecbench run --wait` and audited with `ecbench verify --replay 12`: ok, 336 of 336 verified, 12 replays reproduced, 0 problems. *Host:* cloud container, Xeon @ 2.80 GHz; isolation earned L2 on 231, L1 on 105 (spec asked for L2).

*Measured:* each arm's S (ecbench.gae cost over the floor) per curve, 6 runs each (3 targets × 2 rounds). Per-curve means: TABLES.md, section scaling_h4.

#figure(image("s_by_arm.svg", width: 95%), caption: [Session ECBS1h4a5009a8f716: mean S per arm against log#sub[2] \#E (bars min–max over 6 runs). Records sessions/scaling_h4/records.jsonl.])
