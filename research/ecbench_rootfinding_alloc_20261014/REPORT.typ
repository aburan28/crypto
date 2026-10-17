#set page(margin: 18mm)
#set text(size: 10pt)

= Reusing prime lane collector buffers

The prime relation collector now reuses four vectors between synchronized lane rounds: requests, native left and right points, and outputs. The release stream test passed at widths 4, 8, and 16. The saved baseline and candidate expanded to identical public-target benchmark plans.

#image("timing-ratios.svg", width: 100%)

The chart shows L0 wall-time ratios on one shared macOS arm64 host. Each point is a geometric mean over matched measured rounds; whiskers are target-cluster bootstrap intervals over eight public targets. Blue compares the candidate with the first baseline. Orange compares the repeated baseline with the first baseline. The dashed line is parity.

== Results and scope

The three sealed sessions contain 384/384 verified executions and 288/288 identical measured replays. The corrected cross-session checker found zero mismatches among compared deterministic fields in 128 candidate/baseline record pairs and 128 baseline/baseline pairs. All IC arms had mean cold counted cost 7,952.58 GAE; sealed native field work remains unpriced. The native ecbench comparison reports counted-cost ratio 1.000 over 24 measured batch-8 pairs and refuses a cross-session wall comparison. In the candidate session, batch-8 IC / matched rho has lower-bound counted-cost ratio 6.8747 (interval 5.6689--8.4198).

Custom L0 diagnostics: batch-8 relation-phase wall ratio 0.937 (95% target-cluster interval 0.786--1.122) and whole-solve wall ratio 0.993 (0.803--1.251). The repeated baseline ratios were 1.057 (0.814--1.362) and 1.063 (0.818--1.390), respectively. Every interval includes parity. The run verifies unchanged deterministic behavior and the source-level buffer reuse, while the L0 timings do not establish a wall-time improvement.

The baseline binary SHA-256 is `ab5764921ad616ae8177ca57c125ec8cd7ba09896f1e0644a6999d1fb33b7159`, bound to source commit `431d0b3167027092443705e21d363489b2637a22` by the prior postmerge build receipt. The candidate binary SHA-256 is `db9cf52602e8b3bc19a885c3eae99c899bd6acc541958246f3e540e58e8de743`, with session source commit `316cc8386611d9b6ba7a57c86442b4d59ba9ee58`.

The first release build failed on temporary-storage exhaustion and the first comparison checker classified phase wall time as a deterministic mismatch. Both failed outputs and their corrections are retained. Raw sessions, replay audits, comparison JSON, protocol, and reproduction commands are linked from `REPORT.md`.
