#set page(paper: "a4", margin: 18mm)
#set text(size: 9.5pt)
= Public-only direct Semaev and quotient relation solve
*Run:* 2026-10-15 18:51 UTC. *Host:* Darwin arm64, shared L0. *Protocol:* frozen V3.

The solver used a redacted public fixture; it rejected the answer-bearing file with exit 101. The separate checker proved that the four public objects differed from the private objects only by removed answer fields. Source and quotient S3 root arms each passed 22/22 exhaustive target sets and recovered the scalar in 12/12 paired relation rounds. Separate `u128` arithmetic replayed 123/132 arm-round scalar outputs; nine inherited raw-pair rows exhausted their cap. All 132 semantic rows matched V2 in recovered value, relations, trials, rank and status. V2/V3 timing is not interpreted as a before/after speed comparison.

S_cpu is per-arm pipeline current-thread CPU / calibrated source addition CPU / sqrt(r), with factor-base, polynomial/map, relation, rank/linear solve, target recovery and verification phases inside. Common input parsing, calibration measurement and process startup are outside this timer. At r=19/61/211/811, source-root/rho paired ratios were 13.237/8.254/13.427/21.844; χ finite-scan/source-root ratios were 2.125/1.988/3.864/7.571. Direct rooting improved this bounded polynomial implementation at unchanged relation outcomes; matched rho remained lower.

#image("GRAPH.svg", width: 100%)

This is an answer-separated correctness and shared-host local-cost result. Portable counted S, isolated wall speed, larger subgroup orders and independent-host performance replay remain open. Raw SHA-256: `8a63b7c543392d3dcad9225cf33c0fe1e4adb36b6235dcc8747a697214b1fd55`; replay SHA-256: `024674afd5d78674f807b4a7f20eeed796ed34181fe5e7c2786b08d4b4ef95da`.
