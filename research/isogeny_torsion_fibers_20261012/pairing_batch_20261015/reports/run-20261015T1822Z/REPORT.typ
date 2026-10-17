#set page(paper: "a4", margin: 19mm)
#set text(size: 10pt)
= V7B post-suite pairing, source BSGS and direct rho
*Run:* 2026-10-15 18:22 UTC. *Host:* Darwin arm64, shared L0. *Protocol:* `REPEAT_V7B.md`.

The controlled release library suite exited successfully before V7B: 4,083 passed, zero failed, 117 ignored. The frozen source and binary SHA-256 hashes matched V7. All three arms returned all 66 scalars in each of three rotated rounds, and separate source arithmetic replayed 594/594 outputs. The same public inputs, target order, cap, H selection, and nine-entry source BSGS table policy were used. Each arm's cold CPU includes the common input cost, setup, every solve and final checks.

At N=64, median complete current-thread CPU was 0.378 ms for direct source BSGS, 0.938 ms for direct prime-field rho, and 1.219 ms for pairing. Paired pairing/BSGS ratios were 3.224, 3.696, 3.674; their median was 3.674. Direct BSGS was the lower whole-CPU arm at every prefix in every round. The V7 concurrent-test median ratio was 5.190; the timing variation prevents an isolated effect attribution.

#image("GRAPH.svg", width: 100%)

The BSGS full-list ledger has 9 baby entries, 283 lookups, zero false candidates, 488 additions, 552 doublings, 132 scalar multiplications, 3,726 base-field multiplications, 1,586 squarings and 1,032 inversions. The raw records retain the pairing/rho primitive ledgers, attempts, wall time and RSS. These L0 observations do not establish portable counted S or a growing-size crossover. Raw SHA-256: `9208af97ef9e2c0db554c887e324249a2c69f8991f02543d4c7f2e722b26bf70`; replay SHA-256: `68aa6cc8e9c3ee48ce1d7135da7326c08ac2345dd94d647124606d3c6376994d`.
