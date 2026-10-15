#set page(paper: "a4", margin: 18mm)
#set text(size: 9.5pt)
= Direct Semaev and quotient polynomial relation solve
*Run:* 2026-10-15 18:28 UTC. *Host:* Darwin arm64, shared L0. *Protocol:* frozen `RUN_SPEC_V2.json`.

The source and quotient arms each passed 22/22 exhaustive smallest-cell target sets, then returned the scalar in all 12 paired relation rounds over four public `F_(p²)` cells. Separate `u128` arithmetic replayed 123/132 arm-round scalars; the other nine were raw-pair control cap exhaustions. Both new arms matched the projected-pair and χ finite-scan arms in relation count, useful rank and query trials in every paired round.

For V2 the runner read each fixture's planted answer for a final framework assertion; neither root oracle nor collector received it. A separately frozen V3 repeat removes the answer from the entire solver input and leaves it only with the checker.

The single measured cost unit is S_cpu = per-arm pipeline current-thread CPU / calibrated source addition CPU / sqrt(r). It charges base construction, polynomial/map work, relation collection, rank/linear solve, target recovery and verification. Common fixture reading, calibration measurement and process startup are outside each arm. It is a shared-host L0 diagnostic, not portable counted S.

At r=19/61/211/811, median S_cpu was 16.8/23.5/18.4/6.8 for direct negation rho; 281.7/132.2/141.1/200.9 for projected pairs; 267.5/131.2/174.6/196.1 for direct source S3 roots; and 265.2/242.4/174.5/202.6 for quotient S3 roots. Paired source-root/rho ratios were 15.474/6.014/14.669/30.631. χ finite-scan values were 666.2/263.2/626.5/2067.0, and paired finite-scan/source-root ratios were 2.124/2.030/3.848/7.234. This is a local polynomial-implementation gain at unchanged relation support, while rho remained lower. The raw-pair baseline exhausted its cap in all nine larger-cell slots.

#image("GRAPH.svg", width: 100%)

On p=3319, paired source-root/projected-pair ratios were 1.422, 0.640 and 0.752, with equal relation outcomes; an L0 rho A/A ratio as low as 0.522 on p=271 shows the timing spread. Neither the quotient nor direct root construction produced a relation-yield gain on this fixed base. Portable primitive pricing, isolation, larger orders and independent-host replay remain open. The canonical frontier is unchanged. Raw SHA-256: `4f0c090af106be03d39dac2379bc0c4c968e095a0e7b472a814dcbca0569c0fc`; replay SHA-256: `bbe6ad66572b0741a646c051662f4501319f1dfc61f4fd100717d80adc3af608`.
