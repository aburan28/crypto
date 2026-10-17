#set page(paper: "a4", margin: 19mm)
#set text(size: 10pt)
= V6B pairing batch repeat after release test exit
*Run:* 2026-10-15 17:52 UTC. *Host:* shared Darwin arm64 L0. *Protocol:* `REPEAT_V6B.md`, inheriting frozen V6 inputs, source, caps and arm order.

An initial invocation used the wrong input-cell index and failed before timing; its empty output and stderr are preserved. The corrected three rounds returned all 66 public order-67 scalars in pairing, extension-rho and direct-prime-rho arms. Independent source arithmetic replayed 594/594 scalars with no mismatch.

At 64 distinct targets, median complete cold current-thread CPU was 0.492 ms for direct prime-field rho, 1.315 ms for extension-field rho, and 1.157 ms for pairing with reused H and finite-field table. Paired pairing/direct-rho ratios were 1.944, 2.470 and 1.992. Direct rho was faster in each round; pairing was faster than extension rho in each round.

At full N=66 in round zero, pairing used 9,727 extension multiplications including 5,265 squares and 152 extension inversions. Direct rho charged 67 attempts, 2,334 additions, 1,932 doublings, 14,251 base-field multiplications, 6,125 squarings and 4,124 inversions. Peak process RSS was 7,045,120 bytes. These different primitives lack one counted-S conversion.

#image("GRAPH.svg", width: 100%)

V6 was concurrent with the release tests; V6B started after they exited. Both are shared-host L0 and the difference between versions has no causal attribution. There is no portable counted-S or scaling claim. Raw V6B SHA-256: `d047dc2444f28955023cbb90f0dfbff3afd66f6a28a36e25d0797d7f8b05d3d0`; independent replay SHA-256: `af26814d85377ef24487818b61571cd725972f4563bbf312f4fee825dc0775e9`. The run report links the full receipts and limits.
