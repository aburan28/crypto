# Frozen V3 balanced-order remediation

Date: 2026-10-13 local. Commit before changing or running V3. V1 and V2
records remain unchanged.

V2 passed 640/640 measured solver rows and the independent arithmetic check.
Its thread-CPU A/A ratios were `0.826851`, `0.923947`, `0.754165`,
`0.817084` at subgroup orders `103`, `439`, `1613`, `6427`. The order-1613
row misses V2's frozen `[0.8,1.25]` noise gate. Under V2 rotation, rho
usually preceded its duplicate; cache state is a plausible cause, not a
demonstrated explanation.

V3 changes only arm ordering. For `(round+target_index)` even, execute
`[source, transport, rho, rho_AA]`; for odd, execute
`[rho_AA, rho, transport, source]`. Across eight targets and each round,
each duplicate rho copy executes first equally often. Keep the exact four
pairs, eight targets, five measured rounds, warmup, seeds, solver algorithms,
caps, counters, thread-CPU clock, wall/phase clocks, and checker unchanged.
The V3 output is `raw-v3.tsv` with its own summary and receipts.

Require every measured row verified, matching A/A counts, independent
arithmetic checks, and an A/A thread-CPU total ratio inside `[0.8,1.25]`
at every size. Stop and preserve V3 if it fails. Local thread-CPU costs
remain L0 and do not substitute for L2 wall isolation or second-host replay.
