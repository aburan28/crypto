#set page(paper: "a4", margin: 19mm)
#set text(size: 10pt)
= V7 reusable source BSGS versus pairing and direct rho
*Run:* 2026-10-15 18:05 UTC. *Host:* Darwin arm64, shared L0, with concurrent release tests. *Protocol:* `RUN_SPEC_V7.json`.

On 66 distinct public order-67 targets, cold pairing with reused H and finite-field BSGS was paired with direct F401 rho and a direct F401 BSGS arm that rebuilt its nine-entry baby table in each round. All three arms returned all 66 scalars in three rotated rounds; separate source arithmetic replayed 594/594 outputs. Public input, setup, failed work, per-target solves and verification are included.

At N=64, median complete CPU was 0.303 ms for direct BSGS, 0.577 ms for direct rho, and 1.740 ms for pairing. Paired pairing/direct-BSGS ratios were 5.190, 4.178, 5.749. The median ratio rose from 2.121 at N=1 to 5.225 at N=66. Source BSGS was the cheapest complete arm at each prefix in this shared-host run.

#image("GRAPH.svg", width: 100%)

The direct BSGS arm charged 9 baby entries and 283 lookups over the full list. Its final primitive ledger includes 488 additions, 552 doublings, 132 scalar multiplications, 3,726 base-field multiplications, 1,586 squarings and 1,032 inversions. The raw record retains all three arms' ledgers, wall times, attempts and RSS. Shared-host timing does not establish isolated wall speed, portable counted S, or a growing-size crossover. Raw SHA-256: `7da080bc0e3a1702b58ae9659564fed57f58c65f5662b628ae6d683270271694`; replay SHA-256: `68aa6cc8e9c3ee48ce1d7135da7326c08ac2345dd94d647124606d3c6376994d`.
