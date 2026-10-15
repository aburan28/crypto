#set page(paper: "a4", margin: 18mm)
#set text(size: 10pt)
= Isogeny and torsion fiber next gates
*Results report, 2026-10-15.* Native Rust, public frozen targets, shared-host L0.

The nonfactorable χ polynomial decomposition checked 22/22 exhaustive cases, completed 12/12 paired relation solves, and independently replayed 99 scalar records over its nine-arm run. Projected-pair and χ arms had the same relation counts, rank and query trials in every paired round. The χ CPU-equivalent ratios to negation rho were 29.64, 26.29, 89.42 and 142.80 at the four frozen cells. The formal counted S conversion for the polynomial primitive is incomplete.

An explicit degree-three map over `F_(67²)` links curves of different j, 7 and 11. Source, transported-image and destination-direct relation arms attained rank nine and recovered scalar 41 in three rounds; 15/15 arm-round scalars passed independent arithmetic replay on source and destination. An exhaustive audit mapped 4575 source points onto 1525 image points, with exactly three rational preimages each. All 61 target-subgroup images were distinct and 3721 ordered subgroup-pair homomorphism checks passed. Lean checks the conditional order-61 kernel-to-injectivity implication; the concrete curve and map bridge remains open.

At 64 distinct order-67 targets over F401, V7 priced cold direct source-group BSGS as the missing multi-target comparator. Median complete CPU was 0.303 ms for its nine-entry table and all 64 solves, 0.577 ms for direct prime-field rho, and 1.740 ms for pairing with reused H and finite-field table. Paired pairing/source-BSGS ratios were 5.190, 4.178 and 5.749. Source BSGS was lower cost at every prefix in each round. V7 ran alongside the release tests. A separate V6B repeat after those tests exited still found direct rho below pairing in all three rounds, with paired ratios 1.944, 2.470 and 1.992.

#image("pairing_batch_20261015/reports/run-20261015T1805Z/GRAPH.svg", width: 100%)

All measurements are shared-host L0. Direct resultant/root and S3/quotient polynomial controls, portable counted-S conversion, isolated scaling, independent-host replay, and concrete curve-level Lean proof remain open. The study run reports link frozen specs, hashes, raw records, failures and independent replay. The canonical ECDLP frontier is unchanged. V6, V6B and V7 are retained as separate runs.
