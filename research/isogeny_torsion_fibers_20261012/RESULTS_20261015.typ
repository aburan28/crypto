#set page(paper: "a4", margin: 18mm)
#set text(size: 10pt)
= Isogeny and torsion fiber next gates
*Results report, 2026-10-15.* Native Rust, public frozen targets, shared-host L0.

The finite χ polynomial scan checked 22/22 exhaustive cases, solved 12/12 paired relation rounds, and independently replayed 99 scalar records across nine arms. The later direct source S3 and quotient-root controls each matched the independent signed-pair sets on 22/22 smallest-cell targets and solved 12/12 paired rounds. Their 11-arm run independently replayed 123/132 scalar outputs; nine raw-pair rows exhausted their cap. Projected-pair, χ finite-scan, source-root and quotient-root relation counts, rank and trials were identical in every paired round. The direct source-root/negation-rho paired pipeline CPU ratios were 15.474, 6.014, 14.669 and 30.631 on the four frozen cells; quotient-root ratios were 14.943, 8.328, 17.690 and 30.441. These pipeline CPU costs exclude common input reading and calibration; the portable counted S conversion is incomplete.

An explicit degree-three map over `F_(67²)` links curves with distinct j values, 7 and 11. Source, transported-image and destination-direct relation arms reached rank nine and recovered scalar 41 in three rounds; 15/15 arm-round values passed independent replay. An exhaustive audit mapped 4575 source points onto 1525 image points, with exactly three rational preimages each. All 61 target-subgroup images were distinct and 3721 ordered subgroup-pair homomorphism checks passed. Lean proves conditional subgroup injectivity and relation reflection; the concrete curve/map and S3 polynomial bridges remain open.

#image("chi_polynomial_20261015/reports/run-20261015T1828Z/GRAPH.svg", width: 100%)

At 64 distinct order-67 targets over F401, the post-suite V7B run priced cold direct source-group BSGS against direct prime-field rho and pairing with reused H and a finite-field BSGS table. Median complete CPU was 0.378, 0.938 and 1.219 ms, respectively. Paired pairing/source-BSGS ratios were 3.224, 3.696 and 3.674; their median was 3.674. Independent arithmetic replayed all 594 scalar outputs. Source BSGS was lower at every prefix in all three rounds. V7, measured with release tests running, had a median ratio of 5.190 at 64 targets; the difference is shared-host timing variation.

#image("pairing_batch_20261015/reports/run-20261015T1822Z/GRAPH.svg", width: 100%)

These results leave the direct χ-coordinate resultant, a full-image destination factor base, portable counted primitive conversion, isolated increasing-order comparisons, independent-host replay, and the concrete Lean curve bridges open. The study run reports contain frozen specs, hashes, raw receipts, failed attempts and exact replay. The canonical ECDLP frontier is unchanged.
