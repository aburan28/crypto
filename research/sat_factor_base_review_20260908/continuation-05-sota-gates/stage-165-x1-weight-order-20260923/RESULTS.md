# Stage 165: target-independent sparse fixed-X1 order

Strong target-specific internal engineering on one finite public toy PDP. The order is target-independent, but it was selected post hoc after inspecting this target. This is not generic expected-time evidence, a full solver panel, a full index-calculus/rho crossover, independent reproduction, novelty review, or SOTA.

The fixed-X1 outer loop now orders algebraic-subspace coefficient masks by Hamming weight and then numeric mask. The order depends only on the published polynomial-basis coefficients, not on the target, truth class, witness, subgroup enumeration, factor-base logs, or a known scalar. `PQ_F4_X1_ORDER=ascending` retains the prior same-binary schedule.

The clean blind F4 process takes 21.512089 wall seconds, 21.484441 core-seconds, and 948502528 bytes peak RSS. It visits 85 masks, constructs 39 rational fixed-X1 systems, and charges 16,821,055,616 elimination-and-table word XORs before returning the same exact curve-verified witness.

The clean ascending control takes 34.446599 wall seconds and 34.412474 core-seconds. The selected clean wall speedup is 1.601x; the interleaved three-run development medians give 1.602x.

Direct MITM still wins by 7.17x wall and 7.18x CPU. The order was selected post hoc after inspecting this target, so the result is a target-specific engineering optimization despite the order itself being target-independent; it is not evidence of generic expected-time scaling.

Licensed Magma F4, the complete same-instance solver panel, full index-calculus/rho crossover, and unaffiliated reproduction and novelty review remain open. This is not a Koblitz index-calculus SOTA result.
