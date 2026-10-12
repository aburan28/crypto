# Stage 167: deterministic parallel fixed-X1 F4 batches

Strong target-specific parallel engineering on one finite public toy PDP. Batch width 13 was selected post hoc because this target's first valid relation occurs in rational system 39. This is not fresh-target evidence, a full solver panel, a full index-calculus/rho crossover, independent reproduction, novelty review, or SOTA.

The algebraic fixed-X1 systems are independent. The backend now gathers a deterministic batch of rational fixed-X1 systems, solves the batch through the bounded Rayon pool, waits for every launched system, charges every result, and then inspects results in the original schedule order. `PQ_F4_X1_BATCH=1` and one Rayon thread retain the prior single-thread path.

The clean 13-thread F4 process takes 5.589951 wall seconds, 36.302718 core-seconds, and 3033563136 bytes peak RSS. The same clean binary at batch one takes 19.828459 wall / 19.805719 core-seconds / 904085504 bytes. Parallel wall improves 3.547x while CPU rises 1.833x and RSS rises 3.355x.

Batch 13 was selected post hoc on this target: the known successful system is the 39th rational system, so three batches cover it without speculative systems. The execution mechanism does not use the witness, truth labels, subgroup enumeration, known scalar, or factor-base logarithms, but this batch-width result is target-specific and is not fresh-target evidence.

Direct MITM remains 1.86x faster by wall, 12.13x cheaper by CPU, and 66.17x smaller by RSS. The parallel arm does not pass the full-cost gate.

Licensed Magma F4, the complete same-instance solver panel, fresh-target validation, full index-calculus/rho crossover, and unaffiliated reproduction and novelty review remain open. This is not a Koblitz index-calculus SOTA result.
