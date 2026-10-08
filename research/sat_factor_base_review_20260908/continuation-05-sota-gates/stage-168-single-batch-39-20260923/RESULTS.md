# Stage 168: one deterministic batch of 39 fixed-X1 systems

Strong target-specific parallel scheduling on one finite public toy PDP. Batch size 39 is exactly the known successful prefix and was selected post hoc. This is not a generic stopping rule, fresh-target evidence, a full solver panel, a full index-calculus/rho crossover, independent reproduction, novelty review, or SOTA.

Stage 168 removes the two barriers in Stage 167 by placing the same 39 rational fixed-X1 systems into one deterministic batch. Twelve Rayon workers dynamically execute that complete batch; every system is charged before results are inspected in algebraic schedule order.

The clean F4 process takes 4.987010 wall seconds, 28.179809 total core-seconds, and 2763997184 bytes peak RSS. `single_core_seconds` is null. It completes the same 39 systems, charges 16,821,055,616 word XORs, and returns the same exact witness as the Stage-167 parallel and single-thread arms.

Relative to Stage 167 batch 13, clean wall improves 1.121x, CPU falls to 0.776x, and RSS falls to 0.911x. Three 12-thread development runs have a 4.744985-second median.

Batch size 39 is selected post hoc because it is exactly the successful prefix on this target. This is target-specific scheduling evidence and not a valid generic stopping rule or fresh-target result.

Direct MITM remains 1.66x faster by wall, 9.42x cheaper by CPU, and 60.29x smaller by RSS. The full-cost gate remains false.

Licensed Magma F4, the complete same-instance solver panel, fresh-target validation, full index-calculus/rho crossover, and unaffiliated reproduction and novelty review remain open. This is not a Koblitz index-calculus SOTA result.
