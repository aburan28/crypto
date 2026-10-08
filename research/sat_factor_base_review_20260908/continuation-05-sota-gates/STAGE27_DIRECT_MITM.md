# Stage 27: matched direct MITM matrix

Stage 27 runs the existing exhaustive direct meet-in-the-middle decomposition backend against the exact 160 truth-free Stage 26 inputs. The four 40-instance cells execute concurrently, while every cell restricts its complete packet-verification and decomposition process tree to one Linux CPU.

The backend materializes the algebraic factor points, indexes every unordered pair sum, and searches for the third point needed to reach the public target. A SAT result is accepted only when the retained three-point witness sums to that target. An exhaustive miss is UNSAT. The result reports factor-point count, pair-table entries, group additions, elapsed time, total CPU, and memory. SAT conflicts do not apply to this algorithm.

The workflow reuses the exact backend artifact built and charged in Stage 26 run 34632018379. It authenticates binary SHA-256 `667fc00681e114acaaf662a8c00c214354e1a08e236a6974fb08e52504711ea0` and records the originating artifact identity in every cell plan.

This matched direct decomposition lane remains a finite public PDP benchmark. It is not an end-to-end index-calculus or SOTA result.
