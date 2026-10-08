# Stage 169: preregistered fresh n=59 target holdout

One preregistered same-cell fresh-target holdout for the fixed schedule and batch policy. It cannot establish a full solver panel, expected-time law, end-to-end index-calculus/rho crossover, independent reproduction, novelty, or SOTA.

The holdout was selected from the blind bundle by presentation order and committed before execution. Batch size 39, twelve Rayon threads, a 300-second F4 budget, and no retry after truth were frozen in the preregistration.

Native F4 exhaustively visits all 512 masks, skips 270 non-rational values, and completes 242 F4 systems. It returns UNSAT in 20.539819 wall seconds / 169.094100 core-seconds / 2872541184 bytes peak RSS. Post-run scoring classifies it as a true negative.

Same-host direct MITM also returns UNSAT in 2.953129 wall seconds. Native XOR is inconclusive at 100,000 conflicts, while WDSat and CryptoMiniSat reach their 120-second watchdogs. Licensed Magma remains unexecuted.

This holdout proves that the fixed policy can complete a fresh nondecomposable target and produce a correct exhaustive UNSAT terminal. Its 20.54-second wall cost is materially above the post-hoc 4.99-second positive target and direct MITM, so it rejects any inference that Stage 168 represents generic fresh-target time.

One fresh holdout is finite evidence, not an expected-time distribution or full panel. Licensed Magma, the full native-F4 panel, end-to-end IC/rho crossover, and unaffiliated reproduction and novelty review remain open. This is not a Koblitz index-calculus SOTA result.
