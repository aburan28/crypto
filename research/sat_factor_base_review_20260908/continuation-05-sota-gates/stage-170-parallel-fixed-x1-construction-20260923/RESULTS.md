# Stage 170: parallel fixed-X1 S4 construction

Same-target engineering after the Stage-169 truth was opened. It establishes an exact deterministic construction-stage speedup for the repository-native F4 arm, not a fresh-target distribution, full index calculus, independent reproduction, novelty, or SOTA.

Each rational fixed-X1 S4 system in a deterministic batch is now constructed through the same bounded Rayon pool used by the repository-native F4 solver. Indexed parallel collection preserves algebraic mask order, and `PQ_F4_DISABLE_PARALLEL_CONSTRUCTION=1` restores serial construction in the same binary.

The clean selected batch-39/twelve-thread process completes the known Stage-169 true-negative target in 17.131544 wall seconds / 176.171665 total core-seconds / 2904489984 bytes peak RSS. The same binary with serial construction takes 20.796178 seconds. This is a 1.214x wall speedup; construction falls from 4.629543 to 0.528026 seconds.

Three alternating development pairs have a 1.206x median wall speedup. Every arm preserves the 14,278 equations, 14,515,915 terms, equation fingerprint `02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb`, 242 F4 calls, 99,199,976,264 word XORs, zero roots, and exhaustive UNSAT classification.

An explicitly post-hoc 14-thread/batch-64 clean arm reaches 14.661977 seconds but raises peak RSS and can do extra speculative work on decomposable targets. It is retained as a target-tuned ceiling rather than the selected generic policy.

Direct MITM remains 5.80x faster by wall and 59.75x cheaper by CPU on this target. The build plus selected run costs 186.059675 wall seconds, so the full-cost gate remains false.

This optimization was developed after the Stage-169 truth was opened. It strengthens implementation performance but is not a second fresh holdout, a full target distribution, a complete index-calculus/rho crossover, independent reproduction, novelty review, or SOTA.
