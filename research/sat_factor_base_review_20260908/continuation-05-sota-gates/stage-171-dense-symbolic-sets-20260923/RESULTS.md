# Stage 171: dense symbolic monomial sets in native F4

Same-target engineering after the Stage-169 truth was opened. This stage improves the repository-native Boolean F4 arm; it is not a fresh-target distribution, complete index calculus, rho crossover, independent reproduction, novelty review, or SOTA.

For systems with at most twenty Boolean variables, F4 now represents the symbolic-preprocessing LCM, examined, and no-divisor sets with one bit per complete-domain monomial plus an insertion-order vector. The allocations are reused across steps. `PQ_F4_DISABLE_DENSE_SYMBOLIC_SET=1` restores the same binary's hash-set path. The fixed-X1 batch cap now reaches the complete public `2^ell` domain, while the default remains one and `PQ_F4_X1_BATCH=1` is the serial scheduling control.

The first clean selected run uses twelve Rayon threads and one complete batch of 242 rational systems. It finishes the known Stage-169 true-negative target in 14.992139 wall seconds / 159.504427 total core-seconds / 2932424704 bytes peak RSS. It visits all 512 masks, preserves 14,278 equations and 14,515,915 terms, performs 242 F4 calls and 99,199,976,264 charged word XORs, finds zero roots, and returns exhaustive UNSAT with equation fingerprint `02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb`.

Across three interleaved clean pairs, the median candidate/control ratios are 0.986961 wall, 0.987166 CPU, 0.914624 RSS, and 0.869077 F4 matrix-build time. The scientific reports agree after removing the declared representation, allocation, and timing counters. The selected run is 1.143x faster in wall time than Stage 170's selected run, but it does not beat Stage 170's one-off 14.661977-second post-hoc wall sample.

The exact-commit clean archive/build costs 170.147153 wall seconds, 166.665865 core-seconds, and 1714372608 bytes peak RSS. Build plus selected execution costs 185.139292 wall seconds and 326.170292 core-seconds. Two failed offline setup attempts are preserved and charged.

Direct MITM remains 5.08x faster by wall, 54.10x cheaper by CPU, and 65.13x smaller by RSS on this target. The broader campaign has a 49-component Stage-171 measured lower bound, while complete development cost remains null because some exploratory compilation/test work was not outer-metered and one early profile receipt was overwritten.
