# Stage 172: native F4 single-core completion

This stage fills the missing valid single-core resource row for the repository-native Boolean F4 implementation on the already-opened Stage-169 `n=59, ell=9, m=3` true-negative target. It is not a fresh holdout, complete solver panel, full index-calculus/rho crossover, independent reproduction, novelty review, or SOTA.

With `RAYON_NUM_THREADS=1`, `PQ_F4_X1_BATCH=1`, and every BLAS-style thread variable fixed to one, dense symbolic sets complete exhaustive UNSAT in **107.925193 wall seconds / 107.530709 core-seconds / 417529856 bytes peak RSS**. The backend reports `single_thread_requested=true`, 242 sequential batches, all 242 systems complete, 14,278 equations, 14,515,915 terms, 99,199,976,264 word XORs, zero roots, and equation fingerprint `02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb`.

The same clean binary with `PQ_F4_DISABLE_DENSE_SYMBOLIC_SET=1` takes 113.184864 wall seconds / 112.835489 core-seconds / 398966784 bytes. Dense sets reduce wall by 4.65%, CPU by 4.70%, and F4 build time by 11.69%; RSS rises 4.65% in this one-worker shape.

The same-binary direct-MITM reference completes in 3.136514 wall seconds / 3.118220 core-seconds / 42057728 bytes, leaving single-core F4 34.41x slower by wall and 34.48x more expensive by CPU. The Stage-171 parallel F4 run is 7.20x faster by wall but consumes 1.48x the CPU and 7.02x the memory.

The inherited exact-commit clean build plus this single-core execution costs 278.072346 wall seconds and 274.196574 core-seconds. All three Stage-172 execution arms are outer-metered; the historical campaign total remains null for the inherited unmetered scope.
