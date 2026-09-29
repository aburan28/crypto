# G2 diagnostic, after the registered runs (not a registered measurement)

G2 failed as registered on the frozen `K_0/2^13 m=2` rung in all three arms
(relative differences 2.5e-5, 7.0e-6 and 8.6e-6 against a gate of 1e-6).
To find the cause, the A2 arm of that rung was run twice more, single-threaded:

    RAYON_NUM_THREADS=1 KIC_CHAIN_ORDER=interleaved KIC_LINEAR_ELIM=1 KIC_F4_MULTIPLIERS=support KIC_F4_DROP=complete \
      valgrind --tool=callgrind --toggle-collect='*groebner_decompose*' \
      target/release/examples/groebner_stage_bench --rung 2 --out ...

Both runs: 693,232,587 instructions, identical.  With the default thread pool
the two registered runs gave 693,209,626 and 693,215,570.  The variation is
rayon's scheduling in the parallel Macaulay elimination, which that rung's
root matrices are large enough to take; the single-threaded count is 3.4e-5
above the threaded one.  Every counter is identical in all four runs.
