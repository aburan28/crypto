# Stage 44: deterministic parallel n=53 pair-pair queries

Stage 42 spends 31.132 of its 38.531 charged direct seconds scanning the signed-expanded pair table. The scan is read-only and divides into cursor intervals, but an ordinary parallel `find` would make both the selected witness and charged overshoot scheduler-dependent.

Stage 44 adds `pair_pair_parallel_512`. A query partitions the table's signed cursor sequence into fixed 512-cursor chunks. Each wave executes exactly one chunk per Rayon thread, charges every addition, inversion, x-filter rejection, and exact miss in every chunk, then consumes results in cursor order and keeps the earliest witness. All chunks in the winning wave remain charged. The relation sequence is therefore deterministic for a fixed input and thread count.

An n=13 control with four threads produced the same 34 relation hashes, trial count, and recovered scalar as sequential `pair_pair_256`. It charged 8,874 support queries versus 8,704 sequential queries, exposing the expected wave overshoot.

On the retired n=53 eta-1/128 tuning seed, the parallel mode produced the same 142 relation hashes and recovered scalar. Collection fell from 15.204 to 4.555 seconds and whole direct wall from 21.520 to 10.658 seconds, a 2.019-fold speedup. Process CPU rose 4.2 percent, from 21.481 to 22.387 core-seconds, while peak RSS stayed near 1.563 GB. Charged support queries rose from 59,250,176 to 59,337,883.

On the already-consumed local same-target stream, direct wall fell from 25.373 to 11.906 seconds and the direct/rho ratio fell from 2.973 to 1.395. The relation transcript, target, and recovered scalar were identical. The hosted A/B run will execute sequential direct, four-thread direct, and rho from one clean build on the exact Stage 42 public target and require complete relation-hash equality.

This is a parallel finite-constant experiment. It does not erase additional core-seconds, change the factor base, provide an unknown-scalar n=53 result, satisfy licensed Magma or external review, change the asymptotic exponent, or establish a Koblitz index-calculus SOTA result.

## Hosted result

[Run `34707644672`](https://github.com/aburan28/crypto/actions/runs/34707644672) completed from merge commit `14ceb064a93b1b3c4c6cf3f4318f656311c98055` and passed a fresh replay after download. Sequential and parallel direct arms emitted the same 189 relation hashes, full-rank point, factor-base log solution, target point, and recovered scalar.

Sequential direct took 48.359 seconds wall and 48.352 core-seconds. Four-thread direct took 30.635 seconds wall and 71.702 core-seconds, a 1.579-fold wall speedup at 1.483 times the process CPU. Peak process RSS stayed essentially unchanged at 1.526 GB. Collection fell from 39.310 to 21.650 seconds, a 1.816-fold speedup. Parallel waves charged 84,721,666 support queries versus 84,606,720 sequential queries.

Packed signed-Frobenius rho on the same public target took 5.128 seconds. The parallel direct/rho wall ratio is therefore a 5.975-times loss, improved from the same-run sequential 9.431-times loss and Stage 42's earlier 12.131-times loss. Fresh build plus parallel direct remains 25.282 times rho wall.

The clean build plus sequential direct, parallel direct, and rho used 183.251 seconds sequential wall, 472.972 core-seconds, and 1,709,301,760 bytes maximum sampled process-tree RSS. Parallelism is retained as the n=53 latency record, not as a core-cost or full-cost crossover.
