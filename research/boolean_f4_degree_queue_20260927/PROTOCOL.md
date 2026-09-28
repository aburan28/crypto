# Boolean F4 pair queues grouped by selection degree

## Candidate and reference

The measured batched insertion default from PR #882 still scans and partitions the entire pair vector to select the next minimum degree. The opt-in `F4_F2_DEGREE_BUCKETS=1` candidate appends each pair to its degree's vector, retains insertion order within each degree, and takes the first nonempty bucket without rescanning all pairs. The same-binary reference is `=0`, the original flat vector. Set `F4_F2_BATCH_INSERTS=1` and `F4_F2_BITMAP_SEEN=1` in both arms. A deterministic unit test compares the selected pair sequence under random insertions and filtering.

## Frozen paired measurement

Use `examples/f4_f2_bench.rs`, one repetition of all seven standard Boolean F4 cases per process, primary `n20_m30`. Pair the frozen seed XOR zero and holdouts `1ac0ffee` and `2468ace0`. Pin one Linux CPU and set `RAYON_NUM_THREADS=1`. Warm each arm once per workload, then run five A/A reference pairs and five alternating-order A/B pairs. Retain every process outcome, full output, source and binary hashes, CPU, load, phase timings and counters in the raw receipt. Compute exact five-pair bootstrap 95% intervals. The target phase is `other_ms = wall_ms - build_ms - eliminate_ms`; complete `wall_ms` excludes process launch and fixture construction.

Exact basis fingerprints, matrix dimensions, pair and skip counters, divisor tests and logical word XORs must match. A stage gain requires the frozen primary other-time ratio to reach 1.05 with lower interval bound above one, its complete-call ratio and lower bound above one, both holdout primary gains outside their A/A ranges, and no smaller complete-call median below its own A/A range. If the one-thread gates pass, run the same paired four-thread control before a default change. Preserve negative results and leave the flat default if a gate fails. This is an internal Boolean F4 solver-stage study, not a one-target IC or DLP speedup claim.
