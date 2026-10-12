# Two-block word batching in the shared GF(2) Macaulay kernel

## Frozen question

The matrix-F5 `n24_m24_d4` elimination repeatedly reads its complete matrix after each pivot block. This candidate keeps the next block's pivot-word strip current, updates each selected pivot row on demand, and clears up to two blocks from the same 64-column word while each destination row is resident. It preserves the same RREF output and keeps the prior per-block path as `KIC_GF2_WORD_BATCH=0`. The candidate is `=1`; `KIC_GF2_DEFER_ABOVE=0` in both arms isolates forward batching from the negative deferred-reverse experiment in PR #866.

## Workloads and accounting

Use the unchanged `examples/f4_f2_bench.rs` F5 suite with `n24_m24_d4` as the primary cell and its six smaller cases as controls. Compare the frozen seed XOR zero and holdouts `badc0de1` and `5eed2026`. Each process runs one repetition of all seven cells. Pin the process to one allowed CPU and set `RAYON_NUM_THREADS=1`. Warm each arm once per workload, then run five A/A reference pairs and five alternating-order A/B pairs. The same release binary is used for both modes. Preserve every call, status, output, load and source/binary hashes in a raw receipt. Pair only within one runner and workload. Report exact five-pair bootstrap intervals for reference/candidate ratios.

The primary metric is `reduce_ms`, timed inside the matrix-F5 call. Report complete `wall_ms` beside it. Neither includes process launch or planted fixture construction. The counted word XORs may differ; exact `rows_fp`, rank, pruning, criterion work, and matrix-kernel unit-test RREF must match. This is an internal solver-stage test, not a one-target IC or DLP speedup claim.

## Promotion gate

Keep the candidate opt-in until the frozen primary elimination ratio reaches 1.20 with the paired 95% lower bound above one, complete F5-call ratio and lower bound exceed one, both holdout primary ratios clear their A/A ranges, and no smaller complete-call median regresses below its own A/A range. If one-thread gates pass, repeat the paired protocol on four pinned CPUs before selecting a default. Preserve failures and negative results. This gate is fixed before the CI measurement.
