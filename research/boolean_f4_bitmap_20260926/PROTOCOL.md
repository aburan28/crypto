# Bitmap membership for Boolean F4 symbolic preprocessing

## Hypothesis

Boolean F4 repeatedly tests and records monomial masks in `examined` and `no_divisor` hash sets while building each symbolic-preprocessing matrix and during final interreduction. When monomials use at most 22 variables, a bit indexed by the mask costs at most 512 KiB and gives direct membership. The existing column index already uses a dense array at this width. Replacing these two sets with a bitmap should reduce build work without changing the monomials, reducers, counters or final basis. Other sets and systems above 22 variables keep their existing hashed path.

## Frozen comparison

- Source base: `5f993b6a3ee73a802b5c0a2f601bd80a8a2c1105` (`origin/main` after the negative F5 experiment). Use the same release binary with `F4_F2_BITMAP_SEEN=0` for the reference and `=1` for the candidate; the original path remains the default until promotion.
- Run `examples/f4_f2_bench.rs` with one repetition, maximum 20 variables and family `f4`. The original seed XOR zero is the frozen suite; new holdouts use XORs `000000001ac0ffee` and `000000002468ace0`. The primary cell is `n20_m30`, with `build_ms` the direct target and `wall_ms` the complete F4 stage control. Keep all other F4 cells as correctness and regression controls.
- On a single pinned Ubuntu x86-64 runner with `RAYON_NUM_THREADS=1`, use one warmup per arm, five A/A reference pairs, and five alternating-order A/B pairs for each workload. Record every call including failures, timeouts and OOMs; limit each call to 120 seconds. Check exact basis fingerprints, basis length, steps, divisor-test and word-XOR counters, and matrix dimensions before computing ratios. Record source and binary SHA-256, CPU model/features, memory, load and affinity.
- Report the A/A ratio range, paired median/minimum and five-pair exact bootstrap 95% interval. A stage gain is accepted only if the primary build ratio is at least 1.05 with lower interval bound above 1, the full F4-call ratio and lower bound exceed 1, both holdout primary cells clear their own A/A noise, and no smaller cell regresses beyond A/A noise. Keep a negative or inconclusive result as a raw receipt and note.

This is an internal solver-stage experiment. It does not imply a one-target online IC or DLP speedup, and it leaves the rho reference, IC candidate IDs and scoreboard figures unset.

## Parallel control frozen after the one-thread result

Run `36294892907` cleared the one-thread stage gate. Before changing the default, rerun the same frozen and holdout workloads at `RAYON_NUM_THREADS=4`, pinning the process to four allowed CPUs instead of one. Keep one warmup per arm and five A/A and five A/B pairs per workload. Exact outputs and counters must still match. The bitmap must not regress the complete `n20_m30` F4 call beyond the four-thread A/A noise on any workload; retain the other six cases as regression controls. A four-thread gain is not required or claimed. The runner can differ from the one-thread run, so decide from paired ratios within each run only. If this control passes, enable the bitmap by default up to 22 variables, keep `F4_F2_BITMAP_SEEN=0` as the reference/rollback, and repeat both thread counts with explicit reference and candidate modes before merging.
