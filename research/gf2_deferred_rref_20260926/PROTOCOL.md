# Deferred above-row reduction in the shared GF(2) Macaulay kernel

## Hypothesis

The matrix-F5 degree-4 step on the frozen `f4_f2_bench` n24/m24 system spends much of its elimination time revisiting rows above every pivot block. A forward echelon pass followed by reverse block reduction should preserve the bit-for-bit reduced row echelon form while reducing memory traffic enough to improve the F5 elimination phase. The intended target is at least a 1.2x median stage gain; the estimated 1.5–2x remains unmeasured until this run.

## Frozen inputs and reference

- Source revision: `a71eac61ebc2b28529be2d367a76e0fcc0e56c37` (`origin/main` as locally available on 2026-09-26).
- Benchmark: `examples/f4_f2_bench.rs`, unchanged; planted random quadratic systems use its fixed seeds and emit row fingerprints, rank, phase timings and counted word XORs.
- Primary cell: `f5_n24_m24_d4`; smaller F5 cells are correctness and regression controls. Use the same `RAYON_NUM_THREADS` for both arms.
- Reference: the same binary with `KIC_GF2_DEFER_ABOVE=0`, which runs the prior clear-above-and-below ordering. This control was selected before timing because the unmodified release build could not finish under host contention. The original source revision above remains the source reference.
- A/A: five alternating executions of the reference mode against itself. A/B: at least five alternating reference/candidate pairs (`KIC_GF2_DEFER_ABOVE=1`). Compare median and per-pair ratio; wall time is a stage diagnostic. The candidate remains opt-in until these gates pass.

## Correctness and stop conditions

Both arms must return the same rows fingerprint, rank, row/pruning counts and reduced matrix bit for bit on the kernel's adversarial and random tests. Counted `reduce_word_ops` may differ because the algorithm changes its table construction; report it, never force equality. Stop and reject if any output differs, a relevant test fails, the primary elimination phase regresses beyond A/A noise, or the host is too contended for a useful wall comparison. Preserve negative results.

## Cost boundary

This is an internal Macaulay elimination experiment. Charge table construction, all forward and reverse reductions, criterion, row build and unpack in the F5 step. Record the full F5 step time alongside the elimination phase. No IC candidate ID, one-target DLP, rho comparison or end-to-end speedup is claimed from this solver-stage result.
