# Parallel target extraction at n=83 (2026-10-06) — an honest negative

`examples/koblitz_orbit_dlp_fast.rs` gained `KIC_TARGET_THREADS`: an
order-preserving parallel scan of the online target extraction.  Workers
scan disjoint blocks of the same rotated state order the sequential scan
uses; each stops at its block's first hit; the merge picks the smallest
global position, so the published relation (point indices, x-codes,
intermediates, probe count) is *identical* to the sequential first hit.
`KIC_TARGET_THREADS=1` (the default) preserves the frozen sequential
path.

## What is claimed here

- **Correctness:** on the frozen n=83 base and frozen public target, the
  12-thread run reproduces the frozen R1 target row exactly — same
  recovered scalar, same 4-point relation, same 8,845,441 probe count.
- **No speedup (the finding):** the deterministic first hit sits at
  ~0.1% of the scan space (8.8·10⁶ of ~10¹⁰ probes), inside the first
  worker's block, so the parallel wall (8.38 s) equals the sequential
  wall (7.94 s).  The online stage is **probe-rate bound** (~1.1 M
  probes/s single-core, dominated by the S₃ quadratic solve per probe),
  not scan-length bound.  The parallel extraction bounds only the
  unlucky-deep-hit tail; the real online levers are the per-probe rate
  and higher-arity relations (fewer probes per relation).

## Files

- `ic_parallel_target.jsonl` — target row from the 12-thread run
  (frozen base, frozen target, rank seed 21, parallel rank too).
- `summary_parallel.json` — producer summary (`target_threads: 12`).
- `comparison_vs_R1.json` — field-by-field comparison against the
  frozen R1 row plus the finding above.

## Reproduce

```bash
KIC_RANK_THREADS=12 KIC_TARGET_THREADS=12 \
  target/release/examples/koblitz_orbit_dlp_fast \
  experiments/koblitz-single-target-n83-20261006/base_n83_K600.jsonl \
  experiments/koblitz-single-target-n83-20261006/frozen/target_points.jsonl \
  21 ic_parallel_target.jsonl
```
