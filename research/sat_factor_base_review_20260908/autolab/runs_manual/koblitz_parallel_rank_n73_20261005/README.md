# Parallel guided rank at n=73 (2026-10-05)

`examples/koblitz_orbit_dlp_fast.rs` gained a `KIC_RANK_THREADS` mode:
every orbit column anchors its own relation
(`[a]G − R_j` with a scheduling-independent scalar
`mix(rank_seed, j, attempt)`), extractions run on W worker threads, and
rows insert into the same incremental echelon as they arrive.  The solved
factor-log table is the unique full-rank solution either way, so the
published target row is unchanged by the parallel path; `KIC_RANK_THREADS=1`
(the default) keeps the frozen sequential LCG path byte-identical.

## What is claimed here

- **Correctness (n=73):** the parallel run on the frozen n=73 base
  (`experiments/koblitz-single-target-n73-20261002/base_n73_K600.jsonl`)
  and the frozen public target reproduces the R1 target row exactly — same
  recovered scalar, same deterministic 4-point relation
  (`point_indices`, `x_codes`, `pinned_intermediates`), same probe count —
  with a *different* rank-relation set (per-column deterministic scalars).
- **Wall-time (n=71, uncontended, 12 threads, Apple M4 Pro):** the guided
  rank stage — the entire precompute cost of the vs_rho rungs — drops from
  the sequential 1,485.7 s (pilot R1) to **128.9 s, an 11.5× speedup
  (96% of ideal)**; rank 600, zero failures, zero rows without gain, and
  the published target row again reproduces the frozen n=71 ledger
  relation exactly (deterministic 14,554-probe relation).  The contended
  n=73 run (alongside two sequential IC arms) is retained as the
  correctness artifact; no speedup is claimed from it.

Not claimed: no change to the online interval, no asymptotic change, no
new vs_rho evidence — the online target stage is identical code and
input.  At n=73 the sequential precompute is 5.1 h; the same stage with
`KIC_RANK_THREADS=12` uncontended projects to ~27 min.

## Files

- `ic_parallel.jsonl` — target row from the parallel-rank run (n=73,
  `KIC_RANK_THREADS=10`, rank seed 7, frozen base + frozen target;
  contended).
- `summary_parallel.json` — producer summary (rank policy string,
  threads, timings).
- `comparison_vs_R1.json` — field-by-field comparison against the frozen
  R1 target row (`runs/…R1/ic.jsonl`), non-timing fields must all match.
- `n71/` — the uncontended wall-time benchmark: `ic_parallel.jsonl`,
  `summary_parallel.json`, `benchmark.json` (11.5× at 12 threads), on
  the frozen n=71 base and frozen n=71 public target.
- `README.md` — this note.

## Reproduce

```bash
KIC_RANK_THREADS=10 target/release/examples/koblitz_orbit_dlp_fast \
  experiments/koblitz-single-target-n73-20261002/base_n73_K600.jsonl \
  experiments/koblitz-single-target-n73-20261003/frozen/target_points.jsonl \
  7 ic_parallel.jsonl
```

Sequential control: same invocation with `KIC_RANK_THREADS=1`
(byte-identical to the frozen R1 producer path for the same seed).
