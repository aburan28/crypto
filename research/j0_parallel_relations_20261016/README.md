# Parallel relation collection in the j=0 staged IC driver — 2026-10-16

## Change

`j0_index_calculus_dlp_staged` collected raw relations serially (one
`find_relation_with_psi_counted_with_index` call at a time, retrying up to
`max_relation_attempts`). Relation searches are independent given the factor
base, so they now run under rayon `(0..target_raw_relations).into_par_iter()`
with results kept in **index order**: the relation sequence, per-relation trial
counts and `attempts_exhausted` accounting are exactly the sequential ones.
`RAYON_NUM_THREADS=1` reproduces sequential behaviour. Class: **engineering**
(wall scaling with cores; counted work unchanged). Baseline: `origin/main`
`569b157f3` (best current j=0 pipeline after #1946's orbit-index map hoist).

## Measurement A/B (matched, same workload)

`cargo run --release --example j0_stage_ladder -- <bits> 20261005`,
Apple M4 Pro (14 cores), shared, no isolation receipt. Two interleaved reps
per arm at default thread count; one rep per arm at `RAYON_NUM_THREADS=1`.
Trials per relation are unseeded, so per-run relation-stage wall varies ~±20%
by trial count; all runs recovered the same planted scalar and agreed with the
rho control (`ic_agrees_rho: true`), same relation count (48).

| workload | arm | relations_ms (reps) | median |
|---|---|---|---|
| 16-bit, 14 threads | baseline | 1184.9 / 864.8 | ~1025 |
| 16-bit, 14 threads | candidate | 221.3 / 246.5 | ~234 (**4.4×**) |
| 20-bit, 14 threads | baseline | 32931.0 / 46182.8 | ~39600 |
| 20-bit, 14 threads | candidate | 15412.7 / 11001.6 | ~13200 (**3.0×**) |
| 16-bit, 1 thread | baseline / candidate | 741.2 / 822.8 | +11% (within trial noise) |
| 20-bit, 1 thread | baseline / candidate | 44786.5 / 38435.3 | −14% (within trial noise) |

No 1-thread regression beyond the unseeded trial-count spread. Wall-clock
diagnostics on a shared host; stage improvement (relations phase), not a
frozen-suite end-to-end claim. External-tool comparison: not run (this is an
internal parallelization of our own relation search; the example's rho
agreement control passes on every run).

## Correctness

- `cargo test --release --lib ec_index_calculus_j0`: 15/15, including
  `staged_dlp_recovers_and_splits_stages` and `j0_ic_recovers_dlp`.
