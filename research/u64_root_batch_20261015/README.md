# u64 S3 root-solve batching in the Koblitz orbit fast path — 2026-10-15

## Change

`examples/koblitz_orbit_dlp_fast.rs`, u64 pipeline (field degrees n ≤ 63):
the `extract` scan previously paid one Itoh–Tsujii field inversion (≈10
multiplies + table applications) **per candidate** via
`S3Solver::roots`. This change ports the wide twin's `RootBatch128`
Montgomery-trick batching down to the u64 path:

- `S3Solver::prepare` / `S3Solver::finish` split (mirrors
  `S3Solver128`), with `roots` delegating so the scalar form is unchanged.
- New `RootBatch` (u64): per state, all `n × |state.normal_roots|` candidate
  root solves share **one** inversion over the non-exceptional denominators.
- `extract` collects the state's `(absolute, target_x)` inputs, batch-solves,
  then walks the results in the original candidate order.

Candidate order, published relations, probe counts and recovered scalars are
unchanged — pinned by a new unit test
(`root_batch_tests::batched_roots_match_scalar_roots_in_order`, 4,096 inputs
including the `(x,x)`, `(0,y)`, `(x,0)` exceptional paths, chunked, compared
against scalar `roots`).

Class: **engineering** (constant factor; the counted probe behaviour and
outputs are identical). Baseline: `origin/main` at
`9f7fb66e6b` ("Merge pull request #1956"), the current best verified
index-calculus pipeline for this example.

## Measurement A/B

Matched baseline/candidate whole-process runs on identical inputs. Host:
Apple M4 Pro, 14 cores, macOS, rustc 1.93.1. Shared host (other agents'
jobs ran concurrently; the baseline arm shows large rep-to-rep spread —
83.9/76.4/41.9 s — and one contended spread on the candidate too). Run
order interleaved base/candidate within each repetition.

```
cargo run --release --example koblitz_orbit_dlp_fast -- \
  construct:41:0:26 scalars41.txt 20261015 out.jsonl
```

16 public synthetic known-answer scalar fixtures at n=41, a=0, columns=26,
rank seed 20261015, 1 thread. Both arms solve the same deterministic 11/16
targets (base-size artifact of the small column count, identical on both
arms). Recovered scalars, `point_indices`, `x_codes` and per-target `probes`
are identical on every repetition (diff-verified, `raw/`).

| arm | process_total (s), 3 reps | median | min | rank_stage median (ms) | target_median (ms) |
|---|---|---|---|---|---|
| baseline (main `9f7fb66e6b`) | 83.96 / 76.38 / 41.94 | 76.38 | 41.94 | 30 333 | 1 298 |
| candidate (this branch) | 40.76 / 23.40 / 22.42 | 23.40 | 22.42 | 10 261 | 213 |
| ratio (median) | — | **3.26×** | — | 2.95× | 6.09× |
| ratio (min, least-contended) | — | — | **1.87×** | — | — |

Wall-clock on a shared host without `isolated_bench`; treat these as stage /
whole-process wall diagnostics. The identical-output check (including
per-target probe counts) is the primary correctness gate and passes on all
reps. Under the least-contended pairing the gain is ≥ 1.87× whole-process;
both arms solve and group-verify the same targets.

## Inversion-count control

Per solved state the scalar path pays `2n` inversions minus exceptional
paths; the batched path pays 1. At n=41 that is ~80× fewer Itoh–Tsujii
inversions in the extract scan; the residual scan cost is the S3
prepare/finish multiplies, table applications, canonicalization and root-table
probes — which matches the observed bound (the scan is not 80× of end-to-end
because probes and lifts dominate after the change).

## Correctness

- `cargo test --release --example koblitz_orbit_dlp_fast`: 21/21 (new
  batched-vs-scalar equivalence test included).
- End-to-end n=41 A/B: identical recovered scalars and relations, both arms,
  3/3 reps; `targets_failed` identical.

## Next

- The same batching is available in the wide twin's table-build loop where it
  is not yet applied (`construct_base128` scalar-mul path remains the setup
  bottleneck).
- Precomposed twisted-linear byte tables (`to_poly ∘ rotate_k ∘ to_normal`)
  for both pipelines.
