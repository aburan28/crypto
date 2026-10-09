# Full-width Koblitz IC coefficient gate

## Hypothesis and frozen inputs

The legacy native Koblitz IC relation collector takes the low `u64` limb of
the subgroup order before drawing `a` and `b` for `[a]G + [b]Q`. On
`icv1-f2m83-tm6151469093347-debefd74` (the E₀ m = 83 gate), the
prime subgroup order is `2417851639230796216685689`, so that sampler
cannot draw most valid coefficients. This is a correctness and workload
distribution repair, not a performance experiment.

Freeze `KoblitzCurve::known_n83_k0()` and that order; use seeds 0, 1, 5,
and 0x83. The unmodified code at the parent commit is the reference for
curves with subgroup order at most `u64::MAX`. The changed code must produce
the exact same coefficient sequence there, preserving committed small-curve
replays. For m = 83, deterministic rejection sampling must draw only
`1 <= a,b < r`, reach coefficients above `u64::MAX` over the frozen seeds,
and reproduce the same sequence on a second invocation with the same seed.

## Acceptance and stop conditions

- Pass a native Rust test of the small-order sequence against direct
  `StdRng::gen_range(1..r)` calls, and a full-width test of bounds,
  reproducibility, and high-limb coverage for the frozen m = 83 order.
- Keep `[a]G + [b]Q` and relation-matrix arithmetic unchanged. Preserve
  existing rank and verification checks. An observed recovered scalar must
  still pass `[d]G = Q` independently.
- Stop this patch at coefficient sampling. If another full-solve path still
  truncates the order, or factor-base construction, PDP, matrix work, or
  final descent cannot run at m = 83 within a frozen budget, leave the
  m = 83 IC comparison pending. Do not infer a full solve from sampler
  tests, a relation yield, or the already recorded wide rho step rate.

## Cost and later comparison

This patch records no speedup, `S` change, or new wall figure. The
required later gate uses matched public targets, the same factor-base
policy and budgets for unmodified baseline and candidate, a strong rho
arm on the exact subgroup, complete cold cost including native SAT/PDP
work and memory, failures and budget exits, replay receipts, and the
canonical scoreboard. Counted GAE is a lower bound until foreign work
has a measured conversion; wall time needs L2 or L3 and an A/A floor.

## Execution record

On the combined validation tree with parent `c191b70be`,
`rustc 1.98.0 (88d9e12ae 2026-08-18)`, the native command
`cargo test --release --lib relation_scalar_sampler_ -- --nocapture`
passed both tests (2 passed, 0 failed). The test checks exact narrow-stream
identity over 384 draws and full-range bounds, replay, and high-limb
coverage over 512 m = 83 draws. This is a deterministic correctness check;
it supplies no timing or full-pipeline IC evidence.
