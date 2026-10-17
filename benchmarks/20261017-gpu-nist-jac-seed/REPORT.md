# gpu/nist rho: inversion-free Jacobian seeding ladder

2026-10-17. Optimisation of the **walk-seeding / jump-table-build** path of
the `gpu/nist` Pollard-rho engine (`EcCurve::scalar_mul` /
`double_scalar_mul`). This is the per-walk setup `a·G + b·Q`, not the hot
walk step; the step is unchanged.

## What changed

The seeding ladder was affine double-and-add: every doubling and every
conditional point addition computed its own field inverse. For a
256-bit interleaved `a·G + b·Q` that is on the order of **~750 field
inversions per seed** (≈256 doublings + ≈256+256 conditional adds, one
Itoh–Tsujii / Fermat inverse each).

The NIST curves are `a = -3`, so the field header already carries the
inversion-free `a=-3` Jacobian formulas (`n_double`, `n_madd`). The ladder
now runs in Jacobian coordinates and converts to affine **once** at the
end: **1 field inversion per seed** instead of ~750. A general-`a` curve
(the test's toy) keeps the affine ladder, selected by the curve's
`AP::is_a_minus3()` policy.

### Correctness: the mixed-add dispatch that makes it exact

`n_madd` is the generic-position mixed addition; fed two equal points it
returns garbage with `Z = 0`. In the real walk `Q` is a scalar multiple of
`G`, so the accumulator `α·G + β·Q` coincides with `±G` or `±Q` **often**
(a measured ~2% of random 256-bit seeds), not negligibly. A naive Jacobian
ladder is therefore wrong, not just slow. The new `jmadd` dispatches the
equal case to a doubling and the inverse case to the identity before
calling `n_madd`, so the Jacobian result is **bit-identical** to the affine
ladder on every input (0 mismatches over 3000 random `Q = 7G` seeds, the
hardest dependent case, and over the 400-seed timed run).

## Measured (operation-count-dominated microbenchmark)

Standalone `g++ -O3` microbenchmark, `double_scalar_mul(G, a, Q, b)` with
random full-width `a, b` (source in the study directory). The ratio is
structural — one inverse versus ~750 — so it is reported as a
reduction, corroborated by wall time; it is **not** an isolated,
PMU/IPC-grade timing (shared cloud VM, L0, no isolation receipt).

| curve | ms/seed Jacobian (new) | ms/seed affine (old) | speedup | mismatches |
|:--|--:|--:|--:|--:|
| P-256 | 0.6615 | 11.5295 | 17.43× | 0 |
| P-384 | not timed | not timed | — | 0 (test-verified) |

Field inversions per seed: ~750 → 1 (P-256), ~1150 → 1 (P-384). P-384 is
verified bit-identical by the test but not separately timed; its affine
arm runs minutes per batch on this shared VM and the op-count story is the
same as P-256.

## Scope and honesty

- This speeds **setup only**: the initial seed of each walk and each
  reseed after a distinguished point or `max_steps`. Reseeds are rare
  relative to walk steps (one per ~2^dp_bits steps), so the effect on
  whole-run GPU throughput is small; it is a clear win for the host-side
  jump-table build and for short/high-reseed-rate runs, and it removes
  ~750 serial inversions from each device reseed.
- The change is in the shared header used by both the CPU test and the
  CUDA kernels, so the device table-build benefits too. No GPU timing is
  claimed here; the device throughput figures are in `gpu/nist/RESULTS.md`.
- No hardware-isolated A/B, no bootstrap interval: this is an op-count
  reduction with a corroborating wall-clock ratio, labelled L0.

## Reproduce

`gpu/nist/test_rho.cpp` (suite "group laws + Jacobian seeding ladder ==
affine ladder") asserts the Jacobian path equals the affine one on P-256
and P-384; `make -C gpu/nist test` runs it. The timing microbenchmark
source and raw output are in this directory (`sb2.cpp`, `data.csv`).
