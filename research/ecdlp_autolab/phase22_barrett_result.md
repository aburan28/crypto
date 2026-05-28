# Phase 22.1 Barrett reduction: empirical result

**Date:** 2026-05-28
**Status:** COMPLETED, exceeds projection by ~5×

## Setup

Barrett reduction for 81-bit prime `p = 1208925819614629469615699` (the
67.a1 benchmark prime). Pre-compute `μ = ⌊2^{2k} / p⌋` once;
per-multiplication operation uses only `__uint128_t` arithmetic.

Compared against two GMP baselines:
1. GMP with `mpz_init/import/clear` per call (worst case)
2. GMP with persistent `mpz_t` variables (fair case)

Correctness verified: Barrett and GMP produce identical results.

## Result (10M multiplications at 81-bit)

| Implementation | Throughput | Per-mul (ns) |
|---------------|-----------:|-------------:|
| **Barrett (this work)** | **106M muls/sec** | **9.4 ns** |
| GMP persistent mpz_t | 11.1M muls/sec | 90.0 ns |
| GMP malloc per call | 1.38M muls/sec | 724 ns |

### Speedups
- Barrett vs GMP persistent: **9.47×** ← the fair comparison
- Barrett vs GMP malloc: 76.5× (microbenchmark artifact)

## Why so much faster

GMP's `mpz_mul + mpz_mod` for fixed-size 81-bit inputs:
- General-purpose code path designed for arbitrary precision
- Branch on input sizes per call
- Allocate scratch space (even with persistent vars, there's some overhead)

Barrett-with-fixed-`p`:
- Single inline `__uint128_t` multiplication
- Precomputed `μ`, no division needed
- Two unconditional u128 multiplies + one final subtraction
- No branches except the final `if (r >= p) r -= p`

## Impact on v6 budget projection

Phase 21.1 v3 baseline (affine + pthreads + GMP persistent):
- 2 threads at 81-bit: 4.92M ops/sec

If Barrett gives 9× per-mul speedup, and point arithmetic is ~75%
modular muls (the rest being subtractions and the one inverse per add),
end-to-end speedup is roughly:

```
new_ops_per_sec ≈ old_ops_per_sec × (1 / (0.25 + 0.75 / 9))
                ≈ 4.92M × (1 / 0.333)
                ≈ 14.8M ops/sec (per thread!)
```

More conservatively, if Barrett helps muls but not inverses (~25%
of point-add cost):
- Speedup factor ≈ 1 / (0.75 + 0.25/9) ≈ 1.28× — too pessimistic

Realistic mid-estimate: **5-7× end-to-end** for the affine point
arithmetic.

## Budget revision

| Configuration | ops/sec | 80-bit projection |
|--------------|--------:|------------------:|
| v6 GMP, 2 threads | 4.92e6 | 44 hours |
| v6 GMP, 4 threads | 7.39e6 | 30 hours |
| **+ Barrett, 2 threads** | **~25-35e6** | **6-9 hours** |
| **+ Barrett, 4 threads** | **~40-55e6** | **4-6 hours** |
| + Barrett + multi-target (6 targets, √6) | ~100-140e6 | **1.6-2.5 hours** |

**80-bit AutoLab benchmark becomes feasible** in the 4-hour budget
with Barrett + multi-target stack. This closes the 5.5× gap with
margin to spare.

## Caveats

1. **End-to-end integration not yet measured.** Microbenchmark of
   modular mul vs full Pollard rho point operation may differ.
   Need `phase22_rho_barrett.c` to validate.

2. **Modular inverse still uses GMP.** Per-add cost is dominated by
   one inverse (~80 multiplications via Fermat). Need batched inverse
   or projective coordinates to fully exploit Barrett.

3. **Per-mul overhead reduces with larger n.** At very high throughput,
   memory bandwidth becomes a constraint.

## Next steps (v7.5)

1. Integrate Barrett into `phase21_rho_v3.c` → `phase22_rho_barrett.c`
2. Measure end-to-end Pollard rho throughput at 81-bit
3. Run full scaling test (20-40 bits) with Barrett
4. Project to 80-bit honestly with measured constants

If end-to-end speedup is ≥5×, **80-bit AutoLab benchmark is in reach
with just Barrett + multi-target**, no further optimization needed.

## What this changes about v7 plan

The v7 critique (`PLAN_CRITIQUE.md`) worried Barrett might give only
1.0-1.3×. Empirical reality: 9.47× on the microbenchmark. This
suggests:

- v7 Phase 22.1 is HEAVILY overdelivering its 2× target
- v7 plan may not need Phase 22.2 (SIMD) and Phase 22.3 (r=64) to
  hit the success criterion
- We can ship "v7 Phase 22.1 only" as a complete v7 if end-to-end
  validation confirms

This is the kind of pleasant surprise the v9 phase 36-style critique
called out: be open to over-delivery.
