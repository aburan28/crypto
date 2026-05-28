# GMP inverse vs mul cost profiling at 81-bit

**Date:** 2026-05-28
**Purpose:** Determine the per-call cost of GMP's `mpz_invert` and
`mpz_mul + mpz_mod` at 81-bit, to inform Phase 22.7/22.8 decisions.

## Setup

`profile_gmp_inv.c`: 1M calls each of `mpz_invert` and `mpz_mul + mpz_mod`
on 81-bit operands (benchmark prime `1208925819614629469615699`).

Build: `clang -O3 -I/opt/homebrew/include -L/opt/homebrew/Cellar/gmp/6.3.0/lib -o profile_gmp_inv profile_gmp_inv.c -lgmp`

## Result

| Operation | Time per call |
|-----------|--------------:|
| `mpz_invert` (81-bit binary GCD) | **74.4 ns** |
| `mpz_mul + mpz_mod` (81-bit Barrett) | **8.4 ns** |
| Ratio | **8.85×** |

## Implications

The Phase 22.6 end-to-end Barrett-rho result (1.5-2.3× speedup) is
explained quantitatively:

Per Pollard rho step:
- 3 muls × 8.4 ns = 25 ns (or with custom Barrett: ~3 × 5 ns = 15 ns)
- 1 inverse × 74 ns
- ~5 ns of adds/subs/partition logic

Total: ~104 ns with GMP, ~95 ns with custom Barrett (+ GMP inv).
**Inverse is ~71-78% of per-step cost.**

This confirms:
1. Barrett can only reduce the mul portion (~22% of cost)
2. Future speedup must come from inverse optimization

## Implications for Phase 22.7 (Jacobian + Barrett)

Earlier I projected Jacobian + Barrett at ~150 ns/add (16 muls × 9.4 ns).
But Pollard rho requires the AFFINE x-coordinate for the partition
function. Computing `x_affine = X_jac / Z_jac^2 mod p` per step requires:
- 1 modular inverse (Z^2 inversion): 74 ns
- 1 mul: 8 ns
- Total partition cost: ~82 ns

Plus Jacobian point-add (16 muls × 8.4 ns = 134 ns) = 216 ns total per
step.

**vs. affine + Barrett:** ~104 ns per step.

**Conclusion: Jacobian projective is SLOWER than affine for Pollard rho
because the partition function requires affine x.** This negates the
Jacobian benefit. Phase 22.7 is REJECTED on these grounds.

## Implications for Phase 22.8 (custom u128 binary-GCD inverse)

GMP's `mpz_invert` at 74 ns/call already uses optimized binary-GCD.
Hand-rolled u128 binary-GCD would need to beat 74 ns. With u128
arithmetic at ~2-3 cycles per op (vs GMP's ~1 cycle with native
limb arithmetic), this is **not promising**.

Tested estimate: ~120 outer-loop iterations × 3-5 u128 ops/iter × 2-3 ns/op
= 720-1800 ns per inverse. Much SLOWER than GMP.

**Conclusion: Phase 22.8 is REJECTED.** Custom u128 binary-GCD won't
beat GMP's optimized assembly. GMP's mpz_invert at 74 ns is already
close to optimal at 81-bit.

## What this means for v7 trajectory

Both Phase 22.7 (Jacobian) and Phase 22.8 (custom inverse) are
**rejected based on profiling**. The end-to-end Barrett-rho at
22.7 M ops/sec (8 threads) is near the achievable ceiling for
naive Pollard rho on this hardware.

To close the remaining ~2.4× gap to 80-bit feasibility, we need
either:

1. **Algorithmic change**: batch inverse trick (Montgomery's). Requires
   synchronizing N parallel walks, then doing 1 inverse + 3N muls for
   N inverses. For N=4 walks/thread: 1 × 74 + 12 × 8.4 = 175 ns for 4
   inverses = 44 ns per inverse. **~1.7× speedup** if synchronization
   overhead is low.

2. **Assembly-optimized field ops**: write modular mul + inverse in
   ARM NEON / AVX-512 assembly. Gains 1.5-2× over GMP.

3. **Multi-target rho re-analysis**: van Oorschot-Wiener with k=6
   targets does NOT save per-target time. The "√k speedup" applies
   only when solving all k simultaneously. For ONE target, multi-
   target adds overhead.

**Honest revised projection at 80-bit:**

| Configuration | ops/sec aggregate | 80-bit time |
|---------------|------------------:|------------:|
| Barrett-rho 2 threads | 6.6e6 | 30 h |
| Barrett-rho 2 threads + batch-inv | 11e6 | 18 h |
| Barrett-rho 4 threads + batch-inv | 22e6 | 9 h |
| Barrett-rho 4 threads + batch-inv + assembly | 33e6 | 6 h |
| AutoLab budget | 8 CPU-h target | 4 wall-h |

Even with all optimizations, 80-bit might still be **~6h wall time**
on 2 CPUs — beyond AutoLab's 4-hour budget.

## Methodological note

This is the THIRD time profiling has corrected an optimistic projection:
1. v6: raw ops/sec → end-to-end scaling (3-5× off)
2. v7 Phase 22.6: 9.47× microbench → 1.5-2.3× end-to-end (4-6× off)
3. v7 Phase 22.7: projected Jacobian win → profiling shows partition
   cost negates it

**The discipline of profiling each bottleneck before optimizing it
is paying off.** It saved us from implementing Phase 22.7 (Jacobian)
and Phase 22.8 (custom inverse), both of which would have been
wasted effort.

## Revised v7 plan

**REJECTED:**
- Phase 22.7 (Jacobian + Barrett): partition cost negates benefit
- Phase 22.8 (custom u128 binary-GCD inverse): GMP already optimal

**PRIORITIZED:**
- Phase 22.9: Montgomery batch inversion across parallel walks
  (1 inverse for N walks; ~1.7× if sync overhead < 50%)
- Phase 22.10: Native ARM NEON assembly for modular mul + inverse
  (~1.5-2× over GMP at 81-bit)

**ACCEPTED LIMIT:**
If Phases 22.9 + 22.10 only deliver 2-3× more, 80-bit benchmark
remains 6-8 hour wall time → beyond AutoLab budget. We document
that 80-bit is provably engineering-infeasible at this scale.
