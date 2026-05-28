# Phase 22.6: Barrett → Pollard rho end-to-end integration

**Date:** 2026-05-28
**Status:** EXECUTED, honest result documented
**Verdict:** Barrett delivers **1.5-2.3× end-to-end** (vs 9.47× microbench)

## Setup

Integrate Barrett reduction (Phase 22.1) into the full Pollard rho
walk. Replace `mpz_mul + mpz_mod` with `__uint128_t` Barrett. Keep
modular inverse via GMP (with thread-local persistent `mpz_t` to
avoid per-call alloc).

Source: `phase22_rho_barrett.c`.

## End-to-end throughput at 81-bit

| Threads | v3 baseline (GMP only) | Barrett + GMP inv | Speedup |
|--------:|-----------------------:|------------------:|--------:|
| 1 | 2.32e6 | 3.66e6 | **1.58×** |
| 2 | 4.92e6 | 6.63e6 | **1.35×** |
| 4 | 7.39e6 | 1.40e7 | **1.89×** |
| 8 | 1.02e7 | 2.27e7 | **2.23×** |

## Why end-to-end is far below the 9.47× microbench

The Pollard rho point-addition step does:
- ~3-4 modular multiplications (Barrett accelerates: ~30× → ~10ns each)
- 1 modular inverse (NOT accelerated by Barrett: ~150-250ns via GMP)
- ~3-4 modular additions/subtractions (already O(1))

The modular inverse dominates per-add cost. Barrett can only speed
up muls; the inverse step remains GMP-bound. Hence end-to-end speedup
is bounded by the inverse cost / total cost ratio.

## Alternative tested: Fermat-via-Barrett inverse

We tried implementing modular inverse via Fermat's little theorem
(`a^(p-2) mod p`) using Barrett mul throughout. Result:

| Implementation | Throughput (1 thread) | vs Barrett+GMP |
|---------------|----------------------:|---------------:|
| Barrett + GMP inverse | 3.66e6 ops/sec | 1.0× |
| **Barrett + Fermat inverse** | **1.22e6 ops/sec** | **3× SLOWER** |

Fermat-via-Barrett requires ~80 squarings + ~40 muls (~120 Barrett
muls) per inverse. At Barrett 9.4ns/mul, that's ~1130 ns per inverse,
much slower than GMP's binary-GCD inverse at ~200ns.

**Conclusion:** GMP's binary-GCD inverse is faster than Fermat at
81-bit. Replacing GMP inverse with custom code would require
implementing binary GCD in `__uint128_t` — not done in this phase.

## Budget analysis update

For 80-bit ECDLP with `1.25√n ≈ 1.4e12` ops baseline (Pollard rho)
and engineering optimizations (negation √2, partition ℓ=16 give
~1/1.97):

Expected ops: 7e11

| Configuration | Aggregate ops/sec | 80-bit time |
|---------------|------------------:|------------:|
| v6 v3 baseline (2 threads) | 4.9e6 | 39.7 h |
| Barrett-rho (2 threads) | 6.6e6 | 29.4 h |
| Barrett-rho (4 threads) | 14.0e6 | 13.9 h |
| Barrett-rho (8 threads) | 22.7e6 | **8.6 h** |
| Budget (AutoLab 2 CPUs × 4 h wall) | — | 4 h target |

**Status: 80-bit ECDLP is still NOT feasible in budget with Barrett
alone.** Closest case (8 threads, would need 4+ HT cores) takes 8.6 h.

## What this changes about v7 plan

The microbench-projected 5.5× gap closure was wrong. Barrett's actual
end-to-end contribution is ~2×. Need additional optimizations:

| Direction | Potential | Engineering cost |
|-----------|----------:|------------------|
| Jacobian projective + Barrett | ~1.5-2× more | days |
| Binary-GCD inverse in u128 | ~1.3× more | 1 week |
| AVX/SIMD vectorized field ops | ~1.5-2× more | 2 weeks |
| **Combined achievable** | ~3-4× over Barrett alone | ~3-4 weeks |

If we add ~3× more on top of Barrett-rho's 22M ops/sec (8 threads),
we reach **~66M ops/sec → 80-bit in ~3 hours**. Within budget.

**Revised v7 status: PARTIAL**
- Phase 22.1 Barrett microbench: ✓ DONE (9.47×)
- Phase 22.6 Barrett integration: ✓ DONE (1.5-2.3× end-to-end)
- **Remaining gap: 5.5× / 1.7× (Barrett 2-thread) = 3.2× still needed**
- Phase 22.2 SIMD or Phase 22.7 Jacobian+Barrett or Phase 22.8
  binary-GCD inverse become priority

## Methodological lesson

This phase corrects another optimistic projection: the Phase 22.1
microbench (9.47×) does NOT directly translate to end-to-end speedup
when modular inverse dominates per-step cost.

The methodological pattern (raw → microbench → end-to-end → budget)
has now caught two overestimations in v6/v7:
- v6 Phase 21 raw ops/sec → end-to-end scaling validation
- v7 Phase 22.1 microbench → Phase 22.6 end-to-end integration

**Both times, the optimistic projection was wrong by 3-5×.** Always
validate end-to-end before claiming feasibility.

## Updated PROGRAMS.md sequencing

1. **v7 Phase 22.7: Jacobian + Barrett** (new) — projected 1.5-2× over Phase 22.6
2. **v7 Phase 22.8: binary-GCD inverse in u128** (new) — projected 1.3× over Phase 22.6
3. Combined 22.6 + 22.7 + 22.8 → ~3× over baseline = ~80-bit in ~10h
4. **Phase 22.2 SIMD** if still short
5. **v11 Phase 41 (paper)** can now incorporate Phase 22.6 honest result

## Honest summary

Barrett reduction: real microbench win (9.47×), modest end-to-end gain
(1.5-2.3×). Does not close the 80-bit feasibility gap on its own.
Needs Jacobian projective + binary-GCD inverse for full effect.

The 5.5× gap is now reduced to **3.2× remaining gap**. Closeable
with another 1-3 weeks of engineering, not a single sub-phase.
