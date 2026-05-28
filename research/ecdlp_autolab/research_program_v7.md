# Research Program v7: Engineering to close the 80-bit gap

**Date:** 2026-05-28
**Status:** Phase 22.1 EXECUTED (overdelivered), Phase 22.2-22.5 deprioritized
**Motivation:** v6 demonstrated a 23-86× C-extension speedup but a **5.5×
gap remains** to make 80-bit AutoLab ECDLP solvable in the budget. v7
targets that gap with concrete engineering work, each item independently
testable.

## Major update (2026-05-28): Phase 22.1 overdelivered

Barrett reduction prototype: **9.47× speedup over GMP persistent mpz_t**
(106M vs 11.1M muls/sec). See `phase22_barrett.c` and
`phase22_barrett_result.md`.

This single result, if it carries through to end-to-end Pollard rho,
**closes the 5.5× gap by itself** (projected 80-bit ECDLP in ~2-6
hours on 2 CPUs). The other v7 phases (22.2-22.5) may be unnecessary.

**Revised v7 plan:**

1. **Priority 1 (next):** Phase 22.6 = Barrett integration into
   Pollard rho walk + end-to-end ECDLP validation at 40-60 bits
2. **Priority 2 (if 22.6 insufficient):** Phase 22.2 SIMD, 22.3 r=64
3. **Optional/Deferred:** Phase 22.5 tag walk (small gain)

If Phase 22.6 validates end-to-end at the projected speedup, v7
**ships** with just Phase 22.1 + 22.6. The other sub-phases become
v7.5 / v8 territory.

## Goal

Bring `phase21_rho_v3.c`'s 4.92e6 ops/sec (2 threads) up to **>2.7e7 ops/sec
on 2 CPUs**. That closes the 5.5× gap and makes 80-bit feasible.

## Sub-phases

### Phase 22.1: Montgomery + Barrett reduction (target: 2×)

**Why:** GMP's `mpz_mod` is generic; it doesn't exploit that `p` is fixed
across millions of multiplications. With precomputed Barrett constant
`μ = ⌊2^{2k}/p⌋` and a small fixed-size representation, we save the
generic-prime overhead.

**Concrete experiment:**
1. Implement Barrett reduction for 81-bit `p` with `k = 81`:
   - Precompute `μ` once at startup
   - For each mul, compute `q = (a*b * μ) >> 162`, then `r = a*b - q*p`
2. Use `__uint128_t` for intermediate products
3. Compare to `phase21_rho_v3.c`'s `mpz_mul + mpz_mod` baseline

**Deliverable:** `phase22_barrett.c`. Hypothesis: 2× speedup (~10M ops/sec
single-thread).

### Phase 22.2: 6×6 SIMD multi-target rho (target: 1.5-2× on top of 22.1)

**Why:** AVX-2/NEON SIMD can do 4-8 wide multiplications per cycle. Each
of the 6 LMFDB benchmark targets is independent; pack 6 walks into one
SIMD lane and run them in lockstep.

**Concrete experiment:**
1. Pack 4 (NEON) or 6 (AVX-2) parallel walks into one process
2. Each walk does its own 81-bit point arithmetic in parallel lanes
3. Distinguished-point detection across lanes
4. Per-walk trail collection

**Deliverable:** `phase22_simd.c`. Hypothesis: 1.5-2× over Phase 22.1
single-thread.

### Phase 22.3: Larger partition (r=64 instead of r=16) (target: 1.2×)

**Why:** Brent's analysis shows partition size `r` reduces the expected
walk length by factor `~r/r-2` for large `r`. r=64 gives ~1.2× over r=16.

**Concrete experiment:**
1. Modify `phase21_rho_v3.c` to use r=64 partition table
2. Precompute 64 multipliers `M[0..63] = a_i*P + b_i*Q` for chosen `a_i, b_i`
3. Benchmark walk length to distinguished point at fixed `d_bits`

**Deliverable:** `phase22_r64.c`. Hypothesis: 1.2× speedup, smaller table
cache miss tradeoff to study.

### Phase 22.4: Hyperthreading benchmark (target: 1.3-1.5× on 2-CPU)

**Why:** AutoLab provides 2 physical CPUs. Hyperthreading might give us
4 effective threads with diminishing returns. Benchmark to see.

**Concrete experiment:**
1. On a 2-CPU machine matching AutoLab's, run `phase21_rho_v3` with
   1, 2, 3, 4 threads
2. Plot scaling; calibrate the effective parallelism

**Deliverable:** `phase22_hyperthreading.md` with measurement.

### Phase 22.5: Endpoint-aware walk function (target: 1.1×)

**Why:** Tag-based walk functions (van Oorschot-Wiener) reduce the
expected walk length to distinguished point.

**Concrete experiment:**
1. Implement walk with linear-feedback shift register tag
2. Distinguished point: tag mod 2^d == 0
3. Compare collision rate to current implementation

**Deliverable:** `phase22_tag_walk.c`.

## Phase 22.6: Barrett → Pollard rho integration ✓ EXECUTED

**Status:** Done. Result: **1.5-2.3× end-to-end**, well below the
projected 5× target.

**Why result fell short of microbench:**
- Pollard rho point-add does ~3 muls + 1 inverse + few add/sub
- Barrett accelerates muls (~30× faster)
- Modular inverse (GMP) NOT accelerated, dominates per-step cost
- Net per-add cost dropped from ~430ns to ~270ns (~1.6× win)

**Result table (81-bit):**

| Threads | v3 (GMP only) | Barrett+GMP inv | Speedup |
|--------:|--------------:|----------------:|--------:|
| 1 | 2.32e6 | 3.66e6 | 1.58× |
| 2 | 4.92e6 | 6.63e6 | 1.35× |
| 4 | 7.39e6 | 14.0e6 | 1.89× |
| 8 | 1.02e7 | 22.7e6 | 2.23× |

**Risk hypothesis confirmed:** Barrett benefit was lost in inverse
noise. The risk mitigation kicks in: execute Phases 22.7-22.8 next.

**Failed alternative:** Fermat-via-Barrett inverse is 3× SLOWER than
GMP binary-GCD inverse. Custom u128 binary-GCD needed.

## Phase 22.7: Jacobian projective + Barrett (NEW, PRIORITY 2)

**Why:** Jacobian projective coordinates avoid modular inverse per
point-add. Per-add cost = ~16 muls (no inverse). With Barrett mul at
9.4ns, that's ~150ns per add vs affine's ~270ns → ~1.8× more.

**Concrete experiment:**
1. Adapt `phase21_rho_v2.c` (Jacobian + GMP) to use Barrett mul
2. Benchmark at 81-bit, 1/2/4/8 threads
3. Convert back to affine only at distinguished points (1 inverse
   per trail, amortized)

**Deliverable:** `phase22_rho_jacobian_barrett.c`.

**Acceptance criterion:** ≥1.5× over Phase 22.6 = ≥34e6 ops/sec at
4 threads. Combined with Phase 22.6, gives 3-4× over v6 baseline.

## Phase 22.8: Binary-GCD inverse in u128 (NEW, PRIORITY 3)

**Why:** GMP's binary-GCD inverse at 81-bit takes ~200ns. Hand-rolled
u128 binary-GCD could give ~1.3× by avoiding GMP function-call overhead.

**Concrete experiment:**
1. Implement extended binary GCD for `__uint128_t`
2. Verify correctness against GMP
3. Replace `mod_inv_gmp` in Phase 22.6 with custom

**Deliverable:** Updated `phase22_rho_barrett.c` with custom inverse.

**Acceptance criterion:** Within 1.2× of GMP inverse, ideally faster.

## Combined v7 trajectory (revised)

| Phase | Speedup | Aggregate (4 threads) |
|-------|--------:|----------------------:|
| v3 baseline | 1.0× | 7.4e6 |
| + Phase 22.6 Barrett | 1.9× | 1.4e7 |
| + Phase 22.7 Jacobian | 1.5× | 2.1e7 |
| + Phase 22.8 custom inv | 1.3× | 2.7e7 |
| + multi-target (if 6 targets shared) | 2.4× | 6.5e7 |

**Projected 80-bit time on AutoLab (2 CPUs × 4 hours):**

- 4-thread effective × Barrett + Jacobian + custom inv = 27M ops/sec
- 80-bit requires 7e11 ops with engineering stack
- Time: 7e11 / 27e6 = 26,000 s = 7.2 hours

Still ~1.8× over budget. Need one of:
- Hyperthreading effective 8 threads → 54M → 3.6 hours (within!)
- SIMD for further speedup

## Path to 80-bit feasibility

Most likely combination that works:
1. Phase 22.6 Barrett (done, 1.9×) +
2. Phase 22.7 Jacobian (projected 1.5-2×) +
3. Phase 22.8 custom inverse (projected 1.3×) +
4. Aggressive use of 2-CPU hyperthreading (effective ~4× via HT)

**Estimated total: 5-7× over v6 baseline.**
**Brings 80-bit from 44 hours to 6-9 hours — borderline feasible.**

## Original aggregate plan (now superseded by Phase 22.6)

If each sub-phase delivers its projected speedup, the multiplicative
combination is:
- 22.1 Barrett: 2× *(actual: 9.47× microbench)*
- 22.2 SIMD multi-target: 1.5×
- 22.3 r=64 partition: 1.2×
- 22.4 Hyperthreading: 1.3× (on 2-CPU effective)
- 22.5 Better walk: 1.1×

Combined projection: 2 × 1.5 × 1.2 × 1.3 × 1.1 = 5.15×

**Empirical Phase 22.1 alone delivers 9.47×, exceeding the combined
projection.** Phases 22.2-22.5 become optional.

## Risk assessment

| Sub-phase | Risk | Mitigation |
|-----------|------|------------|
| 22.1 Montgomery | Implementation bugs in mod arithmetic | Compare against GMP on millions of inputs |
| 22.2 SIMD | Platform-dependent (NEON vs AVX-2 vs AVX-512) | Target both via #ifdef |
| 22.3 r=64 | Cache miss might cost more than walk-length win | Benchmark both |
| 22.4 HT | Hyperthreading benefit highly workload-dependent | Just measure |
| 22.5 Tag walk | Theoretical analysis may not match practice | Verify expected collision rate |

## Success criterion

**v7 ships** if `phase22_*.c` aggregate benchmark reaches ≥2.7e7 ops/sec
on a 2-CPU machine matching AutoLab's hardware, AND end-to-end Pollard
rho recovers a known 80-bit ECDLP secret in <4 wall-clock hours.

## What v7 does NOT pursue

- Algorithmic novelty (deferred to v8)
- Multi-machine distributed attacks (out of session scope)
- ASICs/FPGAs (hardware out of scope)
- Quantum (Shor's not classical)

## Cross-references

- Builds on `phase21_rho_v3.c` (v6, affine + pthreads baseline)
- Validates against `phase21_validate_scaling.sage` (end-to-end scaling)
- Lower bound from `yokoyama_lower_bound.md` (algorithm cannot beat √n)

## Estimated effort

- 22.1 Montgomery: 2-3 days hand-rolled C
- 22.2 SIMD: 1 week with platform tuning
- 22.3 r=64: 1 day
- 22.4 HT measurement: 1 hour
- 22.5 Tag walk: 1 day
- **Total: ~2 weeks of focused engineering**

## What comes after v7

If v7 closes the gap: integration with AutoLab `solve.py` via cffi/ctypes
to verify against the actual benchmark.

If v7 falls short: re-evaluate. Either accept that 80-bit is infeasible
in 4-hour budget, or design v8 (algorithmic novelty) seriously.

See `research_program_v8.md` for the algorithmic-frontier plan that
complements v7.
