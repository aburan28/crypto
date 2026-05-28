# Research Program v7: Engineering to close the 80-bit gap

**Date:** 2026-05-28
**Status:** planned, not yet executed
**Motivation:** v6 demonstrated a 23-86× C-extension speedup but a **5.5×
gap remains** to make 80-bit AutoLab ECDLP solvable in the budget. v7
targets that gap with concrete engineering work, each item independently
testable.

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

## Aggregate plan

If each sub-phase delivers its projected speedup, the multiplicative
combination is:
- 22.1 Montgomery: 2×
- 22.2 SIMD multi-target: 1.5×
- 22.3 r=64 partition: 1.2×
- 22.4 Hyperthreading: 1.3× (on 2-CPU effective)
- 22.5 Better walk: 1.1×

Combined: **2 × 1.5 × 1.2 × 1.3 × 1.1 = 5.15×**

That closes the 5.5× gap to within ~10% — borderline feasible. If
any sub-phase overdelivers, we have margin.

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
