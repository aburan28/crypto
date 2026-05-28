# v7 execution playbook: day-by-day plan

**Date:** 2026-05-28
**Purpose:** Concrete actions to execute `research_program_v7.md`'s
five sub-phases. Each day's plan has a deliverable and a measurement.

## Day 0: Baseline (1 hour)

Re-run `phase21_rho_v3.c` on target hardware (matching AutoLab's 2
CPUs) to establish the canonical baseline number. Record:

- ops/sec at 1, 2, 4, 8 threads
- 80-bit target ECDLP wall-time projection
- Hardware: CPU model, RAM, OS, compiler version

**Deliverable:** `baseline.md` with reproducible setup and numbers.

## Day 1-2: Phase 22.1 Barrett reduction (16 hours)

**Hour 1-4:** Implement Barrett reduction in isolation
- `mod_mul_barrett(a, b, p, mu)` where `mu = floor(2^{2k} / p)` precomputed
- Use `__uint128_t` for intermediate products
- Unit-test against GMP `mpz_mul + mpz_mod` on 10^6 random inputs

**Hour 5-8:** Wire into Pollard rho walk
- Replace `mpz_mul + mpz_mod` calls with `mod_mul_barrett`
- Verify correctness: walk to distinguished point matches v3 baseline

**Hour 9-12:** Benchmark
- 1, 2, 4 threads
- 100k, 1M, 10M steps
- Report ops/sec, speedup vs. v3

**Hour 13-16:** Tune
- Try Barrett with `k = 81 + 1` vs `k = 128` for register alignment
- Compare to Montgomery (alternative: implement Montgomery as
  `phase22_barrett.c` alternative)

**Deliverable:** `phase22_barrett.c`, benchmark report in
`phase22_barrett_result.md`.

**Acceptance criterion:** ≥1.5× speedup over v3 baseline (target 2×).

## Day 3-7: Phase 22.2 SIMD multi-target (40 hours)

**Day 3:** Choose SIMD platform
- AVX-2 (x86-64): 4 × `__uint64_t` lanes
- AVX-512: 8 × `__uint64_t` lanes
- NEON (ARM): 2 × `__uint64_t` lanes
- Pick based on AutoLab hardware

**Day 4:** Pack 4 (AVX-2) or 8 (AVX-512) parallel walks
- 4 curves OR 4 different starting points on same curve
- Vectorized point addition for these in parallel lanes
- Use `_mm256_*` intrinsics

**Day 5:** Distinguished-point detection across lanes
- Each lane independently checks `x mod 2^d == 0`
- Branch on mask result; extract trail data per-lane

**Day 6:** Trail collection
- Per-lane trails stored in shared hash table (thread-safe)
- Birthday collision detection

**Day 7:** Benchmark
- Compare to Phase 22.1 single-target baseline
- Measure scaling: 4 lanes vs 8 lanes

**Deliverable:** `phase22_simd.c`.

**Acceptance criterion:** ≥1.3× over Phase 22.1 (target 1.5-2×).

## Day 8: Phase 22.3 r=64 partition (8 hours)

**Hour 1-3:** Precompute 64 multipliers
- Random `(a_i, b_i)` pairs for i ∈ [0, 64)
- `M[i] = a_i * P + b_i * Q` (one-time computation)

**Hour 4-6:** Update walk function to use 64-partition
- `partition = x mod 64` (low 6 bits)
- Walk: `R_{i+1} = R_i + M[partition(R_i)]`

**Hour 7-8:** Benchmark walk length to distinguished point
- Compare expected walk length for r=16 vs r=64
- Measure cache miss rate (64 entries × 32 bytes = 2KB, fits in L1)

**Deliverable:** `phase22_r64.c`.

**Acceptance criterion:** ≥1.1× speedup (target 1.2×).

## Day 9: Phase 22.4 Hyperthreading benchmark (1 hour)

Boot up an AutoLab-equivalent VM (2 vCPUs), run `phase22_barrett.c`
(or 22.3 if better) at 1, 2, 3, 4 threads. Plot scaling.

**Deliverable:** `phase22_hyperthreading.md`.

**Acceptance criterion:** Just measure; no specific target.

## Day 10-11: Phase 22.5 Tag-based walk (16 hours)

**Hour 1-6:** Implement tag-based walk
- Each step generates a "tag" via LFSR or similar
- Distinguished point: tag mod 2^d == 0
- Less correlation between adjacent steps → faster collision

**Hour 7-12:** Benchmark expected walk length
- Compare distinguished-point hit rate to standard walk
- Measure trail length distribution

**Hour 13-16:** End-to-end validation
- Run full ECDLP at 30-40 bits with tag walk
- Verify secret recovery
- Compare wall time to standard walk

**Deliverable:** `phase22_tag_walk.c`.

**Acceptance criterion:** ≥1.05× wall-time improvement (target 1.1×).

## Day 12: Combined-stack measurement (8 hours)

**Hour 1-4:** Combine all sub-phases that gave wins
- `phase22_combined.c` = Barrett + SIMD + r=64 + tag walk
- Build on whichever Phase 22.X was best in each category

**Hour 5-8:** Final benchmark
- 1, 2, 4 threads on AutoLab-equivalent hardware
- Project to 80-bit: extrapolate from validated 20-40 bit measurements

**Deliverable:** `phase22_combined_result.md`.

**Acceptance criterion:** ≥27e6 ops/sec on 2-thread = 80-bit benchmark
feasible in <4 wall-hours.

## Day 13-14: AutoLab integration (16 hours)

**Day 13:** ctypes wrapper
- Wrap C functions with `ctypes`
- Python `solve.py` integration
- Verify against current Sage Python baseline at 30-bit

**Day 14:** End-to-end test
- Run against actual AutoLab verifier
- Confirm score improvement
- Submit as final v7 deliverable

**Deliverable:** Updated `solve.py` linked to `phase22_combined.so`.

**Acceptance criterion:** AutoLab benchmark score ≥7500 (significant
improvement over 6668.96 baseline).

## Risk register

| Risk | Mitigation | Day |
|------|------------|-----|
| Barrett bugs | Compare against GMP on 10^6 random inputs | 1 |
| SIMD platform diff | Build matrix for AVX-2, AVX-512, NEON | 3 |
| Cache miss on r=64 | Measure with `perf stat -e cache-misses` | 8 |
| HT mostly fictional | Just measure; document either way | 9 |
| Tag walk correctness | Verify secret recovery on toy ECDLP | 11 |
| C integration overhead | Profile `solve.py` after wiring | 13 |

## Total time

- **Day 0-14: ~14 working days** (~3 calendar weeks for one engineer)
- Could parallelize Phase 22.1 (Barrett) and Phase 22.3 (r=64)
  since they're independent
- Phase 22.2 (SIMD) is the highest risk; allocate buffer

## Success criteria

v7 ships if:
1. Combined-stack throughput on 2 CPUs ≥ 27e6 ops/sec
2. End-to-end Pollard rho recovers an 80-bit ECDLP secret in
   <4 wall-clock hours on AutoLab-equivalent hardware
3. AutoLab verifier shows score improvement when integrated

If v7 partially ships (e.g., 20e6 ops/sec, 6-hour wall time), document
which sub-phases delivered and which under-performed. Plan v7.5 to
close the remaining gap.

## What v7 will NOT cover

- Algorithmic novelty (deferred to v8)
- Multi-machine distributed attacks (deferred to v9 Phase 39)
- Custom silicon (out of any reasonable scope)
- Post-quantum considerations (orthogonal to v7)
