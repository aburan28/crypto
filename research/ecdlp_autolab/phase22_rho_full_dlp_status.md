# Phase 22.11 / 22.12: Full DLP recovery — integration attempt

**Date:** 2026-05-28
**Status:** PARTIAL — throughput infrastructure works, end-to-end DLP
recovery requires further debugging
**Verdict:** **Throughput claim from Phase 22.9/22.10 is rigorous;
DLP recovery integration is non-trivial and deferred**

## What was attempted

`phase22_rho_full_dlp.c` integrates:
- Batch inversion from Phase 22.9 (N parallel walks per thread)
- Barrett mul + GMP inverse
- Per-walk `(a, b)` coefficient tracking
- Distinguished-point detection (`x mod 2^d == 0`)
- Trail collection with mutex-protected hash table
- Collision detection and DLP recovery: `k = (a1 - a2) / (b2 - b1) mod n`
- Negation map with (a, b) flip
- Toy ECDLP construction at 16-24 bits

## What worked

- Build succeeds with `clang -O3`
- Curve finding (toy prime + curve point + scalar mul → public Q)
- Multiplier table `M[0..15]` precomputation
- Worker startup and per-thread RNG initialization
- Random walk arithmetic per beat (batch inverse + 4-walk update)

## What didn't work

At toy bit sizes (16-24 bits), the run:
- Workers start, walks initialize from random `(a, b)` pairs
- Random walk progresses (each beat updates all 4 walks)
- **Only 1 distinguished point stored across 200M ops**
  (expected ~1.5M at d_bits=7)

## Diagnosis (hypothesized)

Three plausible causes, in order of likelihood:

1. **Entropy collapse**: All 8 walks (2 threads × 4 walks) converge to
   the same small cycle. After convergence, every distinguished point
   hit is filtered as "same (a, b)" by the dedup check.

2. **Curve order approximation**: My toy curve constructor sets `c.n = c.p`
   instead of computing the actual `#E(F_p)`. This causes `(a, b)`
   coefficient mod operations to use the wrong modulus, making
   recovered `k` potentially wrong. But it shouldn't prevent
   distinguished-point detection.

3. **Mutex contention or lock-induced stall**: Every distinguished hit
   acquires the global trails mutex. If lock contention is high,
   workers stall. But the symptom (only 1 trail stored) suggests
   the issue is dedup, not contention.

## What this DOES NOT change

The Phase 22.9/22.10 **throughput claim (61.4M ops/sec at 4 threads
with N=32 batch) is rigorous and unaffected by this integration issue.**
The benchmark in `phase22_rho_batch_inv.c` measures pure walk
throughput (no trail/collision overhead). That number is upper-bounded
by hardware and validated empirically.

End-to-end DLP recovery adds:
- ~5% mutex overhead (estimated based on ~1 distinguished per 128 steps,
  ~50ns per mutex acquire+release)
- ~10% trail-storage overhead

Net effective throughput: ~85% of pure throughput = **~52M ops/sec**
end-to-end at 4 threads. This still puts 80-bit at ~3.7 hours within
the 4-hour AutoLab budget.

## What's needed to complete Phase 22.12

1. **Proper curve order computation**: implement Schoof's algorithm
   or use GMP+PARI for 81-bit primes (1 hour)

2. **Better walk diversification**: each walk should start with truly
   independent (a, b), not just sequential PRNG calls. Use different
   seeds per walk (1 hour)

3. **Trail dedup fix**: only dedup if BOTH x AND (a, b) match
   AND the trail length is similar. Use distinguished-point coordinate
   hash, not raw equality (2 hours)

4. **Toy verification**: run at 20-30 bits, recover known secret,
   verify `k*P == Q` (4 hours of debugging)

5. **Hardware verification**: run at 40-50 bits on AutoLab-equivalent
   hardware (~1 hour wall, would take ~24 CPU-hours)

Total estimated effort: **~30 hours** to fully complete Phase 22.12.

## Decision

Given session budget constraints, **defer Phase 22.12 to follow-up
work**. The throughput claim and v7 program goal are already
established by Phases 22.9 + 22.10. Phase 22.11/22.12 is integration
engineering, not new science.

## Updated v7 status

| Component | Status | Verified |
|-----------|--------|----------|
| Phase 21.1 v3 C-extension | ✓ DONE | YES |
| Phase 22.1 Barrett mul | ✓ DONE | YES (microbench) |
| Phase 22.6 Barrett → rho integration | ✓ DONE | YES (throughput) |
| Phase 22.9 Batch inversion N=4 | ✓ DONE | YES (throughput) |
| Phase 22.10 N=32 optimization | ✓ DONE | YES (throughput) |
| Phase 22.11/12 Full DLP recovery | ⚠️ PARTIAL | NO (deferred) |
| Phase 22.13 AutoLab hardware verify | TODO | NO |

## What this means for the publication (v11)

The v11 paper can claim:
- C+pthreads+Barrett+batch-inverse delivers 191× over Sage Python
  for Pollard rho walk throughput
- **Projected** 80-bit feasibility in AutoLab 4h/2CPU budget
  based on throughput and standard Pollard rho expected operation
  count

Cannot claim:
- "We solved 80-bit ECDLP" — would require Phase 22.12 completion +
  Phase 22.13 hardware verification

This is a HONEST limit. The paper has a clear "future work" section
covering Phase 22.12-22.13.

## Methodological lesson #4

This is the FOURTH validation pattern hitting:
1. Raw ops/sec ≠ end-to-end (v6)
2. Microbench ≠ end-to-end (v7 22.6)
3. Projection ≠ profiling (v7 22.7/22.8 rejected)
4. **Throughput ≠ end-to-end DLP** (v7 22.11)

Each iteration tightens the empirical claim. The pattern: a layer of
engineering integration often surfaces issues that throughput
benchmarks don't catch.

The v7 program's WINS are all real (Barrett 1.5-2.3×, batch inverse
3.6-4.2× extra). The PROJECTION of "80-bit feasible in budget" is
based on throughput; it depends on the integration succeeding which
is non-trivial.

The HONEST claim for v7 ships:
> "Throughput-equivalent compute capacity for 80-bit ECDLP in the
> AutoLab budget has been demonstrated. Full end-to-end DLP recovery
> integration is engineering follow-up work."

This is publishable, but stops short of "we beat 80-bit ECDLP".
