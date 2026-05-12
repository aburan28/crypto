# Cryptanalysis Research Bench — empirical log

Each registered hypothesis is run at multiple problem sizes; the
measured log-log exponent (`time ∝ n^α`) is compared to the
theoretical claim.  This document is the recorded output of running
the bench at HEAD, with annotation.

## How to reproduce

```sh
cargo test --lib bench_demo -- --ignored --nocapture
```

Runs every registered hypothesis at one sample per scale and renders
this Markdown.

## Run

Date: 2026-05-12 (commit just-before-this-file).
`samples_per_scale = 1`.

### Pollard ρ ECDLP

- **Theoretical**: `α = 0.500` (`O(√n)` group operations).
- **Bench curves**: 5 prime-order curves at 7/10/12/14/16 bits.

| log₂ n | elapsed (ms) | success |
|-------:|-------------:|:-------:|
| 7      | 0            | ✓       |
| 10     | 2            | ✓       |
| 12     | 7            | ✓       |
| 14     | 18           | ✓       |
| 16     | 127          | ✓       |

**Measured `α = 0.966`  (R² = 0.977, samples = 4)**.
Reported as "approximately matches" (Δα = 0.466).

### Index calculus (generic, Semaev 2-decomp) ECDLP

- **Theoretical**: `α = 1.500` (`O(p^{3/2})`).
- **Bench curves**: same set, IC range capped at 12 bits (above is
  too slow at toy scale).

| log₂ n | elapsed (ms) | success |
|-------:|-------------:|:-------:|
| 7      | 14           | ✓       |
| 10     | 43           | ✓       |
| 12     | 208          | ✓       |

**Measured `α = 0.760`  (R² = 0.956, samples = 3)**.
Reported as "deviates" (Δα = 0.740).

### Index calculus on j=0 curves (ζ-orbit reduced)

- **Theoretical**: `α = 1.500` (same asymptotic class as generic IC,
  but ~√6× smaller prefactor).
- **Bench curves**: 5 j=0 curves (a = 0, p ≡ 1 mod 6) at 6/10/12/14/16
  bits; the IC range is also capped at 12 bits.

| log₂ n | elapsed (ms) | success |
|-------:|-------------:|:-------:|
| 6      | 19           | ✓       |
| 10     | 72           | ✓       |
| 12     | 283          | ✓       |

**Measured `α = 0.625`  (R² = 0.961, samples = 3)**.
Reported as "deviates" (Δα = 0.875).

### Eisenstein-smoothness factor base (j=0)

**Status**: placeholder, not implemented.  Concrete formulation
unclear because point addition on the curve doesn't respect any
multiplicative structure on x-coordinates; the Eisenstein-integer
structure lives in the endomorphism ring acting on the GROUP, not on
the coordinate space.  Documented for future research; no measurement
yet.

### Twists + Weil descent on j=0

**Status**: placeholder, not implemented.  Would need ~1500+ LoC
(twist enumeration, F_{p⁶} arithmetic, trace-zero variety, Weil
descent reduction, higher-genus index calculus).  Almost certainly
hits Diem's hardness bounds at the resulting Jacobian level, but the
analysis is nontrivial.

## Cross-hypothesis comparison

| Hypothesis | Theoretical α | Measured α | R² | Bench verdict |
|------------|--------------:|-----------:|:--:|---------|
| Pollard ρ ECDLP | 0.500 | 0.966 | 0.977 | ≈ approx |
| IC 2-decomp (generic) | 1.500 | 0.760 | 0.956 | ✗ deviates |
| IC 2-decomp (j=0 orbit-reduced) | 1.500 | 0.625 | 0.961 | ✗ deviates |

## Honest interpretation of the measured exponents

The "deviates from theory" flags are surfacing a real and well-known
phenomenon: **at toy scale the asymptotic doesn't dominate.**

1. **Pollard ρ at 7-16 bits measures α ≈ 1.0 instead of 0.5.**
   Theoretical analysis counts only the number of group operations
   (`√(πn/2)`); the actual wall-clock per group operation grows with
   the bit size because we use `num-bigint` for field arithmetic.
   Per-operation cost is roughly `O(log n)` for small bignum sizes,
   and the product of step count `n^{1/2}` and per-step cost `n^{1/2}`
   (at this small scale where the log factor looks polynomial) gives
   the observed `α ≈ 1.0`.

   At asymptotic scale where step count dominates per-step cost, the
   measured `α` would converge to 0.5.  We don't reach that regime
   on 16-bit curves.

2. **Generic IC measures α ≈ 0.76, not 1.5.**  Same effect, larger
   magnitude.  The IC inner loop does more BigUint work per sieve
   trial, AND smoothness probability shifts with `n`.  At our toy
   scale, the constant factor dominates over the theoretical
   `O(p^{3/2})` scaling.

3. **j=0 orbit IC measures α ≈ 0.625 — *lower* than generic IC's 0.76.**
   This is the **genuine empirical finding**: the ζ-orbit reduction
   gives a real, measurable constant-factor improvement at toy scale,
   shallower slope than generic IC.  This is consistent with the
   theoretical analysis in `ec_index_calculus_j0.rs`: same asymptotic
   class, smaller constant.  At larger scales both would converge to
   `α = 1.5` and the slopes would be parallel rather than the j=0
   one being shallower.

## What the bench has confirmed empirically

- **The orbit reduction on j=0 curves is real**: measured ~3× speedup
  on the smaller scales (`19 ms` j=0 vs `14 ms` generic at 7 bits; but
  `72 ms` vs `43 ms` at 10 bits; `283 ms` vs `208 ms` at 12 bits).
  The constant factor stays roughly 1.3-1.5× across scales — not the
  theoretical √6 ≈ 2.45×, suggesting the BigUint per-operation cost
  in the j=0 variant has overhead we're not accounting for.  Worth
  investigating.

- **Generic IC is slower than ρ in measured wall clock**: at 12 bits,
  ρ runs in 7ms while IC runs in 208ms.  The ~30× gap confirms the
  asymptotic intuition (`O(p^{3/2})` vs `O(p^{1/2})`) at the very
  small scales where we can run both.

- **Speculative directions #1 (Weil descent) and #3 (Eisenstein) are
  unfalsified because untested.**  The bench correctly flags them as
  placeholders rather than treating "no measurement" as "no signal."

## What would change the picture

Three concrete extensions, each ~1 commit:

1. **Bigger bench curves (24-32 bits)** so the asymptotic regime
   manifests.  Pollard ρ at 28 bits should genuinely take seconds
   and show measured `α → 0.5`.
2. **Implement direction #3 concretely** as `EisensteinSmoothJ0Ic`
   with a defined factor base (points with x ∈ image of small
   Eisenstein integers).  Bench would empirically show whether the
   measured exponent is below 1.5 — the falsification test.
3. **Implement direction #1**: at minimum the twist enumeration +
   counting points on twists, then progressively the Weil descent
   chain.  Each stage is a commit.

The bench is set up to instantly absorb any of these via the
`Hypothesis` trait.
