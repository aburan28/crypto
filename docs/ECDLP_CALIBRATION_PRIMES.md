# ECDLP Calibration Primes — NIST-Shape Scale-Down Set

A set of 29 prime fields and matching `y² = x³ − 3x + b` curves designed for
**end-to-end ECDLP attack calibration** against the NIST curves. Each row
preserves the *exact sign-topology* of one of the NIST primes (P-192, P-224,
P-256, P-384, P-521 / Mersenne) while scaling the bit-size down into ranges
that are reachable on commodity hardware, a small cluster, or — at the
ceiling — Certicom-challenge scale.

The intent is that an attack pipeline (Pollard ρ, GHS / Weil descent,
summation-polynomial, lattice / LLL-GS) can be validated on the
sanity-tier rows, scaled on the mid-tier rows, and benchmarked at the
top tier without ever having to change the curve *topology* — only the
bit-size — when finally pointed at the NIST primes themselves.

All curves use the NIST coefficient convention `a = −3`. The `b` value
listed is the smallest positive `b` that yields a near-prime group order;
the cofactor `h` is the small factor and `n′` is the prime subgroup order.

---

## P-192 family — form `2ⁿ − 2ᵏ − 1`

NIST anchor: **P-192** = `2¹⁹² − 2⁶⁴ − 1` (k/n = 1/3).

| n   | prime                  | b  | h | log₂(n′) | log₂(ρ cost) |
|-----|------------------------|----|---|----------|--------------|
| 48  | `2⁴⁸ − 2¹⁶ − 1`       | 4  | 2 | 47.00    | 23.3         |
| 56  | `2⁵⁶ − 2¹⁴ − 1`       | 5  | 2 | 55.00    | 27.3         |
| 64  | `2⁶⁴ − 2²⁹ − 1`       | 78 | 1 | 64.00    | 31.8         |
| 72  | `2⁷² − 2²⁴ − 1`       | 65 | 1 | 72.00    | 35.8         |
| 80  | `2⁸⁰ − 2²⁴ − 1`       | 57 | 1 | 80.00    | 39.8         |
| 96  | `2⁹⁶ − 2³⁷ − 1`       | 29 | 1 | 96.00    | 47.8         |
| 112 | `2¹¹² − 2³⁵ − 1`      | 32 | 1 | 112.00   | 55.8         |
| 120 | `2¹²⁰ − 2³¹ − 1`      | 21 | 2 | 119.00   | 59.3         |

## P-224 family — form `2ⁿ − 2ᵏ + 1`

NIST anchor: **P-224** = `2²²⁴ − 2⁹⁶ + 1` (k/n = 3/7).

| n   | prime                            | b  | h | log₂(n′) | log₂(ρ cost) |
|-----|----------------------------------|----|---|----------|--------------|
| 56  | `2⁵⁶ − 2²⁴ + 1`                 | 87 | 1 | 56.00    | 27.8         |
| 64  | `2⁶⁴ − 2²⁴ + 1`                 | 14 | 3 | 62.42    | 31.0         |
| 70  | `2⁷⁰ − 2⁴⁰ + 1`                 | 7  | 3 | 68.42    | 34.0         |
| 80  | `2⁸⁰ − 2⁴⁸ + 1`                 | 3  | 3 | 78.42    | 39.0         |
| 84  | `2⁸⁴ − 2⁶⁸ + 1`                 | 91 | 3 | 82.42    | 41.0         |
| 96  | `2⁹⁶ − 2³² + 1`                 | 43 | 1 | 96.00    | 47.8         |
| 126 | `2¹²⁶ − 2¹⁸ + 1`                | 26 | 4 | 124.00   | 61.8         |

## P-256 family — form `2ⁿ − 2ᵃ + 2ᵇ + 2ᶜ − 1`

NIST anchor: **P-256** = `2²⁵⁶ − 2²²⁴ + 2¹⁹² + 2⁹⁶ − 1` (exponent ratios 7/8, 3/4, 3/8).

| n   | prime                                          | b  | h | log₂(n′) | log₂(ρ cost) |
|-----|------------------------------------------------|----|---|----------|--------------|
| 64  | `2⁶⁴ − 2⁵⁶ + 2⁴⁷ + 2²⁴ − 1`                  | 35 | 2 | 62.99    | 31.3         |
| 80  | `2⁸⁰ − 2⁷⁰ + 2⁵⁹ + 2³⁰ − 1`                  | 8  | 3 | 78.41    | 39.0         |
| 96  | `2⁹⁶ − 2⁸⁴ + 2⁷³ + 2³⁶ − 1`                  | 7  | 4 | 94.00    | 46.8         |
| 112 | `2¹¹² − 2⁹⁹ + 2⁸⁵ + 2⁴³ − 1`                 | 74 | 6 | 109.41   | 54.5         |
| 128 | `2¹²⁸ − 2¹¹³ + 2⁹⁵ + 2⁴⁸ − 1`                | 13 | 8 | 125.00   | 62.3         |

## P-384 family — form `2ⁿ − 2ᵃ − 2ᵇ + 2ᶜ − 1`

NIST anchor: **P-384** = `2³⁸⁴ − 2¹²⁸ − 2⁹⁶ + 2³² − 1` (exponent ratios 1/3, 1/4, 1/12).

| n   | prime                                          | b  | h | log₂(n′) | log₂(ρ cost) |
|-----|------------------------------------------------|----|---|----------|--------------|
| 80  | `2⁸⁰ − 2²⁷ − 2²² + 2⁷ − 1`                   | 20 | 4 | 78.00    | 38.8         |
| 96  | `2⁹⁶ − 2³⁴ − 2²³ + 2⁸ − 1`                   | 88 | 3 | 94.42    | 47.0         |
| 112 | `2¹¹² − 2³⁷ − 2²⁷ + 2⁹ − 1`                  | 41 | 3 | 110.42   | 55.0         |
| 120 | `2¹²⁰ − 2⁴¹ − 2³⁰ + 2¹⁰ − 1`                 | 10 | 6 | 117.42   | 58.5         |
| 128 | `2¹²⁸ − 2⁴¹ − 2³² + 2¹³ − 1`                 | 17 | 8 | 125.00   | 62.3         |

## P-521 family — pure Mersenne `2ⁿ − 1`

NIST anchor: **P-521** = `2⁵²¹ − 1`. Scale-downs are the genuine Mersenne primes
M₆₁, M₈₉, M₁₀₇, M₁₂₇ (rows 9, 11, 12, and 13 in Wikipedia's Mersenne table).

| n   | prime                                          | b   | h | log₂(n′) | log₂(ρ cost) |
|-----|------------------------------------------------|-----|---|----------|--------------|
| 61  | `2⁶¹ − 1` (M₆₁)                               | 102 | 4 | 59.00    | 29.3         |
| 89  | `2⁸⁹ − 1` (M₈₉)                               | 10  | 4 | 87.00    | 43.3         |
| 107 | `2¹⁰⁷ − 1` (M₁₀₇)                             | 114 | 5 | 104.68   | 52.2         |
| 127 | `2¹²⁷ − 1` (M₁₂₇)                             | 22  | 3 | 125.42   | 62.5         |

---

## Prime hex values

```
# P-192 shape
P-192/48   p = 0xfffffffeffff
P-192/56   p = 0xffffffffffbfff
P-192/64   p = 0xffffffffdfffffff
P-192/72   p = 0xfffffffffffeffffff
P-192/80   p = 0xfffffffffffffeffffff
P-192/96   p = 0xffffffffffffffdfffffffff
P-192/112  p = 0xfffffffffffffffffff7ffffffff
P-192/120  p = 0xffffffffffffffffffffff7fffffff

# P-224 shape
P-224/56   p = 0xffffffff000001
P-224/64   p = 0xffffffffff000001
P-224/70   p = 0x3fffffff0000000001
P-224/80   p = 0xffffffff000000000001
P-224/84   p = 0xffff00000000000000001
P-224/96   p = 0xffffffffffffffff00000001
P-224/126  p = 0x3ffffffffffffffffffffffffffc0001

# P-256 shape
P-256/64   p = 0xff00800000ffffff
P-256/80   p = 0xffc0080000003fffffff
P-256/96   p = 0xfff002000000000fffffffff
P-256/112  p = 0xfff800200000000007ffffffffff
P-256/128  p = 0xfffe0000800000000000ffffffffffff

# P-384 shape
P-384/80   p = 0xfffffffffffff7c0007f
P-384/96   p = 0xfffffffffffffffbff8000ff
P-384/112  p = 0xffffffffffffffffffdff80001ff
P-384/120  p = 0xfffffffffffffffffffdffc00003ff
P-384/128  p = 0xfffffffffffffffffffffdff00001fff

# P-521 shape (Mersenne)
P-521/61   p = 0x1fffffffffffffff
P-521/89   p = 0x1ffffffffffffffffffffff
P-521/107  p = 0x7ffffffffffffffffffffffffff
P-521/127  p = 0x7fffffffffffffffffffffffffffffff
```

---

## Suggested experiment progression

- **Sanity tier (≤ 64 bits)** — P-192/48, P-224/56, P-256/64, P-521/61.
  Full Pollard-ρ ECDLP runs in seconds; validate the ρ harness end-to-end
  across all four NIST topologies before scaling.
- **Single-machine breakable (72–96 bits)** — P-192/72, P-224/96, P-256/80,
  P-384/80, P-521/89. Solvable in hours on one box; calibrate
  distinguished-point parameters and replay-resistance.
- **Cluster tier (107–120 bits)** — M₁₀₇, P-192/112, P-256/112, P-384/112.
  Certicom-challenge scale; verify cluster scaling and the runtime curves
  of the structural attacks (GHS, summation-polynomial).
- **Aspirational (126–128 bits)** — P-192/120, P-224/126, P-256/128,
  P-384/128, M₁₂₇. Beyond casual hardware but exercises every optimization
  needed before attempting P-224 itself.

Each row is structurally a "miniature NIST curve": the same Solinas reduction
topology, the same `a = −3` coefficient choice, and an essentially-prime group
order. Attack code targeting the NIST primes should be drop-in usable on any
calibration row just by swapping `(p, b)`.

## Generation methodology

- Prime search: for each NIST topology, enumerate candidates by scaling the
  canonical exponent ratios to the target bit-size, then drift the inner
  exponents within a small window (≤ n/12) to find the closest prime to
  the true NIST geometry. Drift % is bounded — most rows are within 3 %.
- Curve fitting: for each prime, find the smallest positive `b ∈ [1, 200]`
  such that `#E(F_p)` for `y² = x³ − 3x + b` has cofactor `h ≤ 4`
  (relaxed to `h ≤ 8` at ≥ 112 bits). Order computed via Sage's SEA.
- All primes verified prime by sympy `isprime` (Baillie-PSW + Miller-Rabin).
- Pollard-ρ cost is `√(π · n′ / 4)`; the log₂ figure is `½ log₂(n′) − 0.17`.

Cross-reference: [ECDLP_ATTACK_MATRIX.md](ECDLP_ATTACK_MATRIX.md) for the
attack taxonomy these calibration curves are meant to exercise.

---

# Cryptanalytic structural fingerprint of the calibration set

A structural probe of all 29 (prime, curve) pairs revealed five findings worth
recording. The headline: the Solinas reduction shape that makes the NIST primes
fast also *algebraically forces* their p±1 to be super-smooth — this is benign
for pure ECDLP but is the load-bearing reason these primes should not be reused
as the modulus of any companion multiplicative-group DLP.

## 1. P-224-family primes have algebraically forced super-smooth p−1

For any prime of the form `p = 2ⁿ − 2ᵏ + 1` (the P-224 topology),

```
    p − 1  =  2ⁿ − 2ᵏ  =  2ᵏ · (2^(n−k) − 1)
```

so `p − 1` inherits *every* factor of the cyclotomic value 2^(n−k) − 1 =
∏_{d | (n−k)} Φ_d(2). These factors are systematically small. Direct
factorisations across the P-224 row of the calibration set:

| n    | k   | n−k | factor(2^(n−k) − 1)                                        | largest prime |
|------|-----|-----|-------------------------------------------------------------|---------------|
| 56   | 24  | 32  | `3·5·17·257·65537`                                          | 65 537        |
| 64   | 24  | 40  | `3·5²·11·17·31·41·61681`                                    | 61 681        |
| 70   | 40  | 30  | `3²·7·11·31·151·331`                                        | 331           |
| 80   | 48  | 32  | `3·5·17·257·65537`                                          | 65 537        |
| 84   | 68  | 16  | `3·5·17·257`                                                | 257           |
| 96   | 32  | 64  | `3·5·17·257·641·65537·6700417`                              | 6 700 417     |
| 126  | 18  | 108 | `3⁴·5·7·13·19·37·73·109·87211·246241·262657·279073`         | 279 073       |
| **224**  | **96**  | **128** | **`3·5·17·257·641·65537·274177·6700417·67280421310721`** | **≈ 2⁴⁶** |

**NIST P-224 itself has a p − 1 that is 46-bit smooth.** This is intrinsic to
the form, not a property of any particular bit-size choice.

## 2. Mirror pattern across the other NIST topologies

The same identity, applied to the other Solinas shapes:

| Family form               | Smooth side          | Mechanism                            |
|---------------------------|----------------------|--------------------------------------|
| `2ⁿ − 2ᵏ − 1`             | **p + 1**            | p+1 = 2ⁿ − 2ᵏ = 2ᵏ·(2^(n−k) − 1)      |
| `2ⁿ − 2ᵏ + 1`             | **p − 1**            | (above)                              |
| `2ⁿ − 1` (Mersenne)       | **p + 1 = 2ⁿ**       | fully 2-smooth                       |
| `2ⁿ − 2ᵃ + 2ᵇ + 2ᶜ − 1` (P-256) | neither       | no clean cyclotomic factorisation    |
| `2ⁿ − 2ᵃ − 2ᵇ + 2ᶜ − 1` (P-384) | neither       | no clean cyclotomic factorisation    |

Observed in the probe (`p±1 rough bit-length / n`):

- P-192/72: `p+1 rough = 1b / 72` (p+1 = 2⁷² − 2²⁴, cyclotomic)
- P-224/56–84: `p−1 rough = 1b / 56–84` (all super-smooth)
- P-521/61, 89, 107, 127: `p+1 rough = 1b / n` (fully 2-smooth)
- P-521/61: **both** `p−1` and `p+1` are 1-bit rough (M₆₁ doubly smooth)
- P-256/* and P-384/*: rough parts are 49–116 bits — *no* special structure

## 3. Universally good ECDLP-relevant fingerprints

The probe also confirmed every calibration curve is structurally sound on the
fingerprints that matter directly for ECDLP:

- **Embedding degree k_MOV > 200** for all 29 — MOV / Frey-Rück reduction to
  DLP in F_{p^k}* is computationally infeasible. Pairings on these curves
  are likewise unreachable, which is a security feature here.
- **No anomalous curve** (#E ≠ p in every case) — Smart / Semaev /
  Satoh–Araki additive-group break does not apply.
- **j-invariant ∉ {0, 1728}** for any row — the `a = −3` constraint forbids
  j = 1728, and no row landed at j = 0 either, so no built-in GLV
  endomorphism comes for free.
- **CM discriminants are generically large** (∼ −4p) — these are not
  special-CM curves.

## 4. Two accidentally-near-twist-secure curves

The `b` values were optimised purely for main-curve cofactor; the twist was
ignored. Most twist orders are wildly composite as a result, but two outliers
fell out of the search anyway:

| curve         | b   | h_E | **twist cofactor h_T**       |
|---------------|-----|-----|-------------------------------|
| **P-224/56**  | 87  | 1   | **33** (= 3·11)               |
| **P-256/96**  | 7   | 4   | **20** (= 2²·5)               |

These are practically twist-secure (≤ 5 bits of small-subgroup leakage) and
make natural positive controls when validating invalid-curve attacks against
the twist-vulnerable rows.

Mid-range twist cofactors worth noting:

| curve       | h_T     |
|-------------|---------|
| P-521/107   | 2 983   |
| P-384/112   | 60 337  |
| P-384/96    | 63 259  |
| P-256/112   | 246     |
| P-256/128   | 199 144 |

All other rows have h_T ≥ 10⁷ — strongly twist-vulnerable, well-suited as
targets for invalid-curve research.

## 5. Properly framing the cryptanalytic implications

The smoothness pattern matters for *adjacent* attack surfaces — not for ECDLP
itself:

| Attack                                | Affected by p±1 smoothness? | Relevance to these curves              |
|---------------------------------------|------------------------------|----------------------------------------|
| Pollard ρ on E(F_p)                   | No                           | Depends only on \|n′\|                  |
| Pollard p−1 / Williams p+1            | Only for factoring composites| These primes *are* prime; N/A          |
| MOV / Frey-Rück → DLP in F_{p^k}*     | Indirectly via SNFS          | k > 200, so unreachable                |
| **SNFS in F_p\* if reused for Schnorr/DSA** | **Yes**                | **P-192 / P-224 / P-521 family unsafe**|
| Twist (invalid-curve)                 | Only via twist cofactor      | See section 4                          |
| GHS / Weil descent                    | No (needs subfield structure)| Prime fields — does not apply          |
| Summation-polynomial / Diem index calc | No                          | Curve-structure dependent, not p±1     |

The bottom line: **the calibration primes are weak field choices for any
companion multiplicative-group DLP — and that is the same weakness P-224
itself inherits — but for pure ECDLP they are as hard as their bit-size
suggests.**

## Suggested follow-up experiments

- Run Pohlig–Hellman on **F_p\*** for the P-224-family primes — should be fast
  at scale-down sizes precisely because p−1 is so smooth. Direct evidence
  of the "dual-DLP weakness".
- Targeted twist attack against P-256/96 and P-224/56 — verify the
  small-subgroup-confinement leakage is exactly ⌈log₂ h_T⌉ bits.
- GLV-friendly variant: drop `a = −3` and refit with `a = 0` to get
  j = 0 curves with built-in ζ₃-endomorphism. Lets you separately
  calibrate GLV-accelerated ρ.
- Index-calculus stress test using the highest-smoothness P-224 rows
  against the [summation-poly](../src/cryptanalysis/summation_poly/)
  implementation. The expected result is *no* speedup from the inherited
  field structure; confirming this is itself valuable.
