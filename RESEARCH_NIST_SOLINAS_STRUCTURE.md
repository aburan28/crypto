# Research Note: Cyclotomic Structure in NIST Solinas Primes

**Subject.** A reproducible structural probe of 29 scaled-down NIST-shape
Solinas primes (the calibration set documented in
[docs/ECDLP_CALIBRATION_PRIMES.md](docs/ECDLP_CALIBRATION_PRIMES.md))
that exposes one inherited property — algebraic super-smoothness of p±1
governed by cyclotomic factorisation — and confirms the absence of every
other classical structural weakness (MOV, anomalous, special CM, GLV-trivial,
twist).

The same algebraic identity, applied at full NIST sizes, shows that **P-224
has a p−1 that is 46-bit smooth** as a direct consequence of its form. This
is well known to specialists but not always foregrounded in deployment
discussions. The note records what's actually there.

---

## 1. Motivation

The calibration set was constructed for ECDLP attack work: every row is a
miniature NIST curve preserving the *exact sign-topology* of one of the
NIST primes (P-192, P-224, P-256, P-384, P-521 / Mersenne) at 48–128 bits,
fitted with a NIST-style `y² = x³ − 3x + b` curve of near-prime order.

The question we asked: *what cryptanalytic structure, if any, do these
primes inherit from the NIST topologies they imitate?*

The answer is sharper than expected — the smoothness pattern in p±1 is
not statistical, it's algebraic; it follows from the cyclotomic
decomposition of 2ᵐ − 1, and it applies identically to the NIST primes
themselves.

## 2. The cyclotomic mechanism

Write 2ᵐ − 1 as a product of cyclotomic values at 2:

$$
2^m - 1 \;=\; \prod_{d \mid m} \Phi_d(2)
$$

For two of the four NIST topologies, the Solinas form directly exposes such
a factor in p ± 1:

| Topology                       | NIST anchor                 | Factor exposed                                       |
|--------------------------------|-----------------------------|------------------------------------------------------|
| `2ⁿ − 2ᵏ − 1`                  | P-192 (n=192, k=64)         | **p + 1** = 2ᵏ · (2^(n−k) − 1)                       |
| `2ⁿ − 2ᵏ + 1`                  | P-224 (n=224, k=96)         | **p − 1** = 2ᵏ · (2^(n−k) − 1)                       |
| `2ⁿ − 1`                       | P-521                       | **p + 1** = 2ⁿ (fully 2-smooth)                       |
| `2ⁿ − 2ᵃ + 2ᵇ + 2ᶜ − 1` (5-term) | P-256                     | no clean factorisation                               |
| `2ⁿ − 2ᵃ − 2ᵇ + 2ᶜ − 1` (4-term) | P-384                     | no clean factorisation                               |

In every case where the factor is exposed, p ± 1 is forced to be smooth in
the precise sense that its largest prime factor cannot exceed
max_{d | (n−k)} (largest prime divisor of Φ_d(2)). Since the cyclotomic
values Φ_d(2) are small for small d, this is a strong bound.

## 3. P-224 worked out

For P-224 = 2²²⁴ − 2⁹⁶ + 1:

```
p - 1 = 2^96 · (2^128 - 1)
      = 2^96 · 3 · 5 · 17 · 257 · 641 · 65537 · 274177 · 6700417 · 67280421310721
```

Largest prime factor of p − 1 is **67 280 421 310 721 ≈ 2⁴⁶**. A random
224-bit prime would have its largest factor near 2²²⁴; here it is at most
2⁴⁶, which makes p − 1 fully 10¹⁴-smooth.

## 4. Reproducibility — the probe script

A single Sage script verifies the identity at NIST P-224 and across the
calibration scale-downs:

```python
# verify p - 1 = 2^k * (2^(n-k) - 1) for the P-224 family
cases = [(56,24), (64,24), (70,40), (80,48), (84,68),
         (96,32), (126,18), (224,96)]
for n, k in cases:
    p   = 2^n - 2^k + 1
    lhs = p - 1
    rhs = 2^k * (2^(n-k) - 1)
    assert lhs == rhs
    print(f"n={n:3d} k={k:3d}  factor(2^{n-k}-1) = {factor(2^(n-k) - 1)}")
```

Output (`n−k` column shows the cyclotomic exponent driving the factorisation):

| n   | k  | n−k | factor(2^(n−k) − 1)                                              | largest prime |
|-----|----|-----|------------------------------------------------------------------|---------------|
| 56  | 24 | 32  | `3·5·17·257·65537`                                               | 65 537        |
| 64  | 24 | 40  | `3·5²·11·17·31·41·61681`                                         | 61 681        |
| 70  | 40 | 30  | `3²·7·11·31·151·331`                                             | 331           |
| 80  | 48 | 32  | `3·5·17·257·65537`                                               | 65 537        |
| 84  | 68 | 16  | `3·5·17·257`                                                     | 257           |
| 96  | 32 | 64  | `3·5·17·257·641·65537·6700417`                                   | 6 700 417     |
| 126 | 18 | 108 | `3⁴·5·7·13·19·37·73·109·87211·246241·262657·279073`              | 279 073       |
| **224** | **96** | **128** | **`3·5·17·257·641·65537·274177·6700417·67280421310721`** | **≈ 2⁴⁶** |

The full structural probe (over all 29 (prime, curve) pairs) lives at
`/tmp/probe_structure2.sage` while the work is in flight; the salient
parts of its output are reproduced in section 6 below.

## 5. What this means for cryptanalysis

This is the most important section, because the smoothness is striking
enough to invite over-interpretation.

**Pure ECDLP is unaffected.** ECDLP lives in the additive group E(F_p),
not the multiplicative group F_p\*. The factorisation of p − 1 controls
the multiplicative group only. Pollard ρ on E(F_p) depends on the order
of the chosen subgroup, not on p − 1.

**Companion DLP in F_p\* is broken.** If the same prime were reused as the
modulus of a Schnorr / DSA / ElGamal scheme in the multiplicative group,
Pohlig-Hellman would solve the DLP almost instantly given the smoothness.
This is why deployment guides explicitly say to use distinct primes for
ECC and finite-field DLP — but the *structural reason* P-224 cannot
double as an F_p\* modulus is exactly the cyclotomic identity above.

**Pairing-based crypto is doubly unsafe.** Pairings reduce ECDLP to DLP in
F_{p^k}\*, where k is the embedding degree. SNFS variants for F_{p^k}\*
exploit smoothness in p ± 1. For our calibration set k > 200 everywhere
and pairings are computationally infeasible regardless — so this is a
hypothetical only — but for any future curve over a P-224-shape prime
that did admit small k, the smoothness would be doubly problematic.

**Pollard p − 1 factoring is irrelevant.** These primes are prime; p − 1
factoring solves the wrong problem.

**Twist attacks are independent.** Twist security is governed by #E_twist
= p + 1 + t, not by p ± 1 directly. The structural smoothness has no
bearing here.

So the practical posture is: **the smoothness is benign for ECDLP, lethal
for any companion F_p\* protocol, and a clean fingerprint identifying the
P-192 / P-224 / P-521 topologies in the wild.**

## 6. Structural probe summary (29 calibration curves)

Full probe output is in `/tmp/structure_report.txt` during the session;
the salient cross-row signal:

```
=== Universal fingerprints (all 29 rows) ===
  k_MOV > 200    everywhere    →  pairings unreachable; no Frey-Rück route
  #E ≠ p         everywhere    →  no anomalous curve; no Smart break
  j ∉ {0, 1728}  everywhere    →  no built-in GLV (forced by a = -3)
  CM disc ≈ -4p  everywhere    →  no special CM construction
  h_E ∈ {1..8}                 →  all near-prime as designed

=== Family-correlated smoothness (out of n bits of rough part) ===
  P-192/72 :  p+1 rough =  1b /  72   ←  algebraic
  P-224/56 :  p-1 rough =  1b /  56   ←  algebraic
  P-224/64 :  p-1 rough =  1b /  64
  P-224/70 :  p-1 rough =  1b /  70
  P-224/80 :  p-1 rough =  1b /  80
  P-224/84 :  p-1 rough =  1b /  84
  P-521/61 :  BOTH p-1 and p+1 rough = 1b /  61   ←  M61 doubly smooth
  P-521/89 :  p+1 rough =  1b /  89   ←  Mersenne
  P-521/107:  p+1 rough =  1b / 107   ←  Mersenne
  P-521/127:  p+1 rough =  1b / 127   ←  Mersenne (M127)

  P-256/*  :  rough parts 25–113b / n   ←  no algebraic forcing
  P-384/*  :  rough parts 67–116b / n   ←  no algebraic forcing

=== Twist cofactor outliers (h_T below 1000) ===
  P-224/56     h_T =     33   ←  near-twist-secure by accident
  P-256/96     h_T =     20   ←  near-twist-secure by accident
  P-256/112    h_T =    246
  P-521/107    h_T =   2983
```

## 7. Open questions / follow-ups

1. **Quantify Pohlig-Hellman cost on F_p\*** for each P-224-family row. The
   expected runtime is essentially `√(largest prime factor of p − 1)`
   group ops — so e.g. ~2²³ for P-224/56 and ~2²³ for full P-224 itself.
   A runnable demo would make the dual-DLP weakness concrete.
2. **Targeted twist attack** against P-256/96 (h_T = 20) and P-224/56
   (h_T = 33) — verify the small-subgroup-confinement leakage is exactly
   ⌈log₂ h_T⌉ bits.
3. **GLV-friendly variant**: drop the `a = −3` constraint and refit
   each calibration row with `a = 0` to land on j = 0 curves with the
   built-in ζ₃-endomorphism. Lets you separately calibrate
   GLV-accelerated ρ.
4. **Does any 5-term P-256-shape prime admit a hidden cyclotomic
   factorisation we missed?** The probe found none, but the search
   space is large; an exhaustive check of (a, b, c) triples would
   close the question.
5. **Lift to the actual NIST primes**: re-run the same probe on P-192,
   P-224, P-256, P-384, P-521 themselves and publish the full
   factorisations of p − 1 and p + 1, with attribution to the Solinas
   form that forces them.

---

**Cross-references.**

- [docs/ECDLP_CALIBRATION_PRIMES.md](docs/ECDLP_CALIBRATION_PRIMES.md) —
  the calibration set itself plus a condensed version of these findings.
- [docs/ECDLP_ATTACK_MATRIX.md](docs/ECDLP_ATTACK_MATRIX.md) — taxonomy
  of attacks the calibration set is meant to exercise.
- [src/cryptanalysis/lattice.rs](src/cryptanalysis/lattice.rs),
  [src/cryptanalysis/summation_poly/](src/cryptanalysis/summation_poly/),
  [src/cryptanalysis/nist_p224/](src/cryptanalysis/nist_p224/) — attack
  implementations these calibration curves are meant to validate.
- [RESEARCH_SECP256K1_CM.md](RESEARCH_SECP256K1_CM.md) — companion CM
  audit on secp256k1; the smoothness analysis there pairs naturally
  with the P-192-family observation above (secp256k1 is the 2ⁿ − 2³² − c
  variant and has the same p + 1 smoothness pattern).
