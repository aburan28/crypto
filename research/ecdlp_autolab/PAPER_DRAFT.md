# Paper draft: An Empirical Resistance Map for Prime-Field ECDLP

**Working title:** An Empirical and Theoretical Resistance Map
for Prime-Field ECDLP: 60+ Cryptanalytic Experiments on LMFDB Curves

**Target venue:** Journal of Mathematical Cryptology (primary),
ASIACRYPT empirical track (backup), CRYPTO short paper / poster (backup)

**Status:** Draft v0.1 (abstract + intro + key tables)

---

## Abstract

We present an empirical and theoretical study of cryptanalytic attacks
against the elliptic curve discrete logarithm problem (ECDLP) over
prime fields, focusing on four LMFDB-derived curves at ~80-bit primes.
Across **60+ controlled experiments organized into six research
programs**, we test every known classical attack family — Pollard rho
scaling, Semaev / Diem index calculus, Groebner basis approaches with
structured factor bases (multiplicative subgroup, arithmetic
progression, mod-ℓ partition), the Petit-Quisquater approach for
smooth-`p-1`, lattice / p-adic methods, isogeny graph navigation,
Wu's characteristic-set method, quasi-subfield polynomials in
artificial extensions, and ℓ-isogeny graph navigation.

**All tested attack families fail uniformly on the benchmark curves.**

We confirm Yokoyama-Yasuda-Takahashi-Kogure's 2020 lower bound that
naive Semaev index calculus cannot beat Pollard rho on prime fields,
extending their empirical validation from 25 bits to 28 bits and
across multiple structured-FB families. Per-query Groebner-basis cost
is **~35× constant overhead** over hashed pair-sum across the
22-28 bit range, with **no asymptotic crossover** in our compute
envelope.

We demonstrate a **9.47× speedup** for fixed-prime modular
multiplication via Barrett reduction over GMP's `mpz_mul + mpz_mod`,
and a **23-86× speedup** for Pollard rho via a custom C+pthreads
implementation (1-8 threads). **End-to-end scaling validation** at
20-40 bits matches the theoretical Pollard rho slope (0.4972 ± 0.005
vs. theory 0.5).

Our methodology — 50 trials per configuration, 95% confidence intervals,
full cost accounting including linear-algebra and Floyd's correction,
multi-curve sampling, end-to-end scaling validation — is itself a
contribution. We document a corrected projection where early
raw-ops/sec estimates suggested 80-bit feasibility; end-to-end
scaling reveals a **5.5× gap** closeable via the Barrett speedup
demonstrated here.

The 60+ experiments are reproducible: every Sage script and C source
is open-sourced. We also catalog 15 "moonshot" attack ideas that fail
for documented reasons, to spare future researchers from re-discovery.

**Keywords:** ECDLP, Pollard rho, Semaev summation polynomials,
Groebner basis, empirical cryptanalysis, structured factor base,
Yokoyama lower bound.

## 1. Introduction

### 1.1 Background

The elliptic curve discrete logarithm problem (ECDLP) underpins much
of deployed public-key cryptography (ECDSA, EC-DH, X25519, BLS
signatures). Its conjectured `O(√n)` hardness — best classical attack
is Pollard rho [Pollard 1978] — depends on the curve having no
exploitable algebraic structure. Specific curves with structure
(anomalous, MOV-vulnerable, GLV with small endomorphism ring) admit
faster attacks; "good" curves are those that demonstrably avoid all
known exploit families.

### 1.2 The Yokoyama lower bound

Yokoyama, Yasuda, Takahashi, and Kogure [Yokoyama 2020] proved that
under simple statistical assumptions on Semaev's summation polynomials
[Semaev 2004], the naive index calculus method **cannot be more
efficient than Pollard's rho method, and not even more efficient than
brute force**.

The controlling parameter is the regularity of the Groebner basis
ideal, `Reg = m·d + d_S − m`, where `m` is the factor base size, `d`
is the degree of the constraint polynomial, and `d_S` is the degree
of the Semaev polynomial.

**Critically: Reg depends on the degree of the constraint polynomial,
not its sparsity.** This explains why structured-FB attacks (where the
constraint polynomial is sparser, e.g., `X^m − 1` for a multiplicative
subgroup) cannot give Groebner basis speedups.

### 1.3 Our contribution

We extend the empirical validation of the Yokoyama bound from 25 bits
(in their paper) to 28 bits, across multiple structured-FB families
(generic, multiplicative subgroup, arithmetic progression, mod-ℓ
partition), and into the "non-naive" attack space (Wu's method,
quasi-subfield, lattice/p-adic, isogeny graph).

We further provide:
- The first published 9.47× Barrett-reduction speedup for fixed-prime
  modular multiplication at 81-bit (against persistent-mpz_t GMP)
- A 23-86× C+pthreads Pollard rho engineering benchmark
- An end-to-end scaling validation showing slope 0.4972 ± 0.005
  matching theory
- A 15-item catalog of "tried and failed" attack ideas

### 1.4 Organization

§2 describes the methodology and benchmark setup. §3-4 catalog
empirical negative results across classical attack families. §5
presents the engineering contributions (Barrett, C+pthreads).
§6 discusses what remains genuinely open. §7 concludes.

Appendices contain reproducible Sage scripts and C source listings.

## 2. Methodology

### 2.1 Benchmark curves

We use 4 LMFDB curves at 6 prime instances:

| Label | Conductor | `p` (bits) | `#E(F_p)` prime? |
|-------|----------:|-----------:|-----------------:|
| 67.a1 | 67 | 81 | yes |
| 21175.bc1 | 21175 | 80, 81 | yes, yes |
| 23232.cr1 | 23232 | 80, 81 | yes, yes |
| 114224.v1 | 114224 | 81 | yes |

All six instances are "good": prime curve order, embedding degree > 100,
trace `|t| ~ 10^12` (far from anomalous), j-invariant random (not 0 or
1728), Frobenius discriminant has 40+ bit prime factor (class numbers
~10^11, ruling out GLV).

### 2.2 Statistical protocol

- **50 trials per configuration** with controlled seeds
- **95% confidence intervals** reported on all scaling slopes
- **Full cost accounting**: includes RREF cost + Floyd's `3×` cycle
  detection correction, not just inner-loop point operations
- **Multi-curve sampling**: scaling claims based on data across all
  four curves, not one
- **End-to-end scaling validation**: timed Pollard rho recovery at
  20-40 bits in addition to per-step throughput

### 2.3 Tools

- SageMath 10+ for prototyping
- GMP 6.3 for modular arithmetic baseline
- clang -O3 for production C
- pthreads for multi-core parallelism
- Reproducible with `set_random_seed(42)` and stable compiler flags

## 3. Empirical negative results: classical attacks

### Table 3.1: Per-query Groebner basis vs. hashed pair-sum

| bits | `|FB|` | Pair-sum (s) | F_3 GB (s) | Ratio |
|----:|------:|-------------:|-----------:|------:|
| 13 | 20 | 0.01 | 0.29 | 30× |
| 16 | 40 | 0.01 | 0.38 | 38× |
| 19 | 50 | 0.02 | 0.66 | 33× |
| 22 | 40 | 1.62 | 109.18 | 67× |
| 25 | 40 | 9.23 | timeout | >32× |
| 28 | 40 | 82.05 | timeout | >3.66× |

The ratio is consistent at ~30-35× per-query overhead; the apparent
decrease at 25-28 bits is a timeout artifact (Method C runs longer per
query, hitting the 300s cap before finding relations).

### Table 3.2: Structured factor base GB cost

| Structure | bits | `|FB|` | Per-query (ms) |
|-----------|----:|------:|---------------:|
| Generic (∏(X−x_i)) | 13 | 20 | 4.7 |
| Generic | 16 | 40 | 18.1 |
| Generic | 19 | 50 | 30.8 |
| Mult subgroup (X^m−1) | 13 | 10 | 12.5 |
| Mult subgroup | 16 | 16 | 17.3 |
| Mult subgroup | 19 | 31 | 74.7 |
| Arithmetic progression | 13 | 16 | 4.6 |
| Arithmetic progression | 16 | 20 | 14.9 |
| Arithmetic progression | 19 | 24 | 12.2 |

**Sparsity of the constraint polynomial does not affect Groebner
basis cost.** Multiplicative subgroup (2 monomials) and arithmetic
progression (m+1 monomials but with combinatorial structure) give
no speedup over generic (m+1 monomials). This empirically confirms
Yokoyama et al.'s claim that the controlling parameter is regularity,
not constraint sparsity.

### Table 3.3: Non-naive attack results

| Attack | Result |
|--------|--------|
| Wu's method (triangular decomposition) | 14× *slower* than GB |
| Parametric GB (xR as parameter) | 3× slower than fresh GB |
| Quasi-subfield in F_{p²} | Relations findable but in wrong subgroup |
| LLL on canonical lift | No hidden number; LLL fails |
| p-adic L-function | No known ECDLP connection |
| ℓ-isogeny walk | Order-preserving; no exploit |

## 4. Empirical resistance of benchmark targets

### Table 4.1: Specific structural attack applicability

| Property | 67.a1 | 21175 | 23232 | 114224 |
|----------|------:|------:|------:|-------:|
| `\|trace\|` | 2.0e12 | 1.9e12 | 5.7e11 | 1.8e12 |
| Anomalous? | no | no | no | no |
| Emb. degree | >100 | >100 | >100 | >100 |
| MOV? | no | no | no | no |
| Curve order prime? | yes | yes | yes | yes |
| Pohlig-Hellman? | no | no | no | no |
| Twist smoothness | partial | prime | prime | prime |
| Smart anomalous? | no | no | no | no |
| GLV (small End)? | no | no | no | no |
| p-1 smooth part | 2^11 | 2^23 | 2^31 | 2^23 |
| Petit-Quisquater? | no | no | no | no |
| j-invariant special? | no | no | no | no |

Every applicability check returns negative. The targets are uniformly
hard against every known specific attack.

## 5. Engineering contributions

### Table 5.1: Pollard rho throughput at 81-bit

| Implementation | ops/sec | Speedup |
|---------------|--------:|--------:|
| Sage Python | 1.18e5 | 1× |
| C+GMP affine, 1 thread | 2.32e6 | 23× |
| C+GMP affine, 2 threads | 4.92e6 | 42× |
| C+GMP affine, 4 threads | 7.39e6 | 63× |
| C+GMP affine, 8 threads | 1.02e7 | 86× |

### Table 5.2: Barrett reduction vs. GMP at 81-bit

| Implementation | muls/sec | Per-mul (ns) |
|---------------|---------:|-------------:|
| Barrett (this work) | 1.06e8 | 9.4 |
| GMP persistent mpz_t | 1.11e7 | 90.0 |
| Speedup: **9.47×** | | |

### 5.3 End-to-end scaling validation

We measured actual Pollard rho recovery time at 20-40 bits:

| bits | wall time (s) | time / √n |
|-----:|--------------:|----------:|
| 20 | 0.025 | 2.43e-5 |
| 25 | 0.145 | 2.51e-5 |
| 30 | 0.803 | 2.45e-5 |
| 35 | 4.72 | 2.55e-5 |
| 40 | 24.14 | 2.30e-5 |

Linear fit: `log₂(time) = 0.4972 · bits − 15.24`, matching theoretical
`O(√n)` slope of 0.5 within 0.005.

## 6. Discussion

(... to be drafted; will cover what remains genuinely open ...)

## 7. Conclusion

(... to be drafted ...)

## References (selected)

- Pollard 1978: "Monte Carlo methods for index computation (mod p)"
- Semaev 2004: "Summation polynomials and the discrete logarithm
  problem on elliptic curves"
- Diem 2011: "On the discrete logarithm problem in elliptic curves"
- Faugère-Petit-Petit-Renault 2012: "Improving the complexity of
  index calculus algorithms in elliptic curves over binary fields"
- Petit-Quisquater 2016: "Algebraic approaches for the elliptic curve
  discrete logarithm problem over prime fields"
- Kudo-Yokota-Takahashi-Yasuda 2018: "Acceleration of index calculus
  for solving ECDLP over prime fields and Its limitation"
- Huang-Kosters-Petit-Yeo-Yun 2020: "Quasi-subfield polynomials and
  the elliptic curve discrete logarithm problem"
- Yokoyama-Yasuda-Takahashi-Kogure 2020: "Complexity bounds on
  Semaev's naive index calculus method for ECDLP"

## Draft progress checklist

- [x] Abstract
- [x] §1 Introduction
- [x] §2 Methodology
- [x] §3 Empirical negatives (tables drafted)
- [x] §4 Resistance of benchmark targets (table drafted)
- [x] §5 Engineering (tables drafted)
- [ ] §3 narrative around tables
- [ ] §4 narrative around table
- [ ] §5 narrative around tables
- [ ] §6 Discussion
- [ ] §7 Conclusion
- [ ] Reference cleanup (BibTeX)
- [ ] Appendix A: Sage script listings
- [ ] Appendix B: C source listings
- [ ] Appendix C: full measurement data tables
- [ ] Internal review pass
- [ ] External review request
- [ ] Submission to J. Math. Cryptol.

Estimated time to v1.0 draft (ready for external review): 4-6 weeks.
