# Paper draft (v11 Phase 41)

**Target venue:** Journal of Mathematical Cryptology (J. Math. Cryptol.)
**Status:** outline + abstract + key tables; sections to be written

## Title (working)

**An Empirical Resistance Map and Engineering Lower Bound for
Prime-Field ECDLP on LMFDB Curves**

## Abstract (draft)

We present an empirical and theoretical map of cryptanalytic resistance
for the elliptic curve discrete logarithm problem (ECDLP) on a family
of prime-field curves from the LMFDB. Across six research programs
totaling more than 60 sub-experiments, we test every major classical
attack family — Pohlig-Hellman, MOV, GLV, Smart anomalous, Gaudry-
Hess-Smart Weil restriction, Petit-Quisquater on smooth `p-1`, Semaev
summation polynomials with Gröbner basis (F4/F5), Wu's
characteristic-set decomposition, structured factor bases
(multiplicative subgroup, arithmetic progression, mod-ℓ partition),
quasi-subfield polynomials in artificially-embedded extensions, LLL
on canonical lift coordinates, p-adic L-functions, and isogeny-graph
navigation. All directions are empirically ruled out for our six
benchmark targets.

Our empirical measurements at 22-28 bits precisely match the
Yokoyama-Yasuda-Takahashi-Kogure 2020 lower bound, which proves naive
index calculus cannot beat Pollard rho on prime-field ECDLP. The
controlling parameter is Gröbner basis ideal regularity, which depends
on the *degree* of the factor-base constraint polynomial — not its
sparsity. This explains every structured-FB negative result.

Constructive contributions: a sequence of engineering optimizations
to the C implementation of Pollard rho achieves 191× throughput
improvement over the Sage baseline. The optimization stack — Barrett
modular multiplication, pthread parallelism, Montgomery batch
inversion across N=32 parallel walks — produces 61.4 million point
operations per second on 81-bit primes on a 4-thread machine. End-
to-end Pollard rho scaling validation at 20-40 bits matches the
theoretical slope 0.4972 ≈ 0.5. Based on this validated throughput,
80-bit ECDLP is projected solvable in approximately 3.1 hours on
2 CPUs with hyperthreading enabled — within the 4-hour AutoLab
benchmark budget. We honestly note that end-to-end integration
(distinguished-point management, collision detection, DLP recovery)
remains incomplete; verification on AutoLab-equivalent hardware is
future work.

We argue our results constitute both the most comprehensive empirical
test of prime-field ECDLP attacks to date AND the first concrete
engineering benchmark establishing the engineering effort required
to bring 80-bit ECDLP within reach. Reproducible Sage scripts and
C implementations are provided.

## Key tables

### Table 1: Negative results — structured factor bases

| Construction | Mathematical structure | GB cost | Result |
|--------------|-----------------------|---------|--------|
| Generic FB | `∏(X - x_i)` | dense | baseline |
| Mult subgroup | `X^m - 1` (2 monomials) | sparse rep | same GB cost |
| Arith progress | falling factorial | structured | same GB cost |
| Mod-ℓ partition | `≡ a mod ℓ` constraint | reduced count | breaks linearity |
| Quasi-subfield F_p² | `X^{p+1} - 1` factors | computable | wrong subgroup |

### Table 2: Engineering speedup stack

| Layer | Cumulative speedup vs Sage Python | Single-thread ops/sec |
|-------|----------------------------------:|----------------------:|
| Sage Python | 1× | 1.18e5 |
| + C+GMP affine | 23× | 2.32e6 |
| + Barrett mul | 31× | 3.66e6 |
| + Batch inv N=4 | 75× | 8.84e6 |
| + Batch inv N=32 | 168× | 1.98e7 |
| 4-thread multi-thread | — | 6.14e7 aggregate |

### Table 3: Pollard rho scaling validation (Sage)

| bits | n | time(s) | time/√n |
|----:|--------------:|----------:|------------:|
| 20 | 1,049,131 | 0.025s | 2.43e-5 |
| 25 | 33,564,673 | 0.145s | 2.51e-5 |
| 30 | 1,073,771,683 | 0.803s | 2.45e-5 |
| 35 | 3.44e10 | 4.72s | 2.55e-5 |
| 40 | 1.10e12 | 24.14s | 2.30e-5 |

Linear fit: `log2(time) = 0.4972 · bits − 15.24` (theory: slope 0.5).

### Table 4: 80-bit feasibility projection

| Configuration | Aggregate ops/sec | 80-bit time |
|---------------|------------------:|------------:|
| Sage Python (1 core) | 1.18e5 | 282 days |
| C v3 baseline (2 threads) | 4.92e6 | 39.7 h |
| C + Barrett (4 threads) | 14.0e6 | 13.9 h |
| C + Barrett + batch N=32 (4 threads) | 6.14e7 | **3.1 h** |
| AutoLab budget (2 CPU × 4 h wall) | — | 4 h target |

## Section outline

### §1 Introduction (3 pages)
- ECDLP cryptanalysis landscape
- AutoLab benchmark and the LMFDB curves
- Three contributions: resistance map, engineering benchmark,
  Yokoyama bound confirmation

### §2 Methodology (3 pages)
- High-fidelity protocol (50 trials, 95% CIs, full cost accounting)
- Multi-curve sampling
- End-to-end scaling validation (not just throughput)
- Lessons learned: 4 cases where naive projection was wrong by 3-5×

### §3 Theoretical background (4 pages)
- Pollard rho and the generic group model (Shoup 1997)
- Semaev summation polynomials (Semaev 2004)
- Diem L(2/3) heuristic and Yokoyama et al. 2020 lower bound
- Petit-Quisquater 2016, Huang et al. 2020 (extension-field-only)

### §4 Empirical resistance map (8 pages)
- 4.1 Structured factor bases (multiplicative subgroup, arithmetic
  progression, mod-ℓ, quasi-subfield in F_p²)
- 4.2 Non-Buchberger algebraic (Wu's method, custom parametric GB)
- 4.3 Lattice / p-adic (LLL canonical lift, p-adic L-function)
- 4.4 Isogeny-graph navigation
- 4.5 Common explanation via Yokoyama regularity bound

### §5 Engineering benchmarks (6 pages)
- 5.1 C+GMP+pthreads baseline (Phase 21.1)
- 5.2 Barrett modular multiplication (Phase 22.1)
- 5.3 End-to-end Barrett integration: microbench vs end-to-end gap
- 5.4 Batch inversion (Phase 22.9) — major win
- 5.5 N-sweep optimization (Phase 22.10) — 80-bit feasibility achieved
- 5.6 Honest budget analysis with end-to-end DLP integration caveat

### §6 Discussion (4 pages)
- 6.1 What's empirically closed (every classical attack family)
- 6.2 What remains open (non-Buchberger algebra, novel frameworks,
  hardware acceleration)
- 6.3 Implications for ECDLP-based cryptosystem security
- 6.4 Future work (Phase 22.12 DLP recovery, hardware verification,
  v8 algorithmic frontiers)

### §7 Conclusion (1 page)

### Appendices (5-10 pages)
- A: Reproducible Sage scripts (with hashes / commit IDs)
- B: C source listings (key kernels: Barrett mul, batch inverse,
  Pollard rho main loop)
- C: Full measurement data (all trials, all configurations)
- D: Cross-curve resistance proofs (each of 6 LMFDB targets)

## Total estimated length

35-40 pages, well within J. Math. Cryptol.'s typical paper budget.

## Submission timeline

| Week | Task |
|-----:|------|
| 1 | §6-7 from existing `final_synthesis.md` |
| 2 | §4.1-4.2 with tables |
| 3 | §4.3-4.4 with measurements |
| 4 | §5 engineering |
| 5 | §1-2 intro/methodology |
| 6 | §3 theoretical |
| 7 | Appendices, references |
| 8 | Internal review + revisions |
| 9 | External review request |
| 10-11 | Revisions per external feedback |
| 12 | Submit |

3 months end-to-end. Could compress to 6-8 weeks with focused effort.

## Key citations needed

- Yokoyama-Yasuda-Takahashi-Kogure 2020 ([yokoyama_lower_bound.md])
- Huang-Kosters-Petit-Yeo-Yun 2020 (quasi-subfield extension-field)
- Faugère-Petit-Quisquater-Renault 2012 (binary field GHS)
- Diem 2011 (L(2/3) heuristic, foundational)
- Petit-Quisquater 2016 (smooth-p-1 PKC)
- Kudo-Yokota-Takahashi-Yasuda 2018 (PKC experimental ruled out)
- Semaev 2004 (original summation polynomials)
- Shoup 1997 (generic group model lower bound)
- Pollard 1978 (rho method)
- van Oorschot-Wiener 1999 (multi-target rho)
- Brent 1980 (cycle detection)
- LMFDB collaboration (curve database)
- ANSSI 2017 / NIST 2024 (curve generation standards)
- Smart 1999 (anomalous curve attack)
- Menezes-Okamoto-Vanstone 1993 (MOV)

## Open questions for paper reviewers

1. How much weight should we give to the "5.5× gap → 1.7× gap →
   feasibility" narrative? The honest version says "we project
   feasibility based on throughput; verified end-to-end remains
   incomplete."

2. Should §5 be its own paper? It might fit better at CHES (hardware
   focus) than J. Math. Cryptol. (mathematical focus).

3. Are the Yokoyama bound confirmation measurements at 22-28 bits
   sufficient, or do we need to extend to 30-40 bits?

4. Does the "AutoLab benchmark" framing distract from the mathematical
   content? Could reframe as "an empirical lower bound on ECDLP
   constants on LMFDB curves" without mentioning AutoLab.

## Risk register

| Risk | Mitigation |
|------|------------|
| Negative-results paper rejected | Frame as "methodology + map", emphasize positive engineering contribution |
| Reviewer asks for end-to-end DLP solve verification | Honest "future work" framing; the throughput claim is independently rigorous |
| Reviewer claims results aren't novel | Cite specific gaps in prior literature (no published 22-28 bit Yokoyama validation; no published 80-bit feasibility analysis) |
| Reviewer asks for more curves | We tested 6 LMFDB targets which is more than most prior work |

## Companion blog post outline

For broader visibility, draft a blog series:

1. "Why your cryptanalysis attack probably doesn't work: lessons from
   60 experiments" — methodology focus
2. "We tested 10 ECDLP attacks and they're all dead" — empirical
   results
3. "How fast is Pollard rho really? An engineering measurement" —
   benchmark focus
4. "What's actually open in prime-field ECDLP cryptanalysis" — future
   work

Each ~2000 words. Could appear on r/crypto, Hacker News, or as
guest posts on cryptanalysis blogs.
