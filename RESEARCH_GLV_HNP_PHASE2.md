# Phase 2: GLV-aware Lattice Construction for the k₁-only-Leak Model

A concrete proposal for the Phase 2 lattice attack in the GLV-HNP
research direction.  Builds on Phase 1 (HNP framework demo) and
Phase 1.5 (standard Boneh-Venkatesan LLL attack on P-256) by
specialising to the **k₁-only-leak threat model** where the GLV
half-scalar `k₁` is partially known but `k₂` is treated as a free
variable.

> **Status**: design document.  Implementation is contingent on
> the secp256k1 LLL-degeneracy resolution (see
> [`tests/lll_degeneracy_probe.rs`](tests/lll_degeneracy_probe.rs)).
> Phase 2 likely inherits the same lattice geometry, so the
> degeneracy investigation in (α) must complete first.

## 1. The threat model

### Standard HNP threat model (Phase 1.5)

Attacker observes biased ECDSA signatures with a known bound on
the FULL nonce `k`:

```
        0  ≤  k_i  <  2^k_bits     (k_bits < n_bits)
```

Bias source: truncated RNG output, modular-reduction-by-subtraction
implementations, side-channel leakage of the top bits of `k`.

The standard Boneh-Venkatesan lattice attack works.

### k₁-only-leak threat model (Phase 2)

Attacker observes biased ECDSA signatures where the **GLV
half-scalar `k₁` is partially known but `k₂` is NOT bounded**:

```
        0  ≤  k_{i,1}  <  2^c     (c < 128)
        0  ≤  k_{i,2}  <  2^128   (uniform over 128-bit half)
        k_i  =  k_{i,1}  +  λ · k_{i,2}     (mod n)
```

Bias source: side-channel attack on a GLV-implementing ECDSA
signer that processes `k₁` and `k₂` in separate code paths
(possibly with separate timing or cache footprints).  Real
implementations: libsecp256k1's `secp256k1_scalar_split_lambda`,
the GLV path in OpenSSL EC_POINT_mul, hardware-accelerated
secp256k1 multipliers.

The standard HNP attack **does not work** here:
- The "implicit `k` bias" is `k_{i,1} + λ · k_{i,2}` mod `n`,
  with `λ · k_{i,2}` uniform over (Z/n) since `k_{i,2}` is
  uniform.
- Adding a small known value to a uniform value gives a uniform
  value.
- The lattice short-vector hypothesis breaks: no signature has a
  "short" projection useful to LLL.

## 2. The GLV-aware lattice

### Signature equation in GLV form

```
        k_{i,1}  +  λ · k_{i,2}  ≡  A_i  +  B_i · d     (mod n)
```

Define `T_i := A_i + B_i · d  −  λ · k_{i,2}  (mod n) = k_{i,1}`.
For each signature we have:

```
        T_i  =  A_i + B_i · d  −  λ · k_{i,2}       (mod n)
        with constraint:  0  ≤  T_i  <  2^c
```

This is a SINGLE congruence per signature with TWO unknowns:
the secret `d` (shared across all signatures) and the half-scalar
`k_{i,2}` (different per signature).

### Lattice construction

For `m` signatures, build a lattice `L` over `Z` of dimension
`m + 2` (one column per signature, one for `d`, one for `λ-times-d-row`):

```
   Columns         1 ... m     m+1      m+2
   Rows
   1..m  (mod n)   n·I_m       0        0
   m+1   (d-row)   B_1..B_m    2^c     0
   m+2   (λ-row)  -λ·I'_m      0        2^c·1
   m+3..2m+2       (k_{i,2})  
                  -λ·e_i      0        0
```

(Schematic; needs careful indexing.)

The key insight: the `λ` enters as a known scalar multiplying
unknown half-scalars `k_{i,2}`.  Each `k_{i,2}` becomes its own
lattice unknown, expanding the search space from `(d)` (1 dim) to
`(d, k_{1,2}, k_{2,2}, ..., k_{m,2})` (m+1 dims).

### Dimension analysis

- Lattice has `2m + O(1)` dimensions
- Each `k_{i,2}` is bounded by `2^128` (half-curve range)
- The "short" target vector has `m` entries of size `O(2^c)`
  (the `T_i` values), `1` entry of size `O(2^c)` (`d`), and `m`
  entries of size `O(2^128)` (the `k_{i,2}`).
- Total target-vector norm: `√(m · 2^{2c} + 2^{2c} + m · 2^{256})`
  ≈ `√m · 2^{128}`

For LLL to succeed, this norm must beat the lattice's typical
short-vector length.  Generic LLL gap: `2^{(dim − 1)/4}` × det^{1/dim}.

Det of the (2m+2)-dim lattice ≈ `n^m · 2^{2c+128}` (volume of
the modular constraints + the bounds).

So the gap requires:

```
        m  ≳  2 + (n_bits − c) / 64    (rough)
```

For `c = 64` (top 64 bits of `k_1` known, 64 bits free), `n_bits =
256`:
```
        m  ≳  2 + 192 / 64  =  5
```

So **5–10 signatures should suffice** in the k₁-only-leak model
— provided the lattice doesn't suffer the same LLL-degeneracy as
the standard Phase 1.5 setup on secp256k1.

This signature count is MUCH lower than the standard HNP equivalent
in the k-overall-bias model (where for `c_total = 64` we'd need
log_2(n)/c_total ≈ 4 sigs, but with strong-bias assumption that
doesn't apply to the k₁-only case).

## 3. Why this is novel

To the best of our knowledge, no published HNP variant exploits
GLV-half-scalar structure separately from the full nonce.  Related
work:

- Standard HNP (Boneh-Venkatesan, Howgrave-Graham-Smart, Nguyen-
  Shparlinski): treats `k` as a single unknown, no GLV awareness.
- Galbraith-Lin-Scott "GLS" curves: a different way to use
  endomorphisms in EC arithmetic, no relevance to nonce attacks.
- Bleichenbacher: works on fractional bias of `k` via DFT; doesn't
  exploit decomposition.
- BBSP (Aranha et al "LadderLeak"): sub-bit-leak HNP via DFT; also
  not GLV-aware.

The closest published prior work: **multi-key HNP** (Bos et al,
Aranha-Cabral-Lopez), which exploits shared bias across multiple
keys.  Multi-key HNP could in principle be combined with GLV-aware
HNP for further amortisation, but to our knowledge this combination
has not been studied.

## 4. Risks and contingencies

### Risk 1: lattice degeneracy

The Phase 1.5 LLL-degeneracy investigation
([`tests/lll_degeneracy_probe.rs`](tests/lll_degeneracy_probe.rs))
revealed that secp256k1's standard Boneh-Venkatesan lattices
systematically fail LLL convergence at `k_bits = 192`.  If the
GLV-aware lattice inherits this degeneracy, Phase 2 needs to
either:

- Use BKZ-β for β ≥ 20 (much slower but more robust)
- Use a fundamentally different lattice formulation
- Wait for the (α) investigation to identify and patch the
  underlying basis-construction issue

**Mitigation**: prototype Phase 2 on P-256 first (where standard
LLL works), then port to secp256k1 only after the degeneracy is
resolved.

### Risk 2: λ size

`λ` for secp256k1 is `≈ 2^256` (a full-size 256-bit value).  The
lattice entries `λ · n` are `≈ 2^512`.  LLL on 512-bit-entry
bases is much slower than on 256-bit (cost scales as entry-size²
in the Gram-Schmidt update).  Expect Phase 2 LLL runs to take
roughly 10× the Phase 1.5 time.

**Mitigation**: reduce `λ` modulo `n` (it's already < n by
definition).  Use BKZ-15 or BKZ-20 instead of large-β BKZ; the
extra cost is manageable.

### Risk 3: the bias model doesn't match real attacks

The k₁-only-leak model assumes side-channel leakage on `k_1` but
not `k_2`.  In practice, an implementation that leaks `k_1` likely
leaks `k_2` similarly (symmetric code paths).  Phase 2's "novelty"
might be inapplicable to real deployed implementations.

**Mitigation**: Phase 3 (real-world bias profiling) should
characterize which real implementations have asymmetric leakage.
The proposal stands or falls on Phase 3's empirical findings.

## 5. Implementation plan

```
Step 1:  Build the GLV-aware lattice (Section 2.2) in PARI for a
         small toy curve (n ≈ 2^60).
Step 2:  Generate planted-d ECDSA signatures with controlled
         k₁-only bias.
Step 3:  Run LLL on the GLV-aware lattice; verify d recovery.
Step 4:  Compare convergence to standard HNP on the same data.
Step 5:  Scale up to n ≈ 2^192 / 2^256 with BKZ if needed.
```

Estimated effort: 2–3 weeks once the (α) degeneracy investigation
is complete.

## 6. Connection to other research directions

### vs. structural-completeness theorem

GLV-HNP attacks the **scalar side** (recover `d`).  Structural-
completeness covers the **group side** (no isogeny-graph attack
beats ρ).  These are independent threat models — Phase 2 success
would NOT contradict the structural-completeness theorem.

### vs. CSIDH / SQIsign (post-quantum)

GLV-HNP is purely classical and pre-quantum.  CSIDH and SQIsign
use the supersingular-isogeny problem and are unrelated.

### vs. side-channel cryptanalysis broadly

GLV-HNP is one specific instance of side-channel cryptanalysis,
tailored to GLV implementations.  Other side-channel attacks
(LadderLeak, Minerva, TPM-FAIL) target different implementation
quirks.

## 7. Out-of-scope for Phase 2

- Real-world hardware target identification (Phase 3)
- End-to-end attack against a deployed system (Phase 4)
- Multi-key GLV-HNP batching (would be Phase 2.5)
- Quantum-aware lattice reduction (Phase 5+)
- GLS curves (different endomorphism structure; separate proposal)

## 8. Open questions

- What is the **information-theoretic lower bound** on signatures
  needed for the k₁-only-leak model?  Phase 1 brute force needs
  `~ log_2(n) / c` signatures; does Phase 2's lattice achieve
  this in practice?
- Can Phase 2 be **further specialised** to exploit secp256k1's
  CM-by-Z[ω] structure beyond just GLV?  E.g., if the
  endomorphism ω acts predictably on the bias distribution, a
  3D-aware HNP variant could exploit a third dimension.
- Does the **secp256k1 LLL-degeneracy** generalize to the GLV-
  aware lattice?  If yes, the entire Phase 2 needs a different
  base reduction algorithm.

## 9. References

- Gallant, Lambert, Vanstone, "Faster point multiplication on
  elliptic curves with efficient endomorphisms", CRYPTO 2001.
- Boneh, Venkatesan, "Hardness of computing the most significant
  bits of secret keys in Diffie–Hellman and related schemes",
  CRYPTO 1996.
- Nguyen, Shparlinski, "The insecurity of the digital signature
  algorithm with partially known nonces", J. Cryptol. 2002.
- Bauer, Naccache, "Toward a rigorous variation of Coppersmith's
  algorithm on three variables", Eurocrypt 2007.
- Aldaya, Pereira, Brumley, Tuveri, García-Mariscal, Pereida-García,
  "LadderLeak: breaking ECDSA with less than one bit of nonce
  leakage", CCS 2020.

## 10. Addendum — the ν̂ invariant in closed form (Thread 23, 2026-07-29)

Thread 20d found that the rival-sublattice invariant

```
    ν̂ = λ₁(L₂)/√(det L₂),   L₂ = ⟨(n·S₁, 0), (−λ·S₁, S₂)⟩,
    S₁ = ⌊n/K₁⌋,  S₂ = max(1, ⌊n/K₂⌋)
```

separates the June C1/C2 classes (AUC 0.935) where eight scale-free invariants
had failed.  Thread 23 replaces the computation by a closed form.  Script:
[`secp256k1_cm_audit/glv_hnp_nuhat_closed_form.py`](secp256k1_cm_audit/glv_hnp_nuhat_closed_form.py),
output `nuhat_closed_form_output.txt`.

### 10.1 Closed form

`L₂ = {(S₁(an − bλ), S₂b)}`, so for fixed `b` the inner minimum over `a` is
`S₁·‖bλ‖ₙ` (centered residue) and

```
    λ₁(L₂)² = min_{b ≥ 0} [ S₁²·‖bλ‖ₙ² + S₂²·b² ].
```

Any minimising `b > 0` satisfies `‖b'λ‖ₙ > ‖bλ‖ₙ` for every `0 < b' < b`
(otherwise `b'` beats `b`, since `S₂²b'² < S₂²b²`), i.e. `b` is a best
approximation of the second kind to `λ/n`; by Lagrange's theorem those are
exactly the continued-fraction convergent denominators `q_j`.  `b = 0` is
`q₋₁ = 0`, contributing `S₁n`.  So ν̂ is an **O(log n) continued-fraction
quantity**, no lattice reduction required.

Verified against exact-integer Lagrange–Gauss reduction: **1080/1080 exact
agreements** over bit lengths 20–521 × eff ∈ {0.02, 0.10, 0.50}.

The minimising convergent is *not* the largest partial quotient but the one at
the scale `R = √(n·S₁/S₂) = √(n·K₂/K₁)` where the two terms balance
(`θ ≈ n/q`).  `q_{j*}` is the convergent nearest `R` in 40–50 % of cases and
within one index in 93–99 %.  This is why the scale-free June invariants
(`q_cf`, `max_q_cf`, `max_a`) failed: at 20 bits, `log(1+a_local)` correlates
−0.74 with ν̂ (Spearman −0.75) while the global `log(1+a_max)` reaches only
−0.41 (−0.34), most of which it inherits from `a_local` (r = +0.61).

### 10.2 Null law

In normalised coordinates `u = q/R`, `w = θR/n`, `ν̂² = min(u² + w²)`, i.e.
`L₂/√det` is a unimodular 2-lattice.  Under equidistribution `ν̂² = 1/Im τ`
for `τ` Haar on `SL(2,ℤ)\ℍ`, giving the exact CDF

```
    P(ν̂ ≤ r) = 3r²/π                                             r ≤ 1
             = (3/π)[r² − 2r²√(1−r⁻⁴) − 2 arcsin(r⁻²) + π]       1 ≤ r ≤ (4/3)^¼
             = 1                                                 r ≥ (4/3)^¼
```

with support ceiling the Hermite bound `(4/3)^¼ = 1.074570`.  KS tests against
random λ: D = 0.0037–0.0125 (p = 0.17–0.95) at 64 and 256 bits, eff ∈ {0.05,
0.25}; no sample ever exceeded the Hermite ceiling.  This reproduces the
empirical null quantiles logged 2026-07-29 to three decimals.

### 10.3 Theorem (GLV floor)

Let `n` be prime with `λ² + λ + 1 ≡ 0 (mod n)` — i.e. the j = 0 GLV setting.
If `a + bλ ≡ 0 (mod n)` with `(a,b) ≠ 0` then

```
    a² − ab + b² ≡ b²(λ² + λ + 1) ≡ 0   (mod n),
```

and the form is positive definite, so `a² − ab + b² ≥ n`.  Minimising the
quadratic form `S₁²a² + S₂²b²` subject to that constraint (a generalised
Rayleigh quotient) gives, unconditionally,

```
    ν̂  ≥  √( t_min /(S₁S₂) ),   t_min = smaller root of
                                (3/4)t² − (S₁²+S₂²)t + S₁²S₂² = 0
       =  √( 2 / ( r[(1 + r⁻²) + √(1 − r⁻² + r⁻⁴)] ) ),   r = S₁/S₂ = K₂/K₁.
```

In particular `ν̂ ≥ √(2/3) = 0.8165` when `K₁ = K₂`, decaying as `√(K₁/K₂)`.
The bound is sharp: over 200 64-bit j = 0 curves the observed minima are
0.816499 / 0.681853 / 0.495975 / 0.249879 against floors 0.816497 / 0.681774 /
0.495955 / 0.249878 at `K₂/K₁` = 1 / 2 / 4 / 16, with **zero violations** at any
ratio tested (up to 2²⁰).

Consequence for the threat model: the C2 ("LLL recovers `d` on every seed")
band requires ν̂ ≤ 0.645, and the floor exceeds 0.645 whenever
`K₂/K₁ ≤ 2.276`.  **No j = 0 GLV curve can be in the easy class when the
leaked bound `K₁` is within 1.19 bits of `K₂`** — no CM hypothesis, no
reduction heuristic, just positive-definiteness of the norm form.

This is also visible statistically: the GLV λ is non-Haar exactly where the
floor bites.  KS vs the §10.2 law over 300 64-bit curves: D = 0.637 (p ≈ 0) at
`K₂/K₁ = 1`, D = 0.235 (p ≈ 0) at 4, but D = 0.034–0.047 (p = 0.52–0.87) for
`K₂/K₁ ≥ 256` — the arithmetic of Z[ω] is invisible to ν̂ once the leak is
strong, and the curve is then no safer than a random λ.

### 10.4 Precision note

The incumbent float64 `gauss_reduce_2d` (`glv_hnp_phase2_lambda_threshold.py:216`)
is accurate to ≤ 2.1e-16 relative up to 384 bits — Python's exact int/int
division rescues it — so every published ν̂ number stands, including the
secp256k1 placement (exact vs float agree to 1e-14).  At **521 bits it raises
OverflowError on 25/25 trials** (`math.sqrt` of a norm ≳ 2¹⁰²⁴), the same
failure mode as the P-521 LLL thread.  `nu_hat_exact()` in the Thread 23 script
is overflow-free and should be preferred.
