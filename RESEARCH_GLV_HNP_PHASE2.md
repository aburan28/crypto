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

- ~~What is the **information-theoretic lower bound** on signatures
  needed for the k₁-only-leak model?  Phase 1 brute force needs
  `~ log_2(n) / c` signatures; does Phase 2's lattice achieve
  this in practice?~~  **Answered 2026-08-03 (§10): no.  The binding
  constraint is lattice geometry, not counting.**
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

## 10. Viability criterion and the centring correction (2026-08-03)

Measured over ≈2600 LLL/BKZ reductions on 17-bit and 12-bit j=0 GLV curves.
Scripts: `secp256k1_cm_audit/glv_hnp_phase2_centered{,_addendum,_e3fix}.py`.
Full data and derivation: `RESEARCH_AUTOLAB_LOG.md`, entry 2026-08-03.

### 10.1 When does the Phase-2 lattice work?

Let `eff = K1·K2/n` (K1 = k₁ bound, K2 = k₂ bound ≈ √n), m = signatures.
With the §2 scalings S_K1 = n/K1, S_K2 = n/K2, S_D = 1, S_KANNAN = n:

```
det(L)          = n^(2m+1) · eff^(-m)                     D = 2m+2
E‖v_planted‖²   = m·(K1²/a)·S_K1² + m·(K2²/a)·S_K2²
                  + (n²/3)·S_D² + S_KANNAN²
gap             = GH(det, D) / ‖v_planted‖
```

`a = 3` for limbs drawn as k ∈ [0,K); `a = 12` if they are centred (§10.2).
**Recovery succeeds iff gap ≳ 1** — measured 0.00 success below gap 0.5,
0.47 at gap 0.85–1.00, 0.85 at 1.00–1.20, 1.00 above 1.5.

Asymptotically gap → √a / √(2πe·eff), so the ceiling is `eff ≈ a/(2πe)`:
**0.176 uncentred, 0.703 centred.**  This is far below the counting bound
`eff < 1`, so Phase 2 is limited by lattice geometry, not by information.

Two caveats, both load-bearing:

- gap approaches its asymptote **from below** and only logarithmically in m,
  so at the m ≈ 12 used throughout Phase 2 the realised ceiling is much
  lower — gap = 1 lands near eff ≈ 0.07 uncentred, ≈ 0.21 centred.
- the criterion is calibrated for **LLL at D ≈ 26**.  Above that, LLL's
  approximation factor degrades about as fast as the gap improves, and
  success stops tracking gap: at eff = 0.40, m = 32 (gap 1.05) recovery is
  14/20, while at eff = 0.55, m = 96 (gap 1.05) it is 4/20.  Adding
  signatures past m ≈ 32 bought nothing measurable.

### 10.2 Centring the nonce limbs (correction to §2)

The §2 lattice plants `k1_i ∈ [0,K1)` and `k2_i ∈ [0,K2)` directly.  Absorb
`K1/2 + λ·K2/2` into the constants `A_i` instead, so the planted limbs are
`u_i = k1_i − K1/2` and `v_i = k2_i − K2/2`.  This is the standard HNP
centring step (Boneh–Venkatesan; Nguyen–Shparlinski) — **not a new technique;
the §2 construction simply omitted it.**  It cuts E[limb²] from K²/3 to
K²/12, halving ‖v_planted‖ asymptotically and buying 4× in eff (2.5–2.7× at
m = 12).  Measured, m = 12, 20 curves × 5 seeds:

| eff | §2 as written | centred |
|---|---|---|
| 0.10 | 43/100 | 99/100 |
| 0.15 | 18/100 | 91/100 |
| 0.20 | 10/100 | 65/100 |
| 0.25 | 5/100 | 47/100 |
| 0.40 | 0/100 | 8/100 |

On the (p=2677, n=2647, λ=185) curve the K1 wall moves from ≈4–6 to ≈16–24.
**Centring should be treated as part of the construction from here on.**

### 10.3 What does *not* help

- **Projecting out the d column.**  The trivial vector `n·S_D·e_m` is always
  shorter than the planted vector, so the planted vector is never λ₁.
  Quotienting it out is a **bit-exact no-op** — identical results seed by
  seed, 0/200 paired points differing.  LLL already isolates that vector and
  reduces inside its orthogonal complement, which *is* the projected lattice.
  Being λ₁ was never the binding constraint.
  Trap: the GH gap computed on the *projected* lattice is ≈1.6× optimistic.
  Always score with the unprojected gap.
- **BKZ(20) past the wall.**  At eff 0.70/0.85 (gap 0.57/0.52), LLL 0/40 and
  BKZ(20) 0/40.  Untested in the interesting regime gap ∈ [0.85, 1.05].
- **More signatures**, beyond the logarithmic effect noted in §10.1.
