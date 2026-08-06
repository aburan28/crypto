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

## 8b. Update 2026-08-06 — calibration of the viability condition

Empirical work on toy curves (8–17 bit) has now settled how far this
lattice reaches.  Script:
[`secp256k1_cm_audit/glv_hnp_phase2_centered_cvp.py`](secp256k1_cm_audit/glv_hnp_phase2_centered_cvp.py);
full log entry under `## 2026-08-06` in
[`RESEARCH_AUTOLAB_LOG.md`](RESEARCH_AUTOLAB_LOG.md).

**Center the unknowns.**  The implementation used from 2026-06-15 to
2026-08-05 encoded `k₁ᵢ ∈ [0,K₁)`, `k₂ᵢ ∈ [0,K₂)`, `d ∈ [0,n)` one-sided.
Lattice reduction searches around the origin, so subtracting each
interval's midpoint from the target drops `E[x²]` from `X²/3` to `X²/12`
— a factor 2 in the target vector's length, a factor 4 in tolerable bias
strength.  Measured: pooled recovery on six fresh 17-bit curves (m=12,
5 seeds, 120 trials) rises from **16/120 to 56/120**.

**The viability condition is the Gaussian heuristic, nothing curve-specific.**
With `N = 2m+1`, `det(L) = (n·S_K1)^m · S_D · S_K2^m` and centered
`E‖e‖`, the ratio

    r_GH = E‖e‖ / ( sqrt(N/2πe) · det(L)^(1/N) )

separates success from failure at **r\* = 1.00** with AUC 0.994 over 39
(curve, K₁) cells.  Writing `eff = K₁K₂/n`, the closed form is

    r_GH  =  sqrt(2πe/12) · eff^(m/(2m+1)) · n^(1/(2m+1))
          ≈  1.194 · eff^(m/(2m+1)) · n^(1/(2m+1))

(verified against `gh_ratio()` to 3 decimals: m=8, n=2659, K₁=16, K₂=52 →
1.099 predicted vs 1.098 computed).  Two regimes fall out:

- **m → ∞:** `r_GH → 1.194·sqrt(eff)`, so the ceiling is `eff ≲ 0.70`,
  i.e. `K₁·K₂ ≲ n/1.4`.  The uncentered lattice's ceiling was 4× lower,
  `eff ≲ 0.18`.
- **Finite m:** the `n^(1/(2m+1))` factor is a real penalty, not a
  rounding term.  At m=12 and a 17-bit `n` it is 1.56, which is why the
  measured 17-bit wall sits near `eff ≈ 0.25–0.35` rather than 0.70.  For
  a 256-bit `n` the factor is `2^(256/(2m+1))`, so **`m ≳ 128` signatures
  are needed before the asymptotic ceiling is even approached** — this is
  the dominant cost term for any real-curve instantiation, and it is not
  reflected in the `m ≳ 5` estimate in §2 above, which should be treated
  as superseded.

This answers the first open question above in the negative: the lattice
does **not** reach the information-theoretic line (`K₁K₂ < n`); it stops
at the standard Gaussian-heuristic wall, a fixed constant factor short,
and BKZ-20 did not move it on any tested instance.

**Two hypotheses closed as dead.**  (i) `λ* = min(λ, n−λ)/n` is not a
predictor — the apparent λ*-dependence of the K₁ wall was an artifact of
the uncentered encoding.  (ii) The planted vector is never `λ₁` of this
lattice (the trivial vector `n·S_D·eₘ` is always shorter), but that is
harmless: it lies in the `d mod n` periodicity direction and never
competes for the BDD solution.  Projecting it out and solving CVP
directly performs identically to the centered Kannan embedding.

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
