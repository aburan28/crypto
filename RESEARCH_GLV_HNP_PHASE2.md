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

## 8a. Measured results (2026-07-29)

Scripts: `secp256k1_cm_audit/glv_hnp_phase2_lambda_threshold.py`,
`secp256k1_cm_audit/glv_hnp_phase2_lambda_confirm.py`.

### 8a.1 λ/n is *not* a viability parameter — falsified in both directions

An earlier reading of four data points suggested the Phase-2 attack succeeds
only when λ/n exceeds a threshold in (0.07, 0.34).  A controlled cohort study
falsifies this.  Holding the bias efficiency `eff = K1·K2/n` pinned at ≈0.157
and varying only λ, with 12 seeds and 12/12 required for success:

| curve      | λ/n    | λ/n-rule predicts | measured | recovery tally (wins/12)          |
|------------|--------|-------------------|----------|-----------------------------------|
| n=2551     | 0.0196 | FAILURE           | 12/12 at m=8 | m=4:3 5:6 6:9 7:11 8:12       |
| n=2647     | 0.0699 | FAILURE           | never    | peaks at 3/12 (m=7), 0/12 at m=14 |
| n=2251     | 0.3145 | SUCCESS           | never    | peaks at 4/12 (m=10)              |
| n=2293     | 0.4313 | SUCCESS           | 12/12 at m=8 | m=5:3 6:6 7:10 8:12           |

n=2551 succeeds at λ/n = 0.0196 — a factor 3.6 *below* the lower end of the
claimed threshold — while n=2251 fails at λ/n = 0.3145, comfortably inside the
claimed "safe" band.  Across a 16-curve 12-bit cohort and a 13-curve 20-bit
cohort at pinned `eff`, the success and failure ranges of λ/n overlap almost
completely (12-bit: success [0.020, 0.431] vs failure [0.035, 0.315]).

Two further points make the original claim untenable as *stated*:

1. The two roots of `x²+x+1 ≡ 0 (mod n)` are λ and `n−1−λ`, and the lattice
   entry `−λ·S_K1` is reducible modulo the row `n·S_K1`.  The lattice therefore
   only sees `λ_eff = min(λ, n−λ) ≤ n/2`.  A threshold stated on λ/n ∈ (0,1) is
   not well defined; the "independent" successes at λ/n = 0.66 and λ/n = 0.34
   are the same value of λ_eff/n.
2. Fixing one curve (λ/n = 0.1601, constant) and sweeping only K1 flips the
   outcome from 3/3 at m=4 to no recovery at m≤12 (§8a.2).  λ/n cannot be the
   controlling variable if recovery moves while λ/n is pinned.

Four further candidate predictors were tested and also fail to separate the
success and failure groups: `λ_eff/n`; `μ/‖v‖` where μ is the exact λ₁ of the
scaled 2-D sublattice `{((c·λ mod ±n)·S_K1, c·S_K2)}` (Lagrange-reduced) and
‖v‖ the planted-vector norm; `r_min/K1` with `r_min = min_{0<c<K2}|c·λ mod ±n|`;
and the ambiguity margin `min_c max(|r|/K1, c/K2)`.  All overlap.  This is
consistent with the 2026-06-27/06-28 result that the Phase-2 K1-threshold is a
continuous, curve-specific quantity with no simple algebraic separator, and
extends that negative result from 19-bit to the 12-bit and 20-bit regimes.

### 8a.2 Answer to the open question in §8: LLL needs ≈2× the information-theoretic minimum

Single curve `p=524341, b=22, n=525583, λ=84143` (λ/n = 0.1601 on every row),
K2 = 725, K1 swept.  `m_thresh = ln(n)/ln(1/eff)` is the counting bound;
`m*` is the smallest m with 3/3 recovery.

| K1  | eff     | m_thresh | m*   | m*/m_thresh |
|-----|---------|----------|------|-------------|
| 4   | 0.00552 | 2.53     | 4    | 1.58        |
| 8   | 0.01104 | 2.92     | 4    | 1.37        |
| 16  | 0.02207 | 3.45     | 6    | 1.74        |
| 36  | 0.04966 | 4.39     | 8    | 1.82        |
| 72  | 0.09932 | 5.70     | 12   | 2.10        |
| 144 | 0.19864 | 8.15     | >12  | >1.47       |
| 288 | 0.39727 | 14.27    | >12  | >0.84       |

So Phase 2 does **not** achieve the information-theoretic bound: it needs
roughly `1.4–2.1 × m_thresh` signatures, and the multiplier grows as the bias
weakens.  The mechanism is a race — the lattice dimension is `2m+2`, so demanding
more signatures to compensate for a weaker bias simultaneously degrades LLL's
approximation factor.  Past `eff ≈ 0.2` the two effects cancel and plain LLL
stops recovering at any m in range.  Note this is a statement about the
*attack*, not about the curve: `eff` is set by how much nonce bias the
implementation leaks, so it is an implementation property, not a curve property.
No curve-intrinsic invariant found so far predicts viability.

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
