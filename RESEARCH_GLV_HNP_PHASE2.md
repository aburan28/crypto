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

## 6b. The ν̂ separator and its closed form (Threads 20c/20d/21)

After eight curve-level invariants were falsified as predictors of
whether LLL recovers `d` (δ/n, κ(M), q_cf, max_q_cf, max_a,
a_corn/n, λ/n, and μ = min(λ,n−λ)/n — see the autolab log,
2026-06-27 … 2026-07-29), the first invariant that does separate is
lattice-geometric.

**Definition.** The (2m+2)-dimensional Phase 2 lattice contains, on
each coordinate pair `(i, m+1+i)`, a copy of the *non-planted* 2-D
sublattice

```
L2 = ⟨ (n·S_K1, 0), (−λ·S_K1, S_K2) ⟩,     det L2 = n·S_K1·S_K2.
```

`det L2` does not depend on λ, so

```
ν̂ = λ₁(L2) / sqrt(det L2)
```

isolates exactly the λ-dependence of the geometry in one
Lagrange–Gauss reduction, `O(log n)`.  Empirically (Thread 20d,
120 fresh 20-bit curves, replicated 2026-08-04): **AUC 0.918** against
the June C1/C2 classes at eff = 0.0993, m = 12, versus 0.523 for μ and
0.748 for max_a.

**Sign.** Opposite to the naive guess: *small* ν̂ makes the attack
*easier*.  With `λ₁·λ₂ ≈ det`, a small λ₁ means L2 is skew, so its
second minimum is unusually long and the planted vector is
comparatively short; a balanced L2 (ν̂ → 1) surrounds the planted
vector with 2m rivals of norm ≈ sqrt(det).

**Closed form (Thread 21).**  Strip the scaling: `L2 = D·Λ` with
`D = diag(S_K1, S_K2)` and `Λ = {(x,b) ∈ Z² : x + λb ≡ 0 mod n}`.
For a j=0 curve `n = u² − uv + v²` is an Eisenstein norm and Λ is the
coordinate image of the principal ideal `(π)`, `π = u + vω`, w.r.t.
the basis `{1, ω}`:

```
Λ = A_π Z²,   A_π = [[u, −v], [v, u−v]],   det A_π = n.
```

Multiplication by π is (scale by sqrt(n)) ∘ (rotate by θ = arg π) in
the Minkowski plane, so with `C = [[1, −1/2], [0, sqrt(3)/2]]` the
n-dependence cancels **exactly**:

```
ν̂ = λ₁( G(τ,θ) Z² ) / sqrt(τ),
G(τ,θ) = diag(τ,1) · C⁻¹ · R_θ · C,      τ = S_K1/S_K2 ≈ 1/eff.
```

The curve enters *only* through the CM angle θ = arg(u + vω), and n
does not appear at all.  Verified to 1.1e-14 on 1600 (curve, root,
eff) triples — `secp256k1_cm_audit/glv_hnp_nuhat_closed_form.py`.
Consequences:

- ν̂ is 60°-periodic in θ (the six units of Z[ω]); the two GLV roots
  sit at θ and 60°−θ and give *different* ν̂ (no reflection symmetry).
- θ is equidistributed on [0°,60°) — χ² = 2.9/12.9/5.1/3.5 on 10 bins
  at 20/24/32/48 bits (5% critical value 16.92).
- ν̂ for secp256k1 is **exact**, not a 20/24-bit extrapolation:
  θ = 9.448519°, and ν̂ = 0.8782 / 0.8709 / 0.6639 / 0.5990 / 0.5852
  at eff = 0.02 / 0.05 / 0.0993 / 0.15 / 0.25.

**Why a_corn/n failed but ν̂ works.**  Both come from the same
Cornacchia datum.  `a_corn = |u|` is the *radial* part, but
`|π| = sqrt(n)` for every curve, so a_corn/n is pinned to a ~n^(−1/2)
sliver.  All Eisenstein ideals are similar as lattices (class number
1, hexagonal), so no isotropic invariant can separate them; the
information lives in the angle and only becomes visible under the
anisotropic weighting `diag(τ,1)`.

**Scope limit (Thread 21b).**  The "C2 ceiling" max{ν̂ : curve
recovers on all seeds} is *not* a constant — it depends on the
operating point `(eff, m)`:

| eff | 0.02 | 0.05 | 0.0993 | 0.15 | 0.25 |
|---|---|---|---|---|---|
| ceiling, m=12 | 1.0364 | 1.0151 | 0.6451 | 0.6017 | — (C2 empty) |
| ceiling, m=20 | 1.0364 | 1.0213 | 0.8674 | 0.7743 | 0.6690 |

So ν̂ is a *within-operating-point* separator (AUC 0.74–0.92
throughout), and any claim of the form "secp256k1 is / is not in the
easy tail" is meaningful only against an (eff, m)-matched ceiling.
Against matched ceilings secp256k1 is below at every point tested
except (eff = 0.0993, m = 12) — i.e. it is near the boundary, not
systematically inside or outside it.

**Standing caveat.** All of this is conditional on a *non-standard*
nonce generator `k = k1 + λ·k2 mod n` with `k1` small.  It is not an
attack on ECDSA as deployed, and nothing here bears on the main
theorem of `PAPER_STRUCTURAL_COMPLETENESS.md`.

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
