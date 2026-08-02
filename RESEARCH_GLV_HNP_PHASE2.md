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

### Two defects of the 2026-06-15 construction, and the fix (2026-08-02)

The implemented construction
(`secp256k1_cm_audit/glv_hnp_phase2_20bit.py:262`, `build_glv_lattice`)
carried two avoidable losses.  Both were measured and corrected in
`secp256k1_cm_audit/glv_hnp_phase2_projected.py` (Thread 23); the
numbers below are from
`glv_hnp_phase2_projected_output.txt`.

**(C) The target was never recentred.**  `k₁ᵢ` is uniform on `[0,K₁)`
and `k₂ᵢ` on `[0,K₂)`, so each planted coordinate had second moment
`n²/3` instead of `n²/12`.  Replacing `Aᵢ` by
`Aᵢ − λ·⌊K₂/2⌋ − ⌊K₁/2⌋` recentres both blocks and keeps the target a
*lattice* vector — no CVP step is needed.

**(P) The `d` column was charged to the norm.**  `d` is only defined
mod `n`, so `d²` should not appear in `‖v‖²`.  Quotienting `L` by its
kernel direction `ℤ·(n·S_D·e_m)` gives `π(L)` of rank `2m+1` with
`covol(π L) = covol(L)/(n·S_D)`.  The trivial vector `n·S_D·e_m` —
which 2026-07-29 exp T5 showed is *always* `λ₁` of the baseline
lattice — maps to `0`.  `d` is still recoverable without a
transformation matrix: the projected vector carries every `k₁ᵢ`,
`k₂ᵢ`, so `d = (k₁₀ + λ·k₂₀ − A₀)·B₀⁻¹ mod n`.

Closed form, with `S_K1·K₁ ≈ S_K2·K₂ ≈ n` and `S_D = 1`:

|            | `‖v_planted‖²`     | `covol`                  | dim    |
|------------|--------------------|--------------------------|--------|
| baseline   | `n²(2m/3 + 4/3)`   | `n^(3m+1)/(K₁K₂)^m`      | `2m+2` |
| proj+cent  | `n²(m/6 + 1)`      | `n^(3m)/(K₁K₂)^m`        | `2m+1` |

Since the Gaussian heuristic scales as `K₁^{−m/(2m+2)}`, a target
shorter by a factor `r` moves the critical `K₁` outward by
`r^{(2m+2)/m)}`.  For `m = 12`, `r = 1.763` predicts `3.4×`.

Measured (`m = 12`, 5 seeds, last `K₁` reaching 5/5):

| curve                       | λ\* | baseline | recentred | projected | proj+cent |
|-----------------------------|------|----------|-----------|-----------|-----------|
| `p=2557, n=2659`            | 0.340| `K₁=8`   | `K₁=16`   | `K₁=8`    | `K₁=16`   |
| `p=2677, n=2647`            | 0.070| `K₁=4`   | `K₁=16`   | `K₁=4`    | `K₁=16`   |

Cross-curve, 20 fresh 17-bit `j=0` GLV curves, `m = 12`:

| `eff = K₁K₂/n` | baseline           | proj+cent          |
|----------------|--------------------|--------------------|
| 0.15           | 2/20 curves, 13/100 trials | 14/20 curves, 93/100 trials |
| 0.25           | 0/20 curves,  6/100 trials |  3/20 curves, 53/100 trials |

Consequences for this document:

- The "K₁ wall" reported 2026-07-26 / 2026-07-29 is **not**
  information-theoretic; ~3–4× of it was a formulation artifact.
- The residual λ\* dependence of the wall (a factor ~3 between
  λ\*=0.07 and λ\*=0.34 in the baseline) **disappears** under
  recentring: both curves wall at `K₁ = 16`.  λ\* is not a parameter
  of the corrected attack.
- Projection alone changes nothing measurable (`V0 ≡ V2` in every
  cell above).  It is worth keeping for the cleaner geometry — it
  removes a vector shorter than the target — but the recovery
  condition remains BDD, not SVP: `sv/pv < 1` still holds in most
  cells even for proj+cent.  Making the target `λ₁` outright is
  still open.

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
