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
  **ANSWERED 2026-08-06 (§9): yes, after the d-free reformulation.**
  The reformulated attack succeeds exactly while the Gaussian-heuristic
  count `N_pred` of §9.4 is below 1, i.e. it saturates the
  information-theoretic bound of its own lattice; the §2 lattice did not.
- Can Phase 2 be **further specialised** to exploit secp256k1's
  CM-by-Z[ω] structure beyond just GLV?  E.g., if the
  endomorphism ω acts predictably on the bias distribution, a
  3D-aware HNP variant could exploit a third dimension.
- Does the **secp256k1 LLL-degeneracy** generalize to the GLV-
  aware lattice?  If yes, the entire Phase 2 needs a different
  base reduction algorithm.

## 9. The d-free BDD reformulation (Thread 23, 2026-08-06)

Script: `secp256k1_cm_audit/glv_hnp_phase2_bdd_reform.py`
Output: `secp256k1_cm_audit/glv_hnp_phase2_bdd_reform_output.txt`

### 9.1 Why the §2 lattice was the wrong lattice

In the §2 Kannan lattice (dim `2m+2`) the secret `d` gets its own column with
scale `S_D`.  Because `d` is uniform in `[0,n)`, the planted vector always
carries an `O(n)` coordinate, while the lattice contains the trivial vector
`n·S_D·e_m` of norm exactly `n·S_D`.  Since

```
||v_planted||^2 ~ n^2 (2m/3 + 4/3)   vs   ||n·S_D·e_m||^2 = n^2·S_D^2
```

the trivial vector is shorter for every `m >= 1`, and no choice of `S_D`
changes this — both scale linearly in `S_D`.  The planted vector is therefore
never `lambda_1` (empirically confirmed 2026-07-29, T5: `sv/pv` in
`[0.34, 0.61]`, 100% of the shortest vector's energy in the `d` column).

### 9.2 The reformulation

Eliminate `d` algebraically instead of giving it a column.  With
`c_i = B_i·B_0^{-1} mod n` and `w_i = k1_i + lam·k2_i = A_i + B_i·d (mod n)`:

```
w_i - c_i·w_0  =  A_i - c_i·A_0  =:  C_i   (mod n),    i = 1..m-1
```

which is `d`-free.  The solution set of the homogeneous system is a rank-`2m`
lattice `L0` with an explicit triangular basis (no HNF needed):

```
r_A    : u_0 = 1,   u_i = c_i               (i >= 1)
r_B    : v_0 = 1,   u_i = c_i·lam mod n     (i >= 1)
r_C,i  : v_i = 1,   u_i = -lam mod n        (i = 1..m-1)
r_D,i  : u_i = n                            (i = 1..m-1)
```

(columns `u` scaled by `S_K1`, columns `v` by `S_K2`).  The planted `(k1,k2)`
lies in the coset `t0 + L0` for `t0 = (0, C_1..C_{m-1}, 0_m)`, so recovery is a
genuine CVP; `d` is then read off coordinate 0 of the error via
`d = (k1_0 + lam·k2_0 - A_0)·B_0^{-1} mod n`.  Dimension drops `2m+2 -> 2m`
and the trivial vector is gone by construction.  Coset membership is verified
exactly in E1.

**Centering matters more than the dimension drop.**  `k1_i in [0,K1)` and
`k2_i in [0,K2)` are one-sided, so the raw error has a large mean.  Shifting
the CVP target by the box centre makes entries uniform on `[-n/2, n/2]`,
cutting `E[e_j^2]` from `n^2/3` to `n^2/12`.  Uncentered BDD is *worse* than
the §2 baseline; centered BDD is substantially better (§9.3).

### 9.3 Measured K1 wall (5 seeds each, successes out of 5)

`12-bit/2677`, `n=2647`, `lam*=0.070` — the curve recorded as a hard failure
on 2026-07-26 and as walled at `K1 ~ 4-6` on 2026-07-29:

| method | K1=2 | 3 | 4 | 6 | 8 | 12 | 16 | 24 | 32 |
|---|---|---|---|---|---|---|---|---|---|
| kannan-LLL (§2 baseline) | 5 | 5 | 5 | 1 | 0 | 0 | 0 | 0 | 0 |
| bdd-babai (uncentered) | 5 | 4 | 5 | 0 | 0 | 0 | 0 | 0 | 0 |
| bdd-babai-C (centered) | 5 | 5 | 4 | 5 | 4 | 5 | 1 | 1 | 0 |
| bdd-cvp-C (exact CVP) | 5 | 5 | 5 | 5 | 5 | 5 | 4 | 1 | 0 |

The wall moves from `K1 ~ 4-6` to `K1 ~ 16-24`: a factor of ~4.  On
`12-bit/2557` (`lam*=0.340`) it moves from `K1 ~ 8-12` to `K1 ~ 12-16`.

### 9.4 The new wall is information-theoretic

`||e_found||/||e_planted||` under *exact* CVP is `1.0000` (planted vector is
the closest vector) up to the wall, and drops below 1 past it — meaning a
strictly closer lattice point exists and no reduction algorithm can recover
`d`.  The transition is predicted by the Gaussian heuristic: with
`det(L0) = S_K1^m·S_K2^m·n^{m-1}` and centered `||e|| ~ n·sqrt(2m/12)`,

```
N_pred = vol(ball_{2m}(||e||)) / det(L0)
```

| curve | K1 | N_pred | exact-CVP outcome |
|---|---|---|---|
| 2557 (m=8) | 8 | 2.3e-03 | planted is closest |
| 2557 (m=8) | 16 | 5.9e-01 | LOST |
| 2677 (m=10) | 16 | 1.3e-01 | planted is closest |
| 2677 (m=10) | 32 | 1.4e+02 | LOST |

`N_pred` crosses 1 exactly where the empirical wall sits, on both curves.
The reformulated attack therefore **saturates the information-theoretic bound
of its own lattice**, which the §2 lattice did not.  At the boundary
(`2557`, `K1=16`, `N_pred=0.59`) exact CVP scores *lower* than Babai (1/5 vs
2/5) — precisely the signature of an information-theoretic rather than
algorithmic limit: exact CVP correctly returns a closer-but-wrong vector where
Babai occasionally stumbles onto the planted one.

### 9.5 Consequence for the paper's thesis

The bias budget `eff = K1·K2/n` at which GLV-HNP is viable improves by a
factor of roughly 2-4 (E4, 17-bit curves: `eff = 0.2` goes from 0/3-2/3 to
3/3).  This is a constant-factor improvement in the leakage requirement, not a
change in asymptotics — the attack still needs a per-signature nonce bias, so
it does not bear on the main theorem of `PAPER_STRUCTURAL_COMPLETENESS.md`.
It does mean prior Phase-2 negative results (2026-06-15 through 2026-07-29)
understated the attack by a factor of ~4 in `K1`, because they were all
measured on the `d`-column lattice.

## 10. References

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
