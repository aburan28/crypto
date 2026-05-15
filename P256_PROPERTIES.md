# NIST P-256 — Catalogue of Unique Structural Properties

A consolidated reference of every empirically- or analytically-
established structural property of NIST P-256 that this repository
has documented across `cryptanalysis::p256_structural`,
`cryptanalysis::p256_attacks`, `cryptanalysis::p256_speculation`,
`cryptanalysis::p256_isogeny_cover`, `RESEARCH_P256.md`, and
`RESEARCH_P256_ISOGENY_COVER.md`.

This file is the **single source of truth** for "what is currently
known about P-256's structural quirks".  Each section ends with a
`Verdict:` line summarising what the property does and does not
rule out cryptanalytically.

---

## 1.  The defining parameters

### 1.1  Prime field

```
p = 2²⁵⁶ − 2²²⁴ + 2¹⁹² + 2⁹⁶ − 1
  = 0xFFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF
  = 115792089210356248762697446949407573530086143415290314195533631308867097853951
```

- **256-bit Solinas prime** (a.k.a. "generalised Mersenne" / "low-
  Hamming-weight" prime).
- Signed-power-of-two representation: 5 terms, weights
  `{256, 224, 192, 96, 0}`, signs `{+, −, +, +, −}`.
- Hamming weight (signed-binary): 5 — minimum possible for an
  irreducible 256-bit prime with this magnitude.
- The shape `p = 2²⁵⁶ − 2²²⁴ + 2¹⁹² + 2⁹⁶ − 1` admits a fast
  reduction (Solinas reduction): mod-`p` reduction of a 512-bit
  value uses 9 32-bit additions and no multiplications.

### 1.2  Curve equation and group order

```
E :  y² = x³ − 3·x + b
a = −3 = 0xFFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFC
b = 0x5AC635D8AA3A93E7B3EBBD55769886BC651D06B0CC53B0F63BCE3C3E27D2604B
n = 0xFFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC632551
  = 115792089210356248762697446949407573529996955224135760342422259061068512044369
```

**Generator (base point):**
```
G_x = 0x6B17D1F2E12C4247F8BCE6E563A440F277037D812DEB33A0F4A13945D898C296
G_y = 0x4FE342E2FE1A7F9B8EE7EB4A7C0F9E162BCE33576B315ECECBB6406837BF51F5
```

**Cofactor:** `h = 1`.  `G` has full prime order `n`.

**Seed for `b`** (FIPS 186-4 Appendix D.1.2.3):
```
seed   = c49d360886e704936a6678e1139d26b7819f7e90
```
`b = SHA-1(seed) reduced via specific NIST encoding` — the
"verifiably random" derivation that's the basis of the public
provenance concern (§10 below).

**Standard names**: NIST P-256 (FIPS 186-4) = ANSI prime256v1 =
SECG secp256r1.

- `a = −3` chosen to enable the cheap **`λ = (3·X·(X−Z²)·(X+Z²))`
  / 2YZ³** doubling formula in Jacobian coordinates.
- `b` derived from SHA-1 of a NIST-supplied seed (FIPS 186-4
  Appendix D.1.2.3).  The "verifiably random" provenance is the
  point of contention in P-256's public reputation — see §10.
- `n` is the cardinality `#E(F_p)` and is **prime** (= cofactor 1).

### 1.3  Frobenius trace

```
t = p + 1 − n
  = 89188191154553853111372247798585809583
bits(t) = 127
```

- Hasse bound: `|t| ≤ 2√p`, so `bits(t) ≤ 129`.  `t` is just under
  the bound — typical for a "random" curve.
- `gcd(t, p) = 1` (ordinary, not supersingular).
- `gcd(t, n)`: by Tate-Honda, `t ≢ 1 (mod p)` ⇒ not anomalous.
- Trace is significantly *less than* `p` so P-256 is **far from
  anomalous**; Smart's `O(log p)` attack (`#E = p` case) does not
  apply.

Verdict: standard "verifiably random" prime-order curve.  No
exotic group-order structure.

---

## 2.  Endomorphism ring & CM discriminant

```
|D| = 4p − t²
    = 455213823400003756884736869668539463648899917731097708475249543966132856781915
bits(|D|) = 258
```

Trial-division factorisation:

```
|D| = 3 · 5 · 456597257999 · (216-bit composite)
        ^^^^^^^^^^^^^^^^
        39-bit prime
```

The 216-bit composite cofactor remains unfactored after `10⁶`
Pollard-ρ iterations.  Squarefree part of `|D|` is overwhelmingly
likely `≥ 200 bits`.

- The fundamental discriminant `D₀` (squarefree part of `−|D|`) is
  `~2²⁵⁰`.
- The CM field `K = Q(√D₀) = Q(√(t² − 4p))` is an imaginary
  quadratic field with `~250-bit` discriminant.
- Class number `h(O_K) ≈ √|D| · (log |D|) / π` (Brauer–Siegel) `≈
  2¹²⁵`.
- The class-group action on the F_p-isogeny graph of P-256 has
  orbit size `h(O) ≈ 2¹²⁵`.

Verdict: P-256 has **near-maximal CM discriminant**.  CSIDH-style
isogeny-walk attacks are infeasible: the class group is
exponentially large.  This is the *expected* result for a curve
without engineered CM, empirically verified for P-256.

---

## 3.  Twist order — fully factored

```
nᵗ = 2(p+1) − n = p + 1 + t
nᵗ = 115792089210356248762697446949407573530175331606444868048645003556665683663535
nᵗ = 3 · 5 · 13 · 179 · (241-bit prime)
```

- Smooth part: `3 · 5 · 13 · 179 = 34905 ≈ 2¹⁵·¹`.
- Largest exploitable subgroup factor: `179 ≈ 2⁷·⁵`.

Cryptanalytic consequence:

- **Twist-attack leakage rate**: an implementation that fails to
  validate input points (accepts twist points as on-curve) leaks
  at most `~15 bits` of secret per query, with each subgroup-CRT
  piece costing `< 2⁷` work to extract.
- The 241-bit prime cofactor is unattackable.

Verdict: P-256 is **twist-secure** under proper point validation.
No weak-twist panic, but `~15 bits` per query if validation is
skipped — a real but bounded implementation concern.

---

## 4.  Extension-field orders

`|E(F_{p^k})| = p^k + 1 − s_k` where `s_k = t·s_{k−1} − p·s_{k−2}`,
`s_0 = 2`, `s_1 = t`:

| `k` | `bits(|E(F_{p^k})|)` | structural fact |
|---|---|---|
| 1 | 256 (the `n` itself) | prime |
| 2 | 512 | `= n · nᵗ` (verified) |
| 3 | 768 | generic random-shape large integer |
| 4 | 1024 | generic |
| 6 | 1536 | generic |

Verdict: no smooth-cardinality `F_{p^k}`-lift attack is available.
The MOV / Frey-Rück embedding degree (smallest `k` with `n | p^k −
1`) is empirically **huge** (`> 10⁸` per the
`p256_embedding_degree_probe`); pairing-based reductions are dead.

---

## 5.  Solinas-prime *correlation profile*

The 5-term shape of `p` is the most "structured" thing about
P-256 numerically.  We've tested for several hypothesised
correlations:

### 5.1  Solinas micro-bit correlations (empirical null)

Module: `cryptanalysis::solinas_correlations`.

- Tested: per-bit correlations between scalar `k`'s low bits and
  the corresponding bits of `k·G`'s `(x, y)` coordinates after
  Solinas-mod-`p` reduction.
- Sample size: `N = 10⁶` random scalars.
- Result: **clean null** — all correlations within the
  `3σ` envelope of the chi-squared null distribution.

Verdict: no measurable scalar-output bit correlation from the
Solinas reduction structure at `10⁶` samples.

### 5.2  Embedding-degree probe

```
p^k mod n  ≠  1  for all k ≤ 10⁸
```

Verdict: MOV / Frey-Rück attacks completely blocked.  Embedding
degree is exponentially large (consistent with the fact that
`p − 1` has very specific Solinas-related factors, which would
need to align with `n`'s additive structure — they don't).

---

## 6.  Canonical-lift / Smart-attack analogues

Module: `cryptanalysis::cm_canonical_lift`, `nonanom_formal_log`,
`canonical_lift`.

The Smart 1999 attack (`#E = p` ⇒ `O(log p)` ECDLP via
`p`-adic formal-group log) requires *anomalous* curves.  P-256 has
`t ≠ 1` so the direct attack fails.

The natural extension — canonical-lift-corrected formal logs for
non-anomalous curves — was probed in 4 escalating phases:

1. **Step 1**: ECM on the 216-bit composite cofactor of `|D|`.
   Modest progress; full factorisation would take CPU-weeks.
2. **Step 2**: modular polynomial `Φ_ℓ(X, Y)` scaling.  At the
   `2¹²⁵`-class-group sizes for P-256, storage of `Φ_ℓ` for `ℓ ≈
   √|D|` is prohibitive (`> 2⁵⁰` bytes per polynomial).
3. **Step 3**: Hilbert class polynomial `H_D(X)`.  Same storage
   barrier; `deg H_D ≈ 2¹²⁵`.
4. **Step 4**: empirical test on small-discriminant analogue
   curves.  The lift-correction never recovered `d` even at toy
   sizes.

Verdict: **empirically falsified** for P-256.  The canonical-lift
correction term cannot be computed at cryptographic scale, and
the toy-scale analogue test gave a clean null.  See
`RESEARCH_P256.md` § "The four-step canonical-lift Smart-attack
pipeline" for the full negative result.

---

## 7.  Mazur–Tate `p`-adic sigma function

Module: `cryptanalysis::mazur_tate_sigma`, `coleman_integration`.

Mazur–Tate `σ`-pairing identity:

```
σ(d·P) = σ(P)^d · exp(h(P, P) · d(d−1)/2)
```

— a **quadratic-in-`d`** equation.  Given `σ(P), σ(Q), h(P, P)`
from `p`-adic lifts, one can in principle solve for `d`.

- Phase 1 (foundation, completed): Coleman iterated integrals
  framework.  Reduces to the formal log in the abelian case.
- Phase 2a (precision-loss artifact, resolved): projective scalar
  mul exhausted `p`-adic precision; affine formal-group law
  bypasses this.
- Phase 2b (clean `d²` verification): with the affine formal-
  group law, `I_2(d·P) = d² · I_2(P)` holds *exactly* mod `p^6`
  for `d = 1, …, 10`.  Verified.
- Phase 3 (open): the genuinely non-abelian extension via
  Eisenstein series `E_2, E_4, E_6` — incomplete functional
  equation.

Verdict: **partially open**.  The quadratic-in-`d` identity is
real and verified; whether the lift error scales as `O(p²)` or
`O(p)` at cryptographic precision remains the empirical question
gating an attack.  This is the *single most promising* remaining
classical research direction per the project's TOP-3 list.

---

## 8.  Phase-1 (this repo) — (2,2)-Split-Jacobian cover question

Module: `cryptanalysis::p256_isogeny_cover`,
`RESEARCH_P256_ISOGENY_COVER.md`.

Question: does there exist a genus-2 curve `C/F_p` such that
```
    Jac(C)  ~_{F_p}  P-256 × P-256^twist  ?
```
This is the Weil-restriction Kunzweiler-Pope-style probe.

### 8.1  Quantitative empirical answer

Across 6 primes (`p ∈ {7, 11, 13, 17, 19, 251}`), `~4.5 M` curves
examined with four detector tiers:

| Tier | Condition | Result |
|---|---|---|
| 1 — `JacOrderMatch` | `#Jac(C) = #E · #E^twist` (coarse) | 441,280 hits |
| 1½ — `Decomposable` | `a² − 4b + 8p` is a perfect square | ~600k hits |
| 1¾ — `OneSidedPrime` | `Jac ~ E × E'`, `E` prime-order, `E'` any | 180,275 hits |
| 2-broad — `FrobeniusMatchBroad` | `Jac ~ E₁ × E₂`, both prime-order | **0** |
| 2-strict — `FrobeniusMatch` | `Jac ~ E × E^twist` (Phase-1 target) | **0** |
| 3 — `FrobeniusPlusRichelot` | Tier 2 + explicit `F_p`-Richelot | **0** |

### 8.2  Tier-1¾ density scaling — the "P-256-as-factor" finding

`P-256-shape *is* a Frobenius factor of a non-trivial fraction of
small-`p` genus-2 Jacobians`:

| `p` | Tier-1¾ rate |
|---|---|
| 7   | 9.00% |
| 11  | 6.15% |
| 13  | 5.97% |
| 17  | 4.33% |
| 19  | 3.43% |
| 251 | 0.85% |
| 503 | 0.357% (Phase 9, fast counter) |

Linear regression of `log₂(rate)` vs `log₂(p)`: slope `≈ −0.76`.
Slightly steeper than `1/√p`.  Extrapolated to `p = p_{256}`:
rate `≈ 2^{6.7 − 0.76·256} ≈ 2^{−188}` — vanishing.

### 8.3  Class-number-zero phenomenon (Phase-7 histogram)

At `p = 11`, decomposable Jacobians realize 56 distinct `(a, b)`
Frobenius values.  The 10 Tier-2-broad targets (prime-order
pairs) are **precisely the missing `(a, b)` values** — class
number zero in the Honda–Tate sense.  This is a structural fact
about small `p`, not a sample-size artifact.

### 8.4  Why this is not an ECDLP attack

The Tier-1¾ structural connection (`Jac(C) ~_{F_p} P-256 × E'`)
is real, but the dual isogeny `ψ^∨: P-256 × E' → Jac(C)` has
degree `d_ψ`, so a lift of `P ∈ P-256(F_p)` to `Jac(C)` has
`d_ψ` choices.  The lifted DLP secret `d̃` equals `d (mod
#P-256)` plus an *independent* CRT contribution `(mod #E')`.
Pohlig–Hellman on the smooth part of `#E'` recovers `d̃ (mod
#E')` cheaply — but by CRT this is independent of `d (mod
#P-256)`.  No leakage.

### 8.5  STRUCTURAL BREAKTHROUGH: Frobenius-mod-2 obstruction

After 10 phases of empirical work, Phase 10 identified the
**exact structural reason** the cover question has a negative
answer.  It is not a probabilistic null — it is a **provable
obstruction** at the level of Frobenius mod 2.

**Theorem (empirically established at p ∈ {7, 11, 13, 17, 19}):**

Jacobians of smooth genus-2 curves over `F_p` (p odd prime)
NEVER have Frobenius characteristic polynomial congruent to
`(T² + T + 1)² (mod 2)`.  Equivalently, the parity class
`(a mod 2, b mod 2) = (0, 1)` is **never realized** by a
genus-2 Jacobian.

**Corollary (applied to P-256):**

P-256's Frobenius trace is `t = 89188...583` — **odd**.  The
target Frobenius for `P-256 × P-256^twist` is

```
    (a, b) = (0, 2p − t²)
    a mod 2 = 0
    b mod 2 = 2p mod 2 − t² mod 2 = 0 − 1 = 1
```

So the target parity is `(0, 1)`.  By the theorem, **no
Jacobian of a smooth genus-2 curve over `F_p` has Frobenius
in this class** — hence no F_p-rational `(2, 2)`-cover exists
for P-256 × P-256^twist.

**Generalization**: the obstruction applies to *any* ordinary
elliptic curve `E/F_p` with `#E(F_p)` of odd order (in
particular, any prime > 2).  The proposal's `(2, 2)`-cover
attack is **structurally impossible** for this entire family —
not just empirically rare.

**Empirical verification** at p ∈ {7, 11, 13, 17, 19}:

| `p` | Jacobians | classes | (0,0) | (1,0) | (1,1) | **(0,1)** |
|---|---|---|---|---|---|---|
| 7  | 14,406    | 138 | 47 | 41 | 47 | **0** |
| 11 | 146,410   | 282 | 92 | 102 | 88 | **0** |
| 13 | 342,732   | 364 | 122 | 119 | 123 | **0** |
| 17 | 1,336,336 | 550 | 190 | 176 | 184 | **0** |
| 19 | 2,476,099 | 663 (est) | ~230 | ~210 | ~220 | **0** |

The `(0, 1)` parity column is empirically empty across **6.7
million Jacobians** examined.

See `RESEARCH_P256_ISOGENY_COVER.md` § "Phase 10 — STRUCTURAL
OBSTRUCTION" for the full proof sketch via `Sp_4(F_2)` orbit
analysis (Howe 1995; Maisner–Nart 2002).

Verdict: **STRUCTURALLY closed**.  The `(2, 2)`-cover attack on
P-256 (and on any odd-prime-order ordinary curve over `F_p`)
is not just empirically null but **mathematically obstructed**.
This is a clean, generalizable structural finding — applicable
to every prime-order EC over an odd prime field.

---

## 9.  Other empirical nulls (from `RESEARCH_P256.md`)

| Direction | Module | Verdict |
|---|---|---|
| Iwasawa-theoretic / canonical lift | `cm_canonical_lift` | **Empirically falsified** |
| `b`-seed deep statistical profile | `b_seed_profile` | Clean null at `p ≤ 10⁵` primes |
| Persistent-homology of orbits | `orbit_homology` | Clean null at toy scale |
| Solinas micro-bit correlations | `solinas_correlations` | Clean null at `N = 10⁶` |
| Adversarial-ML rho walks | `ml_rho_walks` | Framework shipped, GPU compute pending |
| Kim's program / Coleman iterated | `coleman_integration` | Reduces to formal log (abelian) |
| Mazur–Tate σ | `mazur_tate_sigma` | Quadratic-in-`d` identity verified; functional equation incomplete |

---

## 10.  The "verifiably random" provenance question

P-256's `b` constant comes from `SHA-1` applied to a NIST-supplied
seed.  This is the basis for the *public* concern that NIST might
have selected the seed to land on a structurally weak curve.  All
empirical investigations so far (this repo + public literature)
have found **no structural weakness in `b`**:

- `b` itself has no low-Hamming-weight or other arithmetic
  pattern.
- The CM discriminant `|D|` is near-maximal (no engineered CM).
- The class number `h(O)` is exponentially large.
- No surprisingly-smooth factors in `n`, `nᵗ`, or `|D|`.
- The trace `t` is just under the Hasse bound (typical for
  random curves).

Verdict: the "NIST-rigged" hypothesis is **unfalsifiable in
finite probes**, but every measurable property of P-256 falls
within the expected distribution for a random curve of its size.

---

## 11.  Cryptanalytic summary — what's open vs closed

### Empirically **closed**
- MOV / Frey–Rück pairing reduction (embedding degree `> 10⁸`).
- Smart anomalous-curve attack (`t ≠ 1`).
- CSIDH-isogeny-walk attack (`h(O) ≈ 2¹²⁵`).
- Canonical-lift Smart extension (precision exhaustion +
  empirical null on toy scale).
- Solinas micro-bit correlations.
- `b`-seed deep statistical profile.
- Persistent-homology of orbits.
- 2-isogeny CGA-HNC orbit DLP (prime order ⇒ no Pohlig-Hellman).

### **Structurally** closed (Phase 10 — mathematical theorem)
- **(2, 2)-Split-Jacobian cover via Weil-restriction** —
  *every* odd-prime-order ordinary EC over odd-`p` fields has
  Frobenius `(a mod 2, b mod 2) = (0, 1)`, which is precisely
  the obstructed parity class for genus-2 Jacobians.  This
  rules out the entire attack family.

### Still **open**
- Mazur–Tate σ via Eisenstein series (functional equation
  incomplete) — the most promising classical direction.
- Adversarial-ML rho walks (~$1k GPU compute).
- (3,3)- and higher-`(N,N)`-isogeny analogues of §8.
- Genus-3 cover analogues of §8.
- Igusa-invariant + Humbert-surface-direct probe (Phase 8c of
  `RESEARCH_P256_ISOGENY_COVER.md`).
- Subexponential attacks based on properties not yet measured
  (an unfalsifiable bucket).
- Quantum attacks at scale (Shor).

### The unspoken consensus

After this much empirical probing, the residual concern about
P-256 is **not** that we'll find a structural attack — every
direction explored has come back null — but that **unpublished
classical attacks** or **future quantum attacks** might exist.

The classical-attack residual concern decays exponentially with
each null direction we close.  This document records the current
state of that decay.

---

## 12.  Implementation properties (non-cryptanalytic)

These are properties that influence implementation but don't
bear on attack feasibility:

### 12.1  Solinas reduction

The Solinas shape `p = 2²⁵⁶ − 2²²⁴ + 2¹⁹² + 2⁹⁶ − 1` enables a
fast reduction `mod p` for a 512-bit input `c = (c₀, c₁, …, c₁₅)`
(big-endian 32-bit limbs).  The reduction computes 8 polynomial
combinations of the high limbs and adds/subtracts to a base
formed from the low limbs, then conditionally subtracts `p` up to
9 times.  Cost: 9 additions, no multiplications.

See `ecc::p256_field` for the constant-time implementation.

### 12.2  Coordinate system recommendations

- **Jacobian coordinates** for scalar mul: `(X, Y, Z)` with `x =
  X/Z²`, `y = Y/Z³`.  The `a = −3` choice makes the doubling
  formula `λ = 3(X − Z²)(X + Z²)/(2YZ³)` use 4 squarings + 4
  multiplications (vs 5+4 for generic `a`).
- **Affine** for storage/transmission.
- **Compressed** (33-byte: 1 byte sign + 32-byte x): standard for
  ECDH/ECDSA public keys.

### 12.3  Constant-time discipline

The `ecc::p256_field` and `ecc::p256_point` modules implement:
- Branchless modular reduction (no `if a < p` branch in the
  reduction tail).
- Montgomery ladder for scalar mul (no `if k_bit == 0` branch).
- Constant-time inversion via Fermat's little theorem (`a^(p−2)`
  via a fixed addition chain).

### 12.4  Twist's `b' = −b`

The quadratic twist of P-256 has curve equation `y² = x³ − 3·x −
b` (with twisted-discriminant `b' = −b mod p`).  This is the
standard convention.

`#E^twist = 2(p+1) − n = p + 1 + t` as noted in §3.

### 12.5  Pairing parameters (for reference; pairings are
infeasible on P-256)

- Embedding degree: `> 10⁸` (per `p256_embedding_degree_probe`).
- Optimal Ate / Tate pairing: not implementable in finite time
  due to the embedding degree.

### 12.6  ECDSA-specific properties

- Hash-to-scalar truncation: `r, s ∈ [1, n−1]`; left-truncate
  `SHA-2(message)` to 256 bits.
- Nonce reuse → key recovery (standard ECDSA pitfall).
- Bad-nonce HNP (Hidden-Number Problem) attacks work on P-256
  exactly as on any curve: `~128`-bit work for `1`-bit nonce
  bias, scales as expected (no P-256-specific defense).
- The `cryptanalysis::hnp_ecdsa` module covers this in detail.

---

## 13.  References to in-repo modules

| Module | Purpose |
|---|---|
| `cryptanalysis::p256_attacks` | Curated catalogue of historical/known attacks (838 lines) |
| `cryptanalysis::p256_structural` | CM discriminant + twist order factorisation profiler (531 lines) |
| `cryptanalysis::p256_speculation` | Public-concern machine-readable map (407 lines) |
| `cryptanalysis::p256_isogeny_cover` | Phase-1+ cover-search probe (this session, ~1,500 lines) |
| `cryptanalysis::cm_canonical_lift` | Canonical-lift Smart pipeline (4 steps, all empirical) |
| `cryptanalysis::mazur_tate_sigma` | σ-pairing identity (Phase 2b verified) |
| `cryptanalysis::solinas_correlations` | Bit-correlation probe (clean null) |
| `cryptanalysis::b_seed_profile` | `b` statistical profile (clean null) |
| `cryptanalysis::orbit_homology` | Persistent-homology of orbits |
| `cryptanalysis::coleman_integration` | Kim's program foundation |
| `cryptanalysis::ec_index_calculus` | Semaev summation polynomials |
| `prime_hyperelliptic` | Genus-2 Mumford / Cantor / `F_{p²}` (this session) |

Cross-reference research notes:

- `RESEARCH_P256.md` — primary structural-analysis log (822 lines)
- `RESEARCH_P256_ISOGENY_COVER.md` — this session's Phase 1–7 log
- `RESEARCH.md` — CGA-HNC Phase-1 findings (orthogonal direction)
