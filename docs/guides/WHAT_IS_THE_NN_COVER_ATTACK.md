# What is the (N, N)-Cover Attack Family?

## A 1-paragraph summary

The **(N, N)-cover attack family** is a class of proposed attacks
against the elliptic curve discrete logarithm problem (ECDLP) that
work by finding a *higher-genus cover* — a smooth genus-2 curve
`C/F_p` whose Jacobian (a 2-dimensional generalized group) is
`F_p`-isogenous to `E × E^twist` (a product of `E` and its twist).
If such a cover exists, the ECDLP on `E` transports into the
*hyperelliptic curve discrete logarithm problem* (HCDLP) on
`Jac(C)`, which is potentially easier to solve via index-calculus
techniques.  The "(N, N)" refers to the kernel structure of the
isogeny: an `(N × N)`-shaped subgroup of size `N²`.

---

## The objects involved

### Elliptic curve `E`

`E` is a smooth curve of the form

```
    y² = x³ + a·x + b   (over F_p, characteristic ≠ 2, 3)
```

The set of `F_p`-rational points `E(F_p)` forms an abelian group
under the chord-tangent law.  The size `#E(F_p)` is between
`p + 1 − 2√p` and `p + 1 + 2√p` (Hasse).

### Twist `E^twist`

The **quadratic twist** of `E` is

```
    y² = x³ + a·d²·x + b·d³   for d a non-square in F_p
```

Geometrically `E` and `E^twist` are isomorphic over `F_{p²}` but
not over `F_p`.  They satisfy `#E + #E^twist = 2(p + 1)`.

### Frobenius characteristic polynomial

The Frobenius endomorphism `π: (x, y) ↦ (x^p, y^p)` satisfies a
characteristic polynomial:

- For `E`:   `π² − t·π + p = 0`   with trace `t = p + 1 − #E`
- For `E × E^twist`:   product of `E` and `E^twist` Frobenius polys:

```
    P_{E×E^twist}(T) = (T² − t·T + p)(T² + t·T + p)
                    = T⁴ − (t² − 2p)·T² + p²        ← a=0, b=2p−t²
```

### Cover Jacobian

`C: y² = f(x)` is a smooth genus-2 hyperelliptic curve with `deg f = 5
or 6`.  Its **Jacobian** `Jac(C)` is a 2-dimensional abelian variety
whose `F_p`-rational points form an abelian group of size
`#Jac(C)(F_p) ∈ [(√p − 1)⁴, (√p + 1)⁴]`.

The Frobenius char poly of `Jac(C)` is

```
    P_J(T) = T⁴ − a·T³ + b·T² − p·a·T + p²
```

By **Tate's isogeny theorem** (1966), two abelian varieties over
`F_p` are `F_p`-isogenous if and only if they have identical
Frobenius char polys.  So:

```
    Jac(C) ~_{F_p} E × E^twist
    ⇔
    (a, b) of Jac(C) equals (0, 2p − t²)
```

### `(N, N)`-isogeny

An isogeny between principally polarized abelian surfaces with
kernel `Z/N × Z/N` (size `N²`).  Special cases:

- `(2, 2)` — the classical **Richelot isogeny** (Richelot 1837),
  obtained by partitioning the 6 roots of `f` into 3 pairs.  The
  kernel has size 4.
- `(3, 3)`, `(5, 5)`, `(7, 7)`, …, `(11, 11)` — generalizations
  by Kunzweiler–Pope (2025) and others.  Kernel sizes 9, 25, 49,
  121 respectively.

Tate's theorem says the **existence** of an `F_p`-isogeny depends
only on Frobenius char polys, not on `N`.  So an `(N, N)` cover
exists for some `N` iff one exists for `N = 2`.

---

## The four steps of the attack

Given `P, Q = d·P ∈ E(F_p)`, the goal is to recover `d`:

### Step 1 — Embed into the product

`E(F_p)` embeds into `E × E^twist(F_p)` via `P ↦ (P, 0)`.  This
gives points `(P, 0)` and `(Q, 0)` in the product abelian variety.

### Step 2 — Lift to the cover via the inverse isogeny

If a cover Jacobian `Jac(C)` is `F_p`-isogenous to `E × E^twist`
via `ψ`, take the inverse isogeny `ψ⁻¹: E × E^twist → Jac(C)`.
The lifts

```
    D_P = ψ⁻¹(P, 0)
    D_Q = ψ⁻¹(Q, 0)
```

are divisors on `Jac(C)`, well-defined modulo `ker(ψ)` (size `N²`).

### Step 3 — Solve the HCDLP

`D_Q = d·D_P + (k, 0)` for some `k ∈ ker(ψ)`.  Modulo the kernel,

```
    d ≡ HCDLP(D_P, D_Q) on Jac(C)
```

For a genus-2 Jacobian, HCDLP can in some cases be solved faster
than ECDLP via Gaudry's index-calculus techniques.

### Step 4 — Recover `d`

`d` modulo `#E` is the answer.  The kernel ambiguity (`(N²)`
choices) is resolved by Pohlig–Hellman on the small kernel
subgroup.

---

## Why the attack works for binary curves but not prime fields

### Binary fields `F_{2^m}` (Teske 2006, GHS)

For curves over binary fields with **composite** `m`, there is a
subfield tower `F_2 ⊂ F_{2^d} ⊂ F_{2^m}`.  The original Teske/GHS
construction performs **Weil descent**: project the curve onto a
smaller-field curve via the subfield structure.  The resulting
cover `C/F_{2^d}` has genus `g = 2^{m_E − 1}` where `m_E` is the
"magic number" determined by the curve's coefficients.

For curves where `m_E` is small (e.g., 2, 3), the descent gives a
low-genus cover, and HCDLP on `Jac(C)` is feasible.  This is
**the original attack family** Teske proposed in 2006.

**Mitigation**: NIST chose `m` PRIME for `sect163k1` (`m = 163`),
`sect233k1` (`m = 233`), etc.  With `m` prime, there's no proper
subfield to descend through, neutralizing the attack.

### Prime fields `F_p` (this work)

For odd prime `p`, there is no proper subfield of `F_p`.  The
Teske/GHS construction cannot be applied directly.  The question
"does an analogous cover construction exist for prime fields?"
was the central question of the 21-phase P-256 investigation.

**The answer (Phase 10)**: For any prime-order ordinary `E/F_p`,
the target Frobenius char poly mod 2 is `(T² + T + 1)²` — a
conjugacy class in `Sp₄(F₂)` (the symplectic group on the
2-torsion of `Jac(C)`) that is **structurally inaccessible** to
smooth genus-2 Jacobians.  This is a known result of Howe (1995)
and Maisner–Nart (2002), computationally verified at scale by this
investigation.

So **no cover Jacobian `Jac(C)` exists over `F_p`** for the
prime-order-target case.  The attack family is structurally
impossible.

---

## Why the obstruction is specifically mod 2

The argument is a parity chain:

```
#E(F_p) odd prime > 2  ⇒  #E ≡ 1 mod 2
t = p + 1 − #E  ≡  (even) − (odd)  ≡  ODD  mod 2
t²  ≡  1  mod 2
b = 2p − t²  ≡  0 − 1  ≡  1  mod 2
(a, b) ≡ (0, 1)  mod 2

⇒ P_J(T) mod 2 = T⁴ + 0·T³ + 1·T² + 0·T + 1 = T⁴ + T² + 1 = (T² + T + 1)²
```

This `(T² + T + 1)² mod 2` is the unique parity class **not**
realized by genus-2 Jacobian Frobenius (Phase 10).

For curves with **even** `#E` (cofactor curves like Curve25519
with cofactor 8), the trace `t` is even, the target parity
becomes `(0, 0)`, and the obstruction doesn't apply.

---

## Exceptions: cofactor curves (Curve25519, Curve448, E-521)

For curves with even cofactor:
- `Curve25519` (cofactor 8): `#E = 8·q`, even
- `Curve448` (cofactor 4): `#E = 4·q`, even
- `E-521` (cofactor 4): same

For these, `t` is even, target parity `(0, 0)` — NOT obstructed by
Phase 10.

**Phase 18-21** showed that cover Jacobians **do exist** at toy
scale for these targets, with density roughly `p^{−3.4}` (Phase
19's strict measurement) or `p^{−0.7}` (Phase 18's broad
measurement).  At cryptographic scale (`p ≈ 2²⁵⁵`), the density
becomes `~2^{−858}` to `~2^{−181}` — vanishingly small, but
**non-zero structurally**.

So Curve25519/Curve448/E-521 are protected by **probabilistic
rarity** rather than the structural-impossibility theorem.

---

## Putting it all together

| Curve family | Protection mechanism | Mechanism type |
|---|---|---|
| Prime-order `F_p`-curves (NIST, secp, brainpool, FRP, SM2, GOST) | Phase 10 theorem | Structural impossibility |
| Cofactor `F_p`-curves (Curve25519, Curve448, E-521) | Phase 18 density argument | Asymptotic rarity |
| Binary curves (sect*) | NIST's prime-`m` design choice | Engineered immunity |

**Every deployed elliptic curve** is protected, just by **different
mechanisms**.  The investigation's contribution is to make these
mechanisms precise and verify them computationally.

---

## References

- **Teske, E.** (2006). *An Elliptic Curve Trapdoor System.*
  J. Cryptology. — The original binary-field attack.
- **Gaudry, P., Hess, F., Smart, N.** (2002). *Constructive and
  destructive facets of Weil descent on elliptic curves.*
  J. Cryptology. — The GHS construction.
- **Howe, E.** (1995). *Constructing distinct abelian varieties
  with isomorphic group structures.* — Sp₄(F₂) Jacobian-orbit
  result.
- **Maisner, D. & Nart, E.** (2002). *Abelian surfaces over
  finite fields as Jacobians.* — Characterizes which Frobenius
  classes contain Jacobians.
- **Tate, J.** (1966). *Endomorphisms of abelian varieties over
  finite fields.* — Tate's isogeny theorem.
- **Kunzweiler, S. & Pope, G.** (2025). *Efficient algorithms for
  the detection of `(N, N)`-splittings and endomorphisms.*
  J. Cryptol. — Modern `(N, N)`-isogeny machinery.

## In-repo references

- `cryptanalysis::p256_isogeny_cover` — Phases 1–21 of this
  investigation
- `cryptanalysis::ec_trapdoor` — Teske binary-curve implementation
- `cryptanalysis::ghs_descent` — GHS Weil descent
- `binary_ecc::*` — Binary-field EC infrastructure
- `prime_hyperelliptic::*` — Prime-field genus-2 Jacobian
  arithmetic (built during this investigation)
- `RESEARCH_P256_ISOGENY_COVER.md` — Full 21-phase research log
- `P256_RESEARCH_REPORT.pdf` — Comprehensive analytical PDF
- `MAINSTREAM_CURVE_AUDIT.pdf` — All 29 standardised curves
- `P256_VISUAL_REPORT.pdf` — 8-figure visual summary
