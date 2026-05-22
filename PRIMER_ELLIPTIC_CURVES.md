# A Visual Primer: CM, Jacobian, Twist, Isogeny, Rank, Embedding Degree

A no-prerequisites explainer for the six terms that show up everywhere
in elliptic-curve cryptography research.  Read this once and the
research notes in this repo (`RESEARCH.md`, `RESEARCH_P256.md`,
`RESEARCH_SECP256K1_CM.md`, etc.) become readable.

Each section: **What it is → Picture → What it looks like in secp256k1**.

---

## 0. The thing everything is about: an elliptic curve

An elliptic curve is the set of solutions `(x, y)` to an equation

```
        y² = x³ + a·x + b
```

over some "world of numbers" — the real numbers, the rationals, a
finite field `F_p`, the complex numbers, whatever.  Plus one extra
"point at infinity" called `O` that serves as the identity element.

**Over the real numbers**, the curve looks like one of these shapes:

```
         y² = x³ − x + 1                y² = x³ − 4x

           ╱⎺⎺⎺╲                          ╱⎺⎺⎺╲
          ╱     ╲                       ╱╲      ╲
   ───── ┃   o   ┃ ─────         ──── ╱   ╲      ┃ ─────
          ╲     ╱                ──── ╲   ╱      ┃
           ╲___╱                       ╲╱      ╱
       (one-piece)                  (two-piece — egg + curve)
```

The magic property: you can **add two points** to get a third point,
and this addition is associative and commutative.  Geometrically,
"P + Q" is the third intersection of the line through P and Q with
the curve, reflected across the x-axis.

```
                     • P+Q       (reflected back to the curve)
                    ╱
                   ╱
        ─────•────────•─────
              ╲ P    Q      ╲
               ╲              ╲
                ╲              ╲
                 ╲              ╲
                  •  (third intersection of line PQ with curve)
```

So **an elliptic curve is a *group***: it has an addition operation,
an identity (`O`), inverses (`−P = (x, −y)`), and associativity.
**Crypto uses this group**: scalar multiplication `k·P = P + P + … + P`
is the operation, and recovering `k` given `P` and `k·P` is the
**discrete log problem (ECDLP)** that secp256k1 / P-256 / Curve25519
rely on for security.

Now the six concepts.

---

## 1. CM — Complex Multiplication

**One-sentence definition:** An elliptic curve has *complex
multiplication* (CM) if it has more symmetries than the obvious ones.

### The obvious symmetries

For any elliptic curve, "multiply a point by an integer" is a
symmetry: `P ↦ 2P`, `P ↦ 3P`, `P ↦ −P`, etc.  These are the
endomorphisms `[n]` for `n ∈ Z`.  The set `End(E)` of all
endomorphisms always contains the integers `Z`.

For a *generic* elliptic curve, that's it — `End(E) = Z`.  No extra
symmetries.

### CM = a bonus endomorphism

A CM curve has an additional endomorphism that *isn't* multiplication
by an integer.  For **secp256k1**, that endomorphism is "multiply
by ω", where `ω` is a primitive cube root of unity (`ω³ = 1`, `ω ≠ 1`).
Concretely it's the map

```
        β(x, y) = (β · x, y)    for β = ω in F_p
```

Verify: cubing both sides, `(βx)³ + 7 = β³·x³ + 7 = 1·x³ + 7 = y²`.
So `(βx, y)` is also on the curve.  And it commutes with point
addition.  So `β` is an endomorphism that *isn't* "multiply by
some integer".

### Why "complex"?

Over the complex numbers, an elliptic curve is literally a torus
(a doughnut).  A torus is `C / Λ` for some lattice `Λ ⊂ C`.

```
   Generic lattice (most curves):       CM lattice (secp256k1):

        •   •   •   •                          •   •   •
                                            ╱  ╱  ╱  ╱
        •   •   •   •                       •   •   •
                                          ╱  ╱  ╱  ╱
        •   •   •   •                       •   •   •
                                          ╱  ╱  ╱  ╱
        •   •   •   •                       •   •   •

    Square-ish, no rotation                Hexagonal — has a
    symmetry except ±1                    60° rotation symmetry
                                          (rotation by ω)
```

The "extra symmetry" is literally **rotation of the complex plane**
by a non-real complex number.  For secp256k1's lattice, rotating
by `ω` (a 60° rotation, since `ω = e^{2πi/3}` is a primitive cube
root of unity) maps the lattice to itself — that's what makes
multiplication by `ω` a well-defined operation on the torus.

### Catalog of "small CM" curves

The only imaginary-quadratic fields `K` with class number 1 (= the
nicest possible CM) are

```
K = Q(√−1), Q(√−2), Q(√−3), Q(√−7), Q(√−11), Q(√−19),
    Q(√−43), Q(√−67), Q(√−163)
```

These nine fields are why "j = 0" (disc = −3) and "j = 1728" (disc =
−4) are special — they correspond to the smallest two of those nine.
**secp256k1 has `End = Z[ω]` where `ω² + ω + 1 = 0`**, coming from
`K = Q(√−3)`.

### What CM buys cryptographically

- **Speed**: The extra endomorphism gives a **2-dimensional scalar
  decomposition** (the GLV trick).  Instead of computing `k·P` for
  a 256-bit `k`, you write `k = k₁ + k₂·λ` with `|k₁|, |k₂| ≈ √k`
  and compute `k₁·P + k₂·β(P)` — half as many doublings.

- **Pollard-ρ speedup**: The factor-6 automorphism group gives a
  `√6 ≈ 2.45×` speedup of the standard ECDLP attack (Duursma–Gaudry–
  Morain 1999), reducing effective security from 128 bits to ~126.

- **No actual ECDLP break**: despite the rich CM structure,
  no one has found a way to exploit CM for sub-exponential ECDLP.
  See `RESEARCH_SECP256K1_CM.md` §9 for the audit.

---

## 2. Jacobian

**One-sentence definition:** The Jacobian `Jac(C)` of a curve `C` is
a higher-dimensional generalisation of the curve's group structure.

### For elliptic curves: it's just the curve

Elliptic curves are *genus 1*.  For these, `Jac(E) = E` itself — there
is no new object.  The Jacobian only becomes interesting at **higher
genus**.

### For genus 2 (hyperelliptic) curves

A genus-2 curve `C` has equation

```
        y² = f(x)        where deg(f) = 5 or 6
```

(So it's like an elliptic curve, but with a quintic / sextic instead
of a cubic on the right.)

**The points of `C` don't form a group.**  You can't naturally "add"
two points on `C`.  But pairs of points up to a certain equivalence
*do* form a group, and that group is `Jac(C)`.

```
   C:   genus-2 curve in the plane
        not a group

   Jac(C):    {pairs (P, Q) of points on C, up to equivalence}
             dimension 2
             IS a group
             |Jac(C)(F_p)| ≈ p²   (instead of ≈ p for elliptic)
```

Visually, `Jac(C)` for a genus-2 curve is a 4-real-dimensional
"surface" (an abelian surface, in geometric terms).  You can think
of it as a "two-torus" `C²/Λ` for a rank-4 lattice `Λ`.

### Why it shows up in the research

The `(N, N)`-cover idea in `RESEARCH_P256_ISOGENY_COVER.md` asks:
is there a genus-2 curve `C/F_p` such that `Jac(C)` is closely
related to `P-256 × P-256^twist`?  If yes, the projection map
`Jac(C) → P-256` might reveal structure invisible from looking at
P-256 alone.

```
        Jac(C)         ⟵ a 2-dimensional object
       ╱      ╲
      ╱        ╲
     ▼          ▼
   P-256    P-256^twist     ⟵ both 1-dimensional
```

This is the only "open" direction in the secp256k1 audit (slice 3).

---

## 3. Twist

**One-sentence definition:** A *twist* of `E/F_p` is another curve
`E'/F_p` that *looks like* `E` over a bigger field but *isn't* `E`
over `F_p`.

### The quadratic twist

Take `E: y² = x³ + 7` over `F_p`.  Pick any non-square `d ∈ F_p`.
The *quadratic twist* is

```
        E^t:  d · y² = x³ + 7

equivalent to

        E^t:  y² = x³ + 7·d³        (after substitution)
```

Over `F_p`, `E` and `E^t` are not isomorphic — they have different
group orders (`p + 1 − t` and `p + 1 + t` respectively, where `t` is
secp256k1's trace).

Over `F_{p²}` (the field where `√d` exists), they *become* isomorphic:
the map `(x, y) ↦ (x, y·√d)` sends one to the other.

```
   F_p:        E   ≄   E^t          (different group orders)
   F_{p²}:     E   ≅   E^t          (isomorphic via √d)
```

### Why j-invariant is the same

Two curves are isomorphic *over the algebraic closure* iff they have
the same **j-invariant**.  For `y² = x³ + ax + b`,

```
        j(E) = 1728 · (4a³) / (4a³ + 27b²)
```

`E` and `E^t` always have the same `j`, because j is invariant under
the twisting.  So **j-invariant classifies curves up to over-the-
algebraic-closure isomorphism**, and curves with the same j over
`F_p` are "twists of each other".

### Sextic twists (j = 0 only)

For most curves there are exactly **two** twists: `E` and its
quadratic twist `E^t`.

For `j = 0` (CM by `Z[ω]`) curves like secp256k1, there are **six**:

```
        E_d:  y² = x³ + d·b      for d ∈ F_p* / (F_p*)^6 ≅ Z/6Z
```

The six twists become isomorphic only over `F_{p⁶}`.  Geometric
picture:

```
   F_p:        E_1   E_2   E_3   E_4   E_5   E_6   (all distinct)
   F_{p²}:    {E_1, E_2}  {E_3, E_4}  {E_5, E_6}   (pairs merge)
   F_{p³}:    {E_1, E_5}  {E_3, ...}  ...          (triples merge)
   F_{p⁶}:                  all isomorphic
```

This is in `RESEARCH_SECP256K1_CM.md` §2: secp256k1 has six twists,
of which the original is the "twist #1" with trace `+L`, the
quadratic twist is "twist #2" with trace `−L`, and the other four
have traces involving `(L ± 3M)/2`.

### Why twists matter for security

If an implementation accepts an x-coordinate without checking it
lies on the actual secp256k1 curve, an attacker can submit `x` of a
point on a *twist* `E^t`.  The implementation does scalar multi-
plication anyway, leaking `k mod ord(P_twist)`.  If a twist has
small smooth-order factors, the secret leaks via Pohlig–Hellman.

For secp256k1, the worst twist (#4) has a 67-bit smooth piece — but
modern implementations validate curve membership, so this surface
is closed in practice.

---

## 4. Isogeny

**One-sentence definition:** An *isogeny* `φ: E → E'` is a map
between elliptic curves that preserves the group structure and has
a finite kernel.

### What it looks like

```
        E ──── φ ────▶ E'

        O    ↦    O                (identity preserved)
        P    ↦    φ(P)
        P+Q  ↦    φ(P) + φ(Q)      (homomorphism)

        ker(φ) = {P ∈ E : φ(P) = O}    finite subgroup of E
```

The **degree** of `φ` is `|ker(φ)|`.  A 2-isogeny has kernel of size
2 (one non-trivial element).  An ℓ-isogeny has kernel of size `ℓ`.

### Picture: identifying points

If `K = ker(φ)` is the kernel, then `E'` is geometrically `E / K`:
take `E`, identify any two points that differ by an element of `K`,
get `E'`.

```
        E:    •─•─•─•─•─•─•─•           (points along the group)
        K:    {O, T}      |T| = 2

        Identify P with P+T everywhere:
        E':   •─•─•─•                     (half as many points,
                                            still a group)
```

### What an isogeny does NOT do for ECDLP

If `Q = d · P` on `E`, then **the same `d` works after isogeny**:
`φ(Q) = d · φ(P)`.  So if you can solve ECDLP on `E'`, you can solve
it on `E`.

```
   E:    (P, Q = d·P)         ⟵ original ECDLP
         ↓ φ
   E':   (φ(P), φ(Q) = d·φ(P))  ⟵ same d, equivalent ECDLP
```

This is why "find an isogeny to an easier curve" is a tempting
attack — but it usually doesn't work, because **F_p-isogenous
curves all have the same group order `n`** (Tate's theorem).  If
the original `n` is a 256-bit prime, every neighbour has 256-bit
prime order too.

### The isogeny graph (the "volcano")

Fix a small prime `ℓ`.  Draw a graph: vertices are j-invariants,
edges are `ℓ`-isogenies.  Local structure looks like a volcano:

```
              ●         ●         ●           ⟵ CRATER
             ╱│╲       ╱│╲       ╱│╲             (curves with maximal
            ● ● ●     ● ● ●     ● ● ●             endomorphism ring)
           ╱ │ │ ╲   ╱ │ │ ╲   ╱ │ │ ╲
          ●  ●  ●  ●  ●  ●  ●  ●  ●  ●        ⟵ FLOOR LEVELS
                                                  (curves with smaller
                                                   endomorphism rings)
                  ...
```

- **Crater** = curves whose endomorphism ring is the *maximal*
  order in their CM field.  Edges *between* crater vertices are
  "horizontal" — they are the **`Cl(O)`-action** of the class group.
- **Floor levels** = curves with non-maximal endomorphism rings.
- **Edges going down** = "descending" isogenies, growing the
  conductor of the endomorphism ring.

For **secp256k1**: the crater has **exactly one vertex** (`j = 0`)
because the class number `h(−3) = 1`.  The "horizontal" CM-action
has nowhere to go — every horizontal `ℓ`-isogeny is a self-loop
back to `j = 0`.  This is why secp256k1 doesn't admit a CSIDH-style
class-group attack.

(`RESEARCH_SECP256K1_CM.md` §4 enumerates these self-loops empirically
for `ℓ ∈ {2, 3, 5, …, 47}`.)

### Different kinds of "isogeny"

- **F_p-rational isogeny**: defined by polynomials over `F_p`, sends
  F_p-points to F_p-points.  Group order preserved (Tate).
- **F_{p^k}-rational isogeny**: needs a bigger field to write down.
- **Horizontal vs descending**: same vs different endomorphism ring.

---

## 5. Rank

**One-sentence definition:** The *rank* of `E(K)` is the number of
"linearly independent infinite-order points" on `E` over the field `K`.

### For elliptic curves over the rationals Q

Over `Q`, the Mordell–Weil theorem says

```
        E(Q) ≅ Z^r ⊕ T
                │     │
                │     └── torsion (finite, ≤ 16 by Mazur 1977)
                │
                └── rank r ≥ 0
```

The rank `r` counts independent "free" generators of infinite order.

**Examples:**

```
   E: y² = x³ − 4x + 4         rank 1:  generator (0, 2)
                                        all of E(Q) is multiples of (0, 2)
                                        plus a small torsion part

   E: y² = x³ − x              rank 0:  only finitely many rational
                                        points — (0, 0), (1, 0), (−1, 0)
                                        and the point at infinity

   E: y² + xy = x³ − ...        rank 28: the current world record
                                        (Elkies 2006).  28 independent
                                        infinite-order rational points.
```

Computing rank is **hard**.  Beating-around-the-Birch-and-Swinnerton-
Dyer-conjecture hard.  The BSD conjecture is one of the Clay
$1M problems.

### For elliptic curves over a finite field F_p (what crypto uses)

**Over `F_p`, `E(F_p)` is finite** — there are no infinite-order
points.  The "rank" in the Mordell–Weil sense is **always 0**.

So when crypto papers talk about secp256k1, they don't usually
mention rank — the concept is for over `Q` (or other infinite
fields).  What crypto cares about is `|E(F_p)|` = the *order* of
the group, which for secp256k1 is the 256-bit prime `n`.

### Two unrelated uses of "rank" you might see

- **Rank of a lattice** (Z-module rank): for the period lattice
  of an elliptic curve over `C`, the rank is always 2 (it's a
  rank-2 lattice in `C`).  This is geometric, not arithmetic.
- **Rank of the endomorphism ring**: for ordinary curves over
  finite fields, `End(E)` is a Z-module of rank 2 (it's an order
  in an imaginary-quadratic field).  This is the "CM" picture.

Don't confuse these with Mordell–Weil rank.

---

## 6. Embedding degree

**One-sentence definition:** The embedding degree of `E(F_p)` with
respect to a subgroup of order `r` is the smallest `k` such that
`r` divides `p^k − 1`.

### Why this number matters

There's a thing called the **Weil pairing** (or Tate pairing) that
takes two points `P, Q ∈ E[r]` and returns an element of `F_{p^k}*`:

```
        e: E[r] × E[r] ──▶ F_{p^k}*
                            ▲
                            │
                       multiplicative group
                       of the extension field
```

This pairing is **bilinear**: `e(aP, bQ) = e(P, Q)^{ab}`.

So if you can solve a discrete log `Q = d · P` in `E(F_p)`, you can
*transport* it to a discrete log in `F_{p^k}*`:

```
        ECDLP on E(F_p)         DLP in F_{p^k}*
        find d s.t.  Q = d·P    find d s.t. e(P,Q) = e(P,P)^d
        cost ≈ √r               cost ≈ subexp(p^k)
                                       via NFS (number field sieve)
```

**If `k` is small**, the finite-field DLP is **much easier** than
the elliptic-curve one — this is the **MOV attack** (Menezes–
Okamoto–Vanstone 1993) and Frey–Rück.

### The trade-off

| Embedding degree `k`     | Cryptographic consequence                                    |
|--------------------------|--------------------------------------------------------------|
| `k = 1` or `k = 2`        | Curve is *broken* — pairing reduces ECDLP to a small DLP    |
| `k ∈ {6, 12, 24, 48}`     | **Pairing-friendly**: usable for pairing-based crypto       |
|                          | (BLS signatures, identity-based encryption, zk-SNARKs).      |
|                          | Curves like BLS12-381, BN254 live here.                     |
| `k ≈ √n` or larger        | "Generic" — MOV reduction gives nothing useful              |
| `k = O(n)`               | **Pairing-resistant** — no pairing attack possible.         |

### Where secp256k1 sits

For secp256k1, the embedding degree is `k ≈ n ≈ 2²⁵⁶`.

- The MOV reduction would map ECDLP to a DLP in `F_{p^{2^256}}*`.
- That field is incomprehensibly larger than the universe.
- Pairing-based attacks on secp256k1 are **completely infeasible**.

This was deliberate: secp256k1 was chosen to be pairing-*hostile*.
Conversely, BLS12-381 was chosen to be pairing-*friendly* (`k = 12`),
because it's used for things like Ethereum 2.0 BLS signatures where
pairings are a feature.

### Picture

```
   secp256k1 (pairing-hostile):

       ECDLP                        DLP in F_{p^{2^256}}*
        ▲                                  ▲
        │     pairing tunnel               │
        │     [length 2^256 — useless]     │
        │ ◀────────────────────────────────│
        │                                  │
       cost = 2^128                       cost ≈ subexp(p^{2^256}) = lol


   BLS12-381 (pairing-friendly):

       ECDLP                        DLP in F_{p^{12}}*
        ▲                                  ▲
        │     pairing tunnel               │
        │     [length 12 — short]          │
        │ ◀────────────────────────────────│
        │                                  │
       cost ≈ 2^128                       cost ≈ 2^126 (NFS)
                                          ⟵ targeted to match!
```

For BLS12-381 the parameters are *chosen* so that the DLP in
`F_{p^{12}}*` is approximately the same difficulty as the ECDLP.
This is what "pairing-friendly" means: efficient pairings, but
not *so* efficient that they break security.

---

> **Part II — supporting vocabulary.**  The next sections cover the
> standard ECC vocabulary the research notes use *without re-defining*:
> j-invariant, Frobenius, group order, endomorphism ring, ordinary
> vs. supersingular, torsion, discriminant, class group / Hilbert
> class polynomial / modular polynomial, point counting, the standard
> attacks, pairings, and the different curve forms (Weierstrass /
> Montgomery / Edwards).  Each is shorter than §0–§6: they're
> reference items, not deep dives.

---

## 7. j-invariant — the curve's "fingerprint"

**One-sentence definition:** `j(E)` is a single number in `F_p` that
classifies elliptic curves up to isomorphism over the algebraic
closure.

For `E: y² = x³ + ax + b`,

```
                       4 · a³
        j(E) = 1728 · ─────────────
                      4 a³ + 27 b²
```

Two curves over `F_p` have the same `j` ⇔ they become isomorphic
over `F_p`-bar.  So:

```
        ────────────────────────────────────────────
        E and E'  →  same  j   ⇔   they are twists
        E and E'  →  same  j   ⇔   E' = E (the case where
                                    the twist is trivial)
        ────────────────────────────────────────────
```

**Picture: j collapses twist families.**

```
   Curves over F_p:           j-invariant in F_p:

   ●────────────●─────╮           [j = 1234]
   │  E         E^t   │  ────▶    one number
   ●────────────●─────╯

   ●──●──●──●──●──●──╮            [j = 0]
   │ six sextic     │  ────▶      one number
   │ twists         │             (secp256k1 case)
   ●──●──●──●──●──●──╯
```

**Two special j-values** (the CM "corners"):

- `j = 0`: curves of the form `y² = x³ + b`.  Endomorphism ring
  contains `Z[ω]`.  secp256k1 lives here.
- `j = 1728`: curves of the form `y² = x³ + a·x`.  Endomorphism ring
  contains `Z[i]`.

These two are *the* "small CM" cases.

---

## 8. Frobenius and its trace `t`

**One-sentence definition:** The Frobenius endomorphism is the map
`π(x, y) = (x^p, y^p)`; it satisfies a quadratic `π² − t π + p = 0`
where `t` is the **trace of Frobenius**.

### Why this matters

Raising to the `p`-th power is the identity on `F_p` (Fermat's little
theorem) but acts non-trivially on `F_p`-bar.  Crucially, `π` sends
points on `E` to points on `E`, and it's a group homomorphism.  So
`π ∈ End(E)`.

Every Frobenius satisfies the equation

```
        π² − t · π + p = 0           (as endomorphisms)
```

The integer `t` is the **trace of Frobenius**.  It's the most
important invariant of `E/F_p` after `p` itself.

### Connection to the group order

```
        #E(F_p) = p + 1 − t
```

So knowing `t` is the same as knowing `|E(F_p)|`.

### Hasse bound

`t` can't be arbitrary.  The Hasse theorem says

```
        |t| ≤ 2 · √p
```

So `#E(F_p) ∈ [p + 1 − 2√p, p + 1 + 2√p]` — concentrated near `p`.

```
                Hasse interval
        ◀────────────────────────────▶
        p+1−2√p                  p+1+2√p
                ┃
                ● #E(F_p) lives somewhere in here
                ┃
                p + 1
```

### For secp256k1

`t = 432420386565659656852420866390673177327` (a 129-bit integer).
That `t` is positive and even-bit-sized while `p` is 256-bit is a
small structural quirk — see `RESEARCH_SECP256K1_CM.md` §1.

---

## 9. Group order `n`, cofactor `h`, and "the subgroup we use"

**One-sentence definition:** `#E(F_p) = h · r` where `r` is a large
prime (used for crypto) and `h` is a small **cofactor**.

A crypto curve has its group order factored as

```
        #E(F_p) = n = h · r
                  ▲    ▲   ▲
                  │    │   │
                  │    │   └── large prime (used for actual crypto)
                  │    └────── small cofactor (1, 4, 8, …)
                  └──────────── total order
```

We only ever do ECDLP-style work in the unique subgroup of order
`r`.  The cofactor `h` is "noise" we strip off.

**Why cofactor matters:**

- `h = 1` (secp256k1, P-256, P-384, P-521, …): no noise. Whole group
  is prime-order.  Cleanest possible setting.
- `h = 4` (Curve448): mild noise.  Implementations must clear the
  cofactor explicitly to avoid small-subgroup confusion.
- `h = 8` (Curve25519): same.  The `clamp-the-low-3-bits` step in
  X25519 is exactly the cofactor-clearing operation.
- `h` large: cofactor attacks become a thing.  **Avoid** — most
  modern curves have `h ≤ 8` for this reason.

For secp256k1: `h = 1`, `r = n` (the 256-bit prime).  Cleanest case.

---

## 10. Endomorphism ring `End(E)`

**One-sentence definition:** `End(E)` is the ring of all isogenies
from `E` to itself (with addition and composition).

`End(E)` always contains the integers `Z` (via the maps `[n]: P ↦
n·P`).  For elliptic curves over a finite field, there are exactly
two possible structures:

```
   Ordinary curve:                       Supersingular curve:

        Z ⊂ End(E)                            Z ⊂ End(E)
            │                                     │
            ├── rank 2 over Z                     ├── rank 4 over Z
            │                                     │
            └── order in imaginary                └── order in
                quadratic field K                     quaternion algebra
                K = Q(√D), D < 0                       (non-commutative!)

        e.g. End = Z[ω] for sec256k1          rare — j ∈ a small finite
                                                  set of values per p
```

### Ordinary vs supersingular: the dichotomy

Every elliptic curve over `F_p` is exactly one of:

- **Ordinary** (`gcd(t, p) = 1`): the common case.  `End(E)` is a
  rank-2 commutative ring.  All deployed curves (secp256k1, P-256,
  Curve25519, etc.) are ordinary.

- **Supersingular** (`p | t`, often just `t = 0`): rare and special.
  `End(E)` is a rank-4 non-commutative ring (an order in a
  quaternion algebra).  The supersingular j-invariants form a finite
  set of size `≈ p/12`.  Used in **isogeny-based post-quantum
  cryptography**: SIDH (broken in 2022 by Castryck–Decru), CSIDH,
  CSIKE.

### Picture of the supersingular world

```
   F_p-bar:    ●───●───●───●───●─── ...     (a finite set of
               │   │   │   │   │              j-invariants)
               └───┴───┴───┴───┴─── ...
                  (richly connected by isogenies)

   All curves in this set are isogenous to each other via
   isogenies of varying degree.  The graph is connected and
   non-trivially structured — basis for SIDH/CSIDH.
```

For everything in this repo (secp256k1, P-256, etc.) the curves
are **ordinary**.  The supersingular world is mentioned only when
contrasting (the (N,N)-cover detector in
`RESEARCH_P256_ISOGENY_COVER.md` was originally invented for the
supersingular case).

---

## 11. Torsion subgroup `E[n]`

**One-sentence definition:** `E[n] = {P ∈ E(F_p\text{-bar}) : n·P =
O}` — the points of order dividing `n`.

### Structure

For `n` coprime to `p`:

```
        E[n] ≅ Z/n × Z/n        (over F_p-bar)
```

So `E[n]` has exactly `n²` points (including `O`).  As a Z-module
it has **rank 2** — *this* is the "rank 2" that everyone mentions,
not the Mordell-Weil rank from §5.

### F_p-rational part

Over `F_p` (instead of `F_p`-bar), only a *subgroup* of `E[n]` is
rational:

```
        E(F_p)[n] ⊆ E[n]            usually proper subgroup
```

The Frobenius `π` acts on `E[n]` as a 2x2 matrix over `Z/n`; the
`F_p`-rational part is the fixed subspace.

### Why we care

- **Kernel of an `ℓ`-isogeny** is a cyclic order-`ℓ` subgroup of
  `E[ℓ]`.  Picking different kernels gives different isogenies.
- **Pairings** (§16 below) live on `E[n]`.
- **Hilbert / modular polynomial computations** (§13) implicitly
  use the structure of `E[ℓ]` as a `Z[π]`-module.

---

## 12. Discriminant

**One-sentence definition:** A number that summarises whether a
curve / order is "non-degenerate" and (in the CM case) records
which imaginary-quadratic field we're in.

### Curve discriminant

For `E: y² = x³ + ax + b`,

```
        Δ(E) = −16 · (4 a³ + 27 b²)
```

`E` is a smooth (= "actually an elliptic curve") iff `Δ(E) ≠ 0`.
If `Δ = 0` the cubic has a repeated root and the "curve" is a
node or cusp.

### Field discriminant (the more important one)

For `E/F_p` *ordinary*, `End(E) ⊗ Q = K` is an imaginary-quadratic
field `Q(√D)` with `D < 0`.  The **fundamental discriminant** of `K`:

- `D = −3`: `K = Q(√−3)` — secp256k1 lives here.
- `D = −4`: `K = Q(√−1)` (Gaussian integers).
- `D = −7, −8, −11, −19, −43, −67, −163`: the other "class number 1"
  fields.

### Order discriminant

Inside `K = Q(√D)`, there are *many* orders: `Z[√D]`, `Z[ℓ·ω]`,
etc.  Each has its own discriminant `D · f²` where `f` is the
**conductor**.

```
        maximal order:                      conductor 1, disc D
        Z + ℓ Z[ω]:                         conductor ℓ, disc D · ℓ²
        Z + ℓ² Z[ω]:                        conductor ℓ², disc D · ℓ⁴
```

For secp256k1 the endomorphism order is the maximal order `Z[ω]`,
disc `−3`, conductor `1`.  The volcano descent in `RESEARCH_SECP256K1_CM.md`
§4 walks to orders of disc `−3·ℓ²` (conductor `ℓ`).

---

## 13. Class group `Cl(O)`, Hilbert and modular polynomials

These three live together: they describe the structure of the set
of elliptic curves with given endomorphism ring, and how isogenies
connect them.

### Class group `Cl(O)` and class number `h(D)`

For an order `O` in an imaginary quadratic field, the **class group**
`Cl(O)` is a finite abelian group that measures "how far `O` is from
being a PID" (a principal-ideal domain).  The **class number** is
`h(D) = |Cl(O)|`.

```
        h(−3) = 1               ← secp256k1 / j=0 case
        h(−4) = 1               ← j=1728 case
        h(−23) = 3              ← smallest "interesting" h
        h(−47) = 5
        h(−71) = 7              ← used in the slice-2 demo
```

The class group **acts on the set of curves with `End = O`**: pick
an ideal class `[a] ∈ Cl(O)`, get an isogeny `E → [a] · E`.  This
action is free and transitive; the orbit has size `h(D)`.

### Hilbert class polynomial `H_D(X)`

```
        H_D(X) = ∏  (X − j(E_i))    over the h(D) curves E_i
                 i                    with End(E_i) = O_D
```

`H_D` has integer coefficients.  Its roots over `F_p` are the
j-invariants of curves with `End = O_D` defined over `F_p` (when
`p` splits in `O_D`).

We used `H_{−71}(X)` in the slice-2 demo in
`RESEARCH_SECP256K1_CM.md` §7: factored `H_{−71}` mod 107, got
seven j-invariants `{19, 30, 46, 57, 63, 64, 77}`.

### Modular polynomial `Φ_ℓ(X, Y)`

```
        Φ_ℓ(j(E), j(E')) = 0   ⇔   E, E' are connected
                                     by a cyclic ℓ-isogeny
```

So `Φ_ℓ` is the "ℓ-isogeny detector".  Setting `Φ_ℓ(0, Y) = 0` and
factoring over `F_p` gave us the secp256k1 isogeny-graph table in
`RESEARCH_SECP256K1_CM.md` §4.

`Φ_ℓ` has integer coefficients, symmetric in `X, Y`, of degree `ℓ + 1`
in each variable.  Coefficients grow fast — `Φ_2` has small (~10-
digit) coefficients; `Φ_100` has coefficients with thousands of
digits.

### How they fit together

```
   Cl(O)              ⟵ the algebra
   ▲
   │ acts on
   │
   {j-invariants}     ⟵ roots of H_D
   ▲
   │ connected by
   │
   ℓ-isogenies        ⟵ zeros of Φ_ℓ
```

---

## 14. Point counting: Schoof / SEA

**One-sentence definition:** Algorithms that compute `#E(F_p)`
efficiently — the engineering problem you face whenever someone
hands you a curve.

Computing `#E(F_p)` naively is `O(p)` — infeasible for crypto-sized
`p`.  Three algorithms:

- **Schoof 1985**: polynomial time `O(log⁸ p)`.  The breakthrough.
  Works by computing `t mod ℓ` for many small primes `ℓ` using the
  ℓ-th division polynomial, then CRT.
- **SEA (Schoof–Elkies–Atkin)**: practical improvement.  For
  primes `ℓ` where Φ_ℓ splits nicely (Elkies primes), uses an
  ℓ-isogeny instead of the full division polynomial.  Reduces
  per-`ℓ` work dramatically.
- **PARI's `ellcard` / `ellsea`**: production implementations.
  Called in the slice-2 PARI script in `secp256k1_cm_audit/`.

### Why this matters here

Every research note in this repo *assumes* you can ask "what is
`#E(F_p)` for this curve?" and get an answer in seconds.  That
capability is Schoof/SEA; without it, none of this would be
computationally feasible.

For curves we *design* (like CM curves built via the Hilbert class
polynomial), we don't even need Schoof — we know `t` already from
the Cornacchia decomposition `4p = a² + |D| b²`.  See
`RESEARCH_SECP256K1_CM.md` §1.

---

## 15. The standard ECDLP attacks

A quick catalog of the algorithms a curve has to survive to be
"secure", in roughly increasing specialness:

### Pollard ρ — the generic attack

Works on **any** group of order `r`.  Cost: `≈ √(π r / 2)` group
operations.  For secp256k1 (`r ≈ 2²⁵⁶`), that's `≈ 2¹²⁸` operations
— the standard "128-bit security" claim.

```
        ECDLP cost (ρ):  √r   ≈  √n
```

Speeds up by `√Aut(E)` if there are automorphisms — `√6 ≈ 2.45×`
for j=0 curves.

### Pohlig–Hellman — works on smooth orders

If `r = ∏ ℓᵢ^eᵢ` with each `ℓᵢ` small, decompose ECDLP into
sub-ECDLPs of size `√ℓᵢ` each.  Total cost: `Σ √(ℓᵢ^eᵢ)`.

```
   If r = 2¹⁰⁰: PH cost  ≈ √2 · 100 = trivial
   If r = (small prime)^k · large prime:
                PH cost  ≈ √large_prime    ⟵ bottleneck
```

This is why crypto curves use **prime-order subgroups** — denies PH
any traction.  The CGA-HNC programme is structured around looking
for PH-amenable factors of `n` across an isogeny orbit; see
`RESEARCH.md`.

### Smart's attack — anomalous curves

If `#E(F_p) = p` exactly (so `t = 1`), the curve is "anomalous" and
Smart 1999 gives a **polynomial-time** attack via the formal
logarithm / canonical lift.  This is why every curve standard
explicitly excludes anomalous curves.

secp256k1's `n = p + 1 − t` with `t ≈ 2¹²⁹` ≠ 1, so this doesn't
apply.

### MOV / Frey–Rück — pairing reduction

Already discussed in §6.  Reduces ECDLP to DLP in `F_{p^k}*`.  Useful
only when `k` is small.

### Index calculus on `F_{p^k}` (Diem, Semaev)

Sub-exponential on `F_{p^k}` for `k ≥ 3` (Diem 2011), but
heuristically slower than ρ for prime-field curves (`k = 1`).
"Prime-field ECDLP looks hard" is the working assumption.

### Weil descent / GHS

Specific to binary-field curves `E/F_{2^n}` with `n` composite.
Doesn't apply to prime-field curves like secp256k1.

### Summary

```
   For secp256k1:
        ρ:                √n   = 2¹²⁸           — best generic attack
        Pohlig-Hellman:   N/A (n is prime)
        Smart:            N/A (n ≠ p)
        MOV:              N/A (k ≈ n)
        Diem:             N/A (over F_p, not F_{p^k})
        ⇒ best known attack: 2¹²⁶ (ρ with √6 GLV speedup)
```

---

## 16. Pairings (Weil, Tate, Ate)

**One-sentence definition:** A *pairing* is a bilinear map
`e: E[n] × E[n] → μ_n` where `μ_n` is the group of `n`-th roots of
unity in some extension field.

```
        e:    E[n] × E[n]  ──▶  μ_n  ⊂  F_{p^k}*

        Bilinear:   e(a·P, b·Q) = e(P, Q)^{ab}
        Non-degenerate: there exist P, Q with e(P, Q) ≠ 1
```

### What pairings let you do

- **MOV / Frey-Rück reduction** (§6, §15): pull ECDLP into a
  finite-field DLP.  Attack — when `k` is small.
- **BLS signatures** (Boneh-Lynn-Shacham): short signatures whose
  verification is a pairing equation.  Used by Ethereum 2.0.
- **Identity-based encryption** (Boneh-Franklin 2001): public key
  *is* the user's identity (email address etc.); pairings make the
  decryption work.
- **zk-SNARKs**: many proof systems (Groth16, Plonk over BLS12-381)
  use pairings.

### The three named pairings

- **Weil pairing**: most symmetric.  `e(P, Q) · e(Q, P) = 1`.
  Slower in practice.
- **Tate pairing**: faster, less symmetric.  The actual one used in
  implementations.
- **Ate / R-ate / optimal Ate**: speed variants of Tate that exploit
  Frobenius eigenstructure.  BLS12-381 implementations use these.

In this repo: `crypto/src/bls12_381/` has a Tate-pairing
implementation (for BLS sig verification, not as an attack tool).

---

## 17. Other curve forms (Weierstrass, Montgomery, Edwards)

The same elliptic curve can be written in different ways.  Each form
has speed / safety / simplicity trade-offs.

### Short Weierstrass

```
        y² = x³ + a x + b
```

The standard form.  Used by NIST P-{256, 384, 521}, secp256k1,
brainpool.  Has 3 distinguished `2`-torsion points (`(α, 0)` for
each root of the cubic) when they're `F_p`-rational.

### Montgomery

```
        B · y² = x³ + A · x² + x
```

Used by **Curve25519** (ECDH only).  Designed by Bernstein for
constant-time x-coordinate-only scalar multiplication (the
"Montgomery ladder").  Same curve as Edwards form below, in different
coordinates.

```
   Standard ECDH with Montgomery + ladder:
       no y-coordinate ever computed
       no branching on key bits
       ⇒ side-channel resistant by construction
```

### (Twisted) Edwards

```
        Edwards:           x² + y²  = 1 + d · x² · y²
        twisted Edwards:   a x² + y² = 1 + d · x² · y²
```

Used by **Ed25519** (EdDSA signatures), **Ed448**.  Beautifully
*complete* addition law: no special case for "point at infinity",
no special case for `P + (−P) = O`.  This eliminates a whole class
of implementation bugs.

```
        Weierstrass addition law:
            if P == −Q: result is point at infinity
            else: chord-and-tangent
                ⇒ requires branching, side-channel risk

        Edwards addition law:
            (x₁,y₁) + (x₂,y₂) = (numerator/denom, …)
            same formula always works
                ⇒ no branching, no infinity special case
```

### Equivalences

A short-Weierstrass curve is birationally equivalent to a Montgomery
curve iff it has a single F_p-rational 2-torsion point.  Montgomery
is birationally equivalent to twisted-Edwards iff it has F_p-rational
points of order 4.  So `Curve25519` (Montgomery) and `Ed25519`
(Edwards) are the *same elliptic curve* in different costumes — chosen
to support both ECDH (Montgomery ladder) and EdDSA (complete Edwards
addition).

---

## Cheat-sheet (one-line summaries)

| Term                       | What it is                                                              | secp256k1 case                                  |
|----------------------------|-------------------------------------------------------------------------|-------------------------------------------------|
| **Curve**                  | `y² = x³ + ax + b`, a group under "chord-and-tangent" addition          | `y² = x³ + 7` over `F_p` (256-bit `p`)          |
| **CM**                     | Extra endomorphism beyond multiplication by integers                    | Multiplication by `ω` (cube root of 1)          |
| **Jacobian**               | Higher-genus generalisation of the curve-as-group                       | `Jac(E) = E` (genus 1, nothing new)             |
| **Twist**                  | A different `F_p`-curve with the same j-invariant                       | Six sextic twists, distinct over `F_p`          |
| **Isogeny**                | Group-preserving map with finite kernel                                 | F_p-isogenies preserve `n` exactly              |
| **Rank**                   | Number of independent infinite-order points (only for `E(Q)`-type)      | Doesn't apply to `F_p`-curves                   |
| **Embedding degree**       | Smallest `k` with `n \| p^k − 1`                                        | `k ≈ n ≈ 2²⁵⁶` — pairing-hostile                |
| **j-invariant**            | Single number in `F_p` that classifies `E` over `F_p`-bar               | `j = 0` (the "small CM disc=−3 corner")        |
| **Frobenius `π`**          | The `(x, y) ↦ (x^p, y^p)` endomorphism; satisfies `π² − tπ + p = 0`     | `π = a + bω` in `Z[ω]`, see §3 of audit         |
| **Trace `t`**              | The integer in the Frobenius equation; `#E(F_p) = p + 1 − t`            | `t = 4.32 × 10³⁸` (129-bit)                    |
| **Hasse bound**            | `\|t\| ≤ 2√p` — the only constraint on `t`                              | Easily satisfied                                |
| **Cofactor `h`**           | Small factor of `#E` left after extracting the large prime               | `h = 1` (whole group is prime-order)            |
| **Endomorphism ring End(E)**| All isogenies from `E` to itself — usually `Z`, sometimes more         | `End = Z[ω]` (rank 2)                           |
| **Ordinary**               | `gcd(t, p) = 1` — `End` is in an imaginary quadratic field              | Yes, ordinary                                   |
| **Supersingular**          | `p \| t` — `End` is in a quaternion algebra; used in SIDH/CSIDH         | No                                              |
| **Torsion `E[n]`**         | Points of order dividing `n`; `E[n] ≅ Z/n × Z/n` over `F_p`-bar         | n² points each                                  |
| **Discriminant** (curve)   | `Δ = −16(4a³ + 27b²)`; curve smooth ⇔ `Δ ≠ 0`                          | Non-zero (computable)                           |
| **Discriminant** (field)   | `disc(K)` for `K = End ⊗ Q`; classifies the CM field                    | `disc(Q(√−3)) = −3`                            |
| **Class group `Cl(O)`**    | Finite abelian group acting freely+transitively on the End=O orbit      | `Cl(Z[ω]) = trivial`, `h(−3) = 1`              |
| **Hilbert class poly `H_D`**| Polynomial whose roots are j-invariants with `End = O_D`               | `H_{−3}(X) = X` — single root `j=0`            |
| **Modular poly `Φ_ℓ`**     | `Φ_ℓ(j, j') = 0` ⇔ `E, E'` are `ℓ`-isogenous                            | Used in audit for `ℓ ≤ 47`                      |
| **Schoof / SEA**           | Polynomial-time `#E(F_p)` computation                                   | The reason audits run in seconds                |
| **Pollard ρ**              | Generic `√r` ECDLP attack                                               | `≈ 2¹²⁸` operations (with √6 speedup: `2¹²⁶`)  |
| **Pohlig-Hellman**         | Attack on smooth-order subgroups                                        | N/A: `n` is prime                               |
| **Smart's attack**         | Polynomial-time when `#E = p` (anomalous)                               | N/A: `n ≠ p`                                    |
| **MOV / Frey-Rück**        | Pairing reduction; effective only if embedding deg small                | N/A: `k ≈ n`                                    |
| **Pairing**                | Bilinear `e: E[n] × E[n] → F_{p^k}*`                                    | Exists; useless for attack (`k` huge)          |
| **Weierstrass form**       | `y² = x³ + ax + b` — the standard                                       | secp256k1's form                                |
| **Montgomery form**        | `B y² = x³ + A x² + x` — fast ladder, side-channel safe                | Curve25519 uses this                            |
| **Edwards form**           | `x² + y² = 1 + d x² y²` — complete addition, no special cases          | Ed25519 uses this                               |

## Where to read next (in this repo)

In rough order of mathematical depth:

1. **[crypto/RESEARCH_SECP256K1_CM.md](RESEARCH_SECP256K1_CM.md)** —
   the recent structural audit of secp256k1.  Now should be readable.
2. **[crypto/RESEARCH.md](RESEARCH.md)** — the CGA-HNC programme.
   Builds on the isogeny + class-group concepts above.
3. **[crypto/RESEARCH_P256.md](RESEARCH_P256.md)** — broader P-256
   cryptanalysis programme.  Touches on many of these concepts.
4. **[crypto/RESEARCH_P256_ISOGENY_COVER.md](RESEARCH_P256_ISOGENY_COVER.md)**
   — the `(N, N)`-cover direction.  Needs Jacobian + isogeny + twist
   all together.

## External references (textbook level)

- **Silverman, *The Arithmetic of Elliptic Curves*** — the standard
  textbook.  Chapters III (group law), V (over finite fields),
  III.10 (endomorphism rings / CM).
- **Washington, *Elliptic Curves: Number Theory and Cryptography*** —
  more crypto-oriented; good chapter on Weil pairing + embedding
  degree.
- **Sutherland's MIT course 18.783** — slides + notes free online.
  Lectures 14–20 cover isogeny volcanoes, CM, and class-group action
  in a hands-on computational style.  Highly recommended.
