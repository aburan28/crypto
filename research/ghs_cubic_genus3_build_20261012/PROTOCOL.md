# Protocol: an explicit genus-3 curve over `F_q` for an elliptic curve over `F_{q³}`

Frozen 2026-10-12, before `genus3_build.rs` was run. Conductor T-2267.
**Structural stage diagnostic on toy fields.** It builds a curve and checks
its Frobenius. It solves no discrete logarithm, runs no index calculus and
claims no speedup. Follows `research/ghs_cubic_cover_census_20261012/`
(PR #1759), which counted these covers without building them.

## Construction under test (the `m = 2, r = 6` type)

`K = F_{q³}`, `σ` the `q`-Frobenius. Take a monic quadratic `u ∈ K[t]` whose
six conjugate roots (those of `u`, `u' = u^σ`, `u'' = u^{σ²}`) are distinct,
with `u, u', u''` linearly independent. Put

    E0 : y0² = u·u'',   E1 : y1² = u·u',   E2 : y2² = u'·u''

so `E1 = E0^σ`, `E2 = E0^{σ²}`. Let `C` be the curve with function field
`K(t, y0, y1)`, and `y2 = y0·y1/u`.

- **Descent datum (derived).** `σ~ : t ↦ t, y0 ↦ y1 ↦ y2 ↦ y0`, acting as
  `σ` on `K`, has order 3. So `C` has a model over `F_q`.
- **Plane model (derived).** With `X = 1/y0, Y = 1/y1, Z = 1/y2` one has
  `(u : u' : u'') = (Z² : X² : Y²)`. If `(1, t, t²)ᵀ = R·(u, u', u'')ᵀ`,
  then `C` is the plane quartic
  `G = (R₁·w)² − (R₀·w)(R₂·w)`, `w = (Z², X², Y²)`.
  `σ~` permutes `X → Y → Z`, so in the coordinates
  `W_j = ω_j X + ω_j^σ Y + ω_j^{σ²} Z` (`ω_j = 1, g, g²`) the quartic has
  coefficients in `F_q` up to one scalar.
- **Transfer (derived).** `φ = (1 + Fr + Fr²)∘π0^* : E0(K) → Jac(C)(F_q)`,
  with `Fr(t, y0, y1, y2) = (t^q, y2^q, y0^q, y1^q)`. Since
  `π0_* ∘ φ = [2]`, the kernel of `φ` lies in `E0[2]`.

## Pre-registered predictions

- **P1 (descent).** Every coefficient of the quartic in `W` coordinates,
  after dividing by the first non-zero one, lies in `F_q`.
- **P2 (smooth).** The quartic is smooth (Ciani conditions on the `K`
  model), hence `C` is not hyperelliptic.
- **P3 (Frobenius).** Counting points of the `F_q` quartic over `F_q`,
  `F_{q²}`, `F_{q³}` gives `N1 = q+1`, `N2 = q²+1`,
  `N3 = q³+1−3t_E` with `t_E = q³+1−#E0(K)`, i.e. L-polynomial
  `1 − t_E T³ + q³T⁶` and `#Jac(C)(F_q) = #E0(K)`.
- **P4 (transfer divisors).** For every `P = (t0, y0) ∈ E0(K)` with
  `u·u'(t0)` a non-zero square, the six points of `π0^*(P)` and its two
  Frobenius images lie on the `F_q` quartic, and the coordinate-wise
  `q`-power of each equals the model's `Fr` image (the divisor is
  `F_q`-rational).

## Falsifiers

P1: a coefficient outside `F_q`. P2: a failed Ciani condition. P3: any
`N_k` off, or `#Jac(C)(F_q) ≠ #E0(K)`. P4: any point off the curve or any
Frobenius mismatch. A falsified P1 or P3 means the census of PR #1759
counted covers that do not descend as assumed.

## Frozen inputs

`u = t² + b t + c`, `(b, c)` in lexicographic order of field index. Per
field the first cases with `#E0(K) = 4ℓ`, `ℓ` prime:

| field | split `u` (roots in `K`) | irreducible `u` |
| --- | --- | --- |
| q = 7 | first 3 | first 2 |
| q = 5 | first 1 | first 1 |
| q = 11 | first 1 | first 1 |

## Not measured

Linear equivalence of divisors on the quartic (no class arithmetic on a
plane quartic is implemented), hence no numerical test that `φ` is a
homomorphism or injective on `E0[ℓ]`: both rest on the derivation. No
relation collection, no DLP, no cost. The `m = 3, r = 5` hyperelliptic
type is not built.

## Language

Rust (`genus3_build.rs`, dependency-free, `rustc -O`). Python only plots.
