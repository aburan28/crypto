---

## 5. Worked example (toy field, every number from `fiber_probe.json` / `pullback_factor.out`)

**Field and curve.** p = 211. E₁: y² = x³ + x + 3.
#E₁(F_p) = N = 228 = 2²·3·19, trace t = −16, target prime r = 19, cofactor
h = 12. E₁(F_p) has exactly one rational point of order 2, T₂ = (168, 0), and a
rational point of order 3, T₃ = (77, 208).

**2-isogeny.** φ₂: E₁ → E₂ = E₁/⟨T₂⟩: y² = x³ + 113x + 97 (Vélu).
Measured: #E₂(F_p) = 228; image φ₂(E₁(F_p)) has 114 points (index 2 = |M^π|);
rational-fiber size histogram over E₂(F_p): {2: 114 points, 0: 114 points}.
Homomorphism and on-curve checks: 50/50.

**3-isogeny with rational kernel.** φ₃: E₁ → E₃ = E₁/⟨T₃⟩:
y² = x³ + 205x + 178. Image has 76 points (index 3); histogram
{3: 76, 0: 152}. For a point Q outside the image, the cubic numerator of
X(x) − x(Q) has no root in F_p; adjoining a root z gives F_{p³}, and
f(z) = z³ + z + 3 is a square there, so P = (z, √f(z)) ∈ E₁(F_{p³}) with
φ₃(P) = Q. The Frobenius orbit of P has length 3 and equals {P, P+T₃, P−T₃};
the cocycle π(P) − P is T₃ exactly. This is Prop. 2.4's "single orbit of
length ℓ over F_{q^ℓ}" case, observed.

**Pullback factorisation (Thm. 4.2).**
* ℓ = 2: numerator of S₃^{E₂}(X(x₁), X(x₂), X(x₃)) has degree 4 in each
  variable, 112 terms; divisible by S₃^{E₁} with zero remainder; cofactor of
  degree (2,2,2), 27 terms, symmetric, vanishing on 50/50 forced
  non-degenerate triples with P₁+P₂+P₃ = T₂.
* ℓ = 3: degree 6 per variable, 330 terms; divisible by S₃^{E₁}; cofactor of
  degree (4,4,4), 113 terms, symmetric, vanishing on 50/50 forced triples.

**Verified lifted relation.** Factor base F₁ = {P ∈ E₁(F_p) : x(P) < 12},
|F₁| = 14. Random target R = (4, 156), φ₃(R) = (27, 169). A 3-term relation
on E₃ over φ₃(F₁): φ₃(R) = 3·φ₃((1, 65)). On E₁: 3·(1, 65) − R = (77, 208) = T₃,
i.e. the relation holds only modulo the kernel (kernel index 1 in
{O, T₃, −T₃}). After cofactor projection: [12]·3·(1, 65) = [12]·R, equality
verified on E₁. So log_{P₀}(R) ≡ 3·log_{P₀}((1, 65)) (mod 19) is a valid
relation in G obtained from a decomposition on the *isogenous* curve, with
the kernel ambiguity resolved at zero cost (Prop. 4.1).

**What the example shows and does not show.** It confirms the torsor/affine
Frobenius structure, the index formula, the pullback factorisation, and
zero-cost lifting. It does not — and at this size cannot — show any cost
difference between curves; Prop. 3.2 says there is none in the system shape.
