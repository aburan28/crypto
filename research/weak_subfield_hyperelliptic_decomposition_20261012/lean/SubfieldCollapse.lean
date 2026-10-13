/-!
# Subfield collapse of the weak-curve factor base

Bhati–Singh choose `μ` so that `ψ = μ² · f̃` has coefficients in the subfield
`k = F_q` of `K = F_{q^n}`, and collect relations from `λ² − ψ` with `λ` over
`k`. This file proves, with no Mathlib, the three steps of the argument that
their factor base then spans only what the subfield curve `y'² = ψ(x)` sees.

1. `point_map`: `y ↦ μ·y` carries points of `y² = f̃(x)` to points of
   `y'² = ψ(x)`.
2. `relation_point_rational`: a relation point (`x ∈ k`, `μ·y = λ(x)`,
   `λ` preserving `k`) has both coordinates in `k` in the `ψ` model.
3. `span_le_image`: in any additive structure, if every generator is in the
   image of a map `φ` that preserves `0`, `+` and `−`, so is everything the
   generators span. With `φ : J_{C0}(k) → J_C(K)` this is the collapse.

Only multiplication's associativity and commutativity are used for (1) and
(2), so they hold in every commutative ring, in particular in `K`.

Not formalized (the bridge to the actual objects): that `φ` exists as the
base-change map of Jacobians and sends the class of a `k`-rational point of
`C0` to the class `(P) − (∞)` on `C`; that `#J_{C0}(F_q) ≤ (1 + √q)^{2g}`.
-/

namespace SubfieldCollapse

/-- The multiplicative laws used. Every commutative ring is an instance. -/
structure CommMul (K : Type) where
  mul : K → K → K
  mul_assoc : ∀ a b c, mul (mul a b) c = mul a (mul b c)
  mul_comm : ∀ a b, mul a b = mul b a

variable {K : Type} (S : CommMul K)

theorem mul_mul_mul_comm (a b c d : K) :
    S.mul (S.mul a b) (S.mul c d) = S.mul (S.mul a c) (S.mul b d) := by
  rw [S.mul_assoc a b (S.mul c d), ← S.mul_assoc b c d, S.mul_comm b c,
    S.mul_assoc c b d, ← S.mul_assoc a c (S.mul b d)]

/-- Step 1. If `μ² · f̃(x) = ψ(x)` for all `x`, then `(x, y)` on `y² = f̃(x)`
gives `(x, μ·y)` on `y'² = ψ(x)`. -/
theorem point_map (μ : K) (ft ψ : K → K)
    (hψ : ∀ x, S.mul (S.mul μ μ) (ft x) = ψ x)
    (x y : K) (hy : S.mul y y = ft x) :
    S.mul (S.mul μ y) (S.mul μ y) = ψ x := by
  rw [mul_mul_mul_comm S μ y μ y, hy, hψ x]

/-- Step 2. A relation point of the paper (`x` in the subfield, cut out by
`μ·y = λ(x)` with `λ` over the subfield) is rational over the subfield in
the `ψ` model, and lies on it. -/
theorem relation_point_rational (k : K → Prop) (μ : K) (ft ψ lam : K → K)
    (hψ : ∀ x, S.mul (S.mul μ μ) (ft x) = ψ x)
    (hlam : ∀ x, k x → k (lam x))
    (x y : K) (hx : k x) (hy : S.mul y y = ft x) (hrel : S.mul μ y = lam x) :
    k x ∧ k (S.mul μ y) ∧ S.mul (S.mul μ y) (S.mul μ y) = ψ x :=
  ⟨hx, by rw [hrel]; exact hlam x hx, point_map S μ ft ψ hψ x y hy⟩

/-- An additive structure: just the operations, no laws needed below. -/
structure AddOps (G : Type) where
  zero : G
  add : G → G → G
  neg : G → G

/-- The set spanned by generators under `0`, `+`, `−`. -/
inductive Span {G : Type} (A : AddOps G) (gens : G → Prop) : G → Prop
  | gen (g : G) : gens g → Span A gens g
  | zero : Span A gens A.zero
  | add (a b : G) : Span A gens a → Span A gens b → Span A gens (A.add a b)
  | neg (a : G) : Span A gens a → Span A gens (A.neg a)

/-- Step 3. If every generator of a subgroup of `G` is in the image of a map
`φ : G0 → G` preserving `0`, `+`, `−`, the whole span is in that image. -/
theorem span_le_image {G0 G : Type} (A0 : AddOps G0) (A : AddOps G) (φ : G0 → G)
    (h0 : φ A0.zero = A.zero)
    (hadd : ∀ a b, φ (A0.add a b) = A.add (φ a) (φ b))
    (hneg : ∀ a, φ (A0.neg a) = A.neg (φ a))
    (gens : G → Prop) (hgens : ∀ g, gens g → ∃ g0, φ g0 = g) :
    ∀ g, Span A gens g → ∃ g0, φ g0 = g := by
  intro g hg
  induction hg with
  | gen g hg => exact hgens g hg
  | zero => exact ⟨A0.zero, h0⟩
  | add a b _ _ iha ihb =>
    obtain ⟨a0, ha⟩ := iha
    obtain ⟨b0, hb⟩ := ihb
    exact ⟨A0.add a0 b0, by rw [hadd, ha, hb]⟩
  | neg a _ iha =>
    obtain ⟨a0, ha⟩ := iha
    exact ⟨A0.neg a0, by rw [hneg, ha]⟩

/-- Consequence for a discrete logarithm: an element outside the image of
`φ` is not in the span of the factor base, so no relation among factor-base
elements expresses it. -/
theorem not_in_span_of_not_in_image {G0 G : Type} (A0 : AddOps G0) (A : AddOps G)
    (φ : G0 → G) (h0 : φ A0.zero = A.zero)
    (hadd : ∀ a b, φ (A0.add a b) = A.add (φ a) (φ b))
    (hneg : ∀ a, φ (A0.neg a) = A.neg (φ a))
    (gens : G → Prop) (hgens : ∀ g, gens g → ∃ g0, φ g0 = g)
    (D : G) (hD : ¬ ∃ g0, φ g0 = D) : ¬ Span A gens D :=
  fun h => hD (span_le_image A0 A φ h0 hadd hneg gens hgens D h)

end SubfieldCollapse

#print axioms SubfieldCollapse.point_map
#print axioms SubfieldCollapse.relation_point_rational
#print axioms SubfieldCollapse.span_le_image
#print axioms SubfieldCollapse.not_in_span_of_not_in_image
