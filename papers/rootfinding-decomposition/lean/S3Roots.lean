import Mathlib.Tactic

/-!
# The roots of Semaev's `S₃(x_R, x_F, X)` are `x(R + F)` and `x(R − F)`

Companion to `papers/rootfinding-decomposition/rootfinding_decomposition.tex`,
Proposition 1.

Setting: a field `K` and the short-Weierstrass curve `y² = x³ + a x + b`.
For affine points `R = (xR, yR)` and `F = (xF, yF)` on the curve with
`xR ≠ xF`, the chord through `R` and `F` has slope
`λ₊ = (yF − yR)/(xF − xR)` and the chord through `R` and `−F = (xF, −yF)`
has slope `λ₋ = (−yF − yR)/(xF − xR)`; the affine addition law gives
`x(R + F) = λ₊² − xR − xF` and `x(R − F) = λ₋² − xR − xF`.
(These are the chord formulas of the affine group law; in Mathlib's
`WeierstrassCurve.Affine` they are `slope` and `addX` with `a₁ = a₂ = 0`.)

Semaev's third summation polynomial for this curve is
`S₃(x₁, x₂, x₃) = (x₁ − x₂)² x₃² − 2[(x₁ + x₂)(x₁x₂ + a) + 2b] x₃
                 + [(x₁x₂ − a)² − 4b(x₁ + x₂)]`.

Main theorem `S3_factor`: as a function of `X`,
`S₃(xR, xF, X) = (xR − xF)² · (X − x(R+F)) · (X − x(R−F))`.
Since `(xR − xF)² ≠ 0`, this says the roots of `X ↦ S₃(xR, xF, X)` are
exactly `x(R+F)` and `x(R−F)`, counted with multiplicity (`S3_root_iff`).
No characteristic assumption is needed for the identity itself.
-/

namespace SemaevS3

variable {K : Type*} [Field K]

/-- Semaev's third summation polynomial for `y² = x³ + a x + b`,
evaluated at `(x₁, x₂, x₃)`. -/
def S3 (a b x₁ x₂ x₃ : K) : K :=
  (x₁ - x₂) ^ 2 * x₃ ^ 2
    - 2 * ((x₁ + x₂) * (x₁ * x₂ + a) + 2 * b) * x₃
    + ((x₁ * x₂ - a) ^ 2 - 4 * b * (x₁ + x₂))

/-- Slope of the chord through `(xR, yR)` and `(xF, yF)` (the affine group-law slope). -/
def slope (xR yR xF yF : K) : K := (yF - yR) / (xF - xR)

/-- `x`-coordinate of the sum of two affine points with distinct `x`
(the affine chord formula `λ² − x₁ − x₂`). -/
def addX (xR yR xF yF : K) : K := slope xR yR xF yF ^ 2 - xR - xF

/-- `x(R + F)`. -/
abbrev xPlus (xR yR xF yF : K) : K := addX xR yR xF yF
/-- `x(R − F)`, i.e. the sum with `−F = (xF, −yF)`. -/
abbrev xMinus (xR yR xF yF : K) : K := addX xR yR xF (-yF)

/-- Clearing the denominator in the chord formula. -/
lemma d2_mul_addX (xR yR xF yF : K) (hx : xF - xR ≠ 0) :
    (xF - xR) ^ 2 * addX xR yR xF yF = (yF - yR) ^ 2 - (xF - xR) ^ 2 * (xR + xF) := by
  unfold addX slope
  rw [div_pow]
  field_simp
  ring

/-- **Proposition 1 (factorization form).** On the curve, with `xR ≠ xF`,
`S₃(xR, xF, X) = (xR − xF)² (X − x(R+F)) (X − x(R−F))` for every `X`. -/
theorem S3_factor (a b xR yR xF yF X : K) (hx : xR ≠ xF)
    (hR : yR ^ 2 = xR ^ 3 + a * xR + b) (hF : yF ^ 2 = xF ^ 3 + a * xF + b) :
    S3 a b xR xF X
      = (xR - xF) ^ 2 * (X - xPlus xR yR xF yF) * (X - xMinus xR yR xF yF) := by
  have hd : xF - xR ≠ 0 := sub_ne_zero.mpr (Ne.symm hx)
  have hd2 : (xF - xR) ^ 2 ≠ 0 := pow_ne_zero 2 hd
  -- Clear the two chord denominators.
  have hp := d2_mul_addX xR yR xF yF hd          -- d² · x(R+F)
  have hm := d2_mul_addX xR yR xF (-yF) hd       -- d² · x(R−F)
  -- Multiply the goal by d² ≠ 0 and reduce to a polynomial identity.
  apply mul_left_cancel₀ hd2
  have expand :
      (xF - xR) ^ 2 * ((xR - xF) ^ 2 * (X - xPlus xR yR xF yF) * (X - xMinus xR yR xF yF))
        = ((xF - xR) ^ 2 * X - (xF - xR) ^ 2 * addX xR yR xF yF)
          * ((xF - xR) ^ 2 * X - (xF - xR) ^ 2 * addX xR yR xF (-yF)) := by
    unfold xPlus xMinus; ring
  rw [expand, hp, hm]
  unfold S3
  -- The remaining identity holds modulo the two curve equations; the
  -- cofactors were computed symbolically (sympy) and are checked by `ring`.
  linear_combination
    (2 * X * xF ^ 2 - 4 * X * xF * xR + 2 * X * xR ^ 2 - a * xR - b + 2 * xF ^ 3
        - 2 * xF ^ 2 * xR - 2 * xF * xR ^ 2 + xR ^ 3 + 2 * yF ^ 2 - yR ^ 2) * hR
    + (2 * X * xF ^ 2 - 4 * X * xF * xR + 2 * X * xR ^ 2 - a * xF + 2 * a * xR + b
        + xF ^ 3 - 2 * xF ^ 2 * xR - 2 * xF * xR ^ 2 + 4 * xR ^ 3 - yF ^ 2) * hF

/-- **Proposition 1 (root form).** On the curve, with `xR ≠ xF`,
`S₃(xR, xF, X) = 0` iff `X = x(R+F)` or `X = x(R−F)`. -/
theorem S3_root_iff (a b xR yR xF yF X : K) (hx : xR ≠ xF)
    (hR : yR ^ 2 = xR ^ 3 + a * xR + b) (hF : yF ^ 2 = xF ^ 3 + a * xF + b) :
    S3 a b xR xF X = 0 ↔ X = xPlus xR yR xF yF ∨ X = xMinus xR yR xF yF := by
  rw [S3_factor a b xR yR xF yF X hx hR hF]
  have hd2 : (xR - xF) ^ 2 ≠ 0 := pow_ne_zero 2 (sub_ne_zero.mpr hx)
  constructor
  · intro h
    rcases mul_eq_zero.mp h with h1 | h2
    · rcases mul_eq_zero.mp h1 with h0 | h1'
      · exact absurd h0 hd2
      · exact Or.inl (sub_eq_zero.mp h1')
    · exact Or.inr (sub_eq_zero.mp h2)
  · rintro (h | h)
    · rw [h]; ring
    · rw [h]; ring

end SemaevS3
