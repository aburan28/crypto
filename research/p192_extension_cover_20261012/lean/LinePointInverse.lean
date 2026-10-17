import Mathlib.Tactic.LinearCombination
import Mathlib.Tactic.Ring

/-!
Algebraic core of the generic plane-point inverse used by the P-192 Type I
line-support screen. Concrete roots, the plane quartic, polynomial gcd, and
elliptic quotient maps are checked by Sage and PARI separately.
-/

namespace P192ExtensionCover

variable {F : Type*} [Field F]

/-- The three conjugate quartic right-hand sides use each of the six branch
    factors twice, so their product is the square of the norm product. -/
theorem three_quartics_product_is_norm_square
    (X z₀ z₁ z₂ b₀ b₁ b₂ : F) :
    ((X - z₀) * (X - z₂) * (X - b₀) * (X - b₂)) *
    ((X - z₁) * (X - z₀) * (X - b₁) * (X - b₀)) *
    ((X - z₂) * (X - z₁) * (X - b₂) * (X - b₁)) =
    ((X - z₀) * (X - z₁) * (X - z₂) *
      (X - b₀) * (X - b₁) * (X - b₂)) ^ 2 := by
  ring

/-- The squared first-quotient equation and the three-factor norm equation
    determine the common projective scale used to recover the cover point. -/
theorem scale_fixed_by_square_and_norm
    (lam A₀ A₁ A₂ f₀ N : F)
    (hsquare : lam ^ 2 * A₀ ^ 2 * f₀ = 1)
    (hnorm : lam ^ 3 * A₀ * A₁ * A₂ * N = 1)
    (hden : A₁ * A₂ * N ≠ 0) :
    lam = A₀ * f₀ / (A₁ * A₂ * N) := by
  apply (eq_div_iff hden).2
  linear_combination (A₀ * f₀) * hnorm - (lam * A₁ * A₂ * N) * hsquare

end P192ExtensionCover

#print axioms P192ExtensionCover.three_quartics_product_is_norm_square
#print axioms P192ExtensionCover.scale_fixed_by_square_and_norm
