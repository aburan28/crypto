import Mathlib.Tactic.Ring

/-!
An algebraic certificate that three distinct abscissae on a quadratic
interpolant with nonzero leading coefficient cannot lie on an affine line.
The P-192 coefficients and the genus-three Riemann--Roch bridge are checked
and explained separately.
-/

namespace P192ExtensionCover

variable {F : Type*} [Field F]

/-- A nonzero quadratic interpolation coefficient prevents three points with
    distinct x-coordinates from lying on one nonvertical affine line. -/
theorem quadratic_triple_noncollinear
    (a b c x₁ x₂ x₃ : F) (ha : a ≠ 0)
    (h₁₂ : x₁ ≠ x₂) (h₁₃ : x₁ ≠ x₃) (h₂₃ : x₂ ≠ x₃) :
    ¬ ∃ m n : F,
      a * x₁ ^ 2 + b * x₁ + c = m * x₁ + n ∧
      a * x₂ ^ 2 + b * x₂ + c = m * x₂ + n ∧
      a * x₃ ^ 2 + b * x₃ + c = m * x₃ + n := by
  rintro ⟨m, n, h₁, h₂, h₃⟩
  have h₁₂zero : (x₁ - x₂) * (a * (x₁ + x₂) + b - m) = 0 := by
    calc
      _ = (a * x₁ ^ 2 + b * x₁ + c - (m * x₁ + n)) -
          (a * x₂ ^ 2 + b * x₂ + c - (m * x₂ + n)) := by ring
      _ = 0 := by rw [sub_eq_zero.mpr h₁, sub_eq_zero.mpr h₂]; ring
  have h₁₃zero : (x₁ - x₃) * (a * (x₁ + x₃) + b - m) = 0 := by
    calc
      _ = (a * x₁ ^ 2 + b * x₁ + c - (m * x₁ + n)) -
          (a * x₃ ^ 2 + b * x₃ + c - (m * x₃ + n)) := by ring
      _ = 0 := by rw [sub_eq_zero.mpr h₁, sub_eq_zero.mpr h₃]; ring
  have h₁₂' : a * (x₁ + x₂) + b - m = 0 :=
    (mul_eq_zero.mp h₁₂zero).resolve_left (sub_ne_zero.mpr h₁₂)
  have h₁₃' : a * (x₁ + x₃) + b - m = 0 :=
    (mul_eq_zero.mp h₁₃zero).resolve_left (sub_ne_zero.mpr h₁₃)
  have h₂₃zero : a * (x₂ - x₃) = 0 := by
    calc
      _ = (a * (x₁ + x₂) + b - m) - (a * (x₁ + x₃) + b - m) := by ring
      _ = 0 := by rw [h₁₂', h₁₃']; ring
  exact (mul_eq_zero.mp h₂₃zero).elim ha (fun h => h₂₃ (sub_eq_zero.mp h))

end P192ExtensionCover

#print axioms P192ExtensionCover.quadratic_triple_noncollinear
