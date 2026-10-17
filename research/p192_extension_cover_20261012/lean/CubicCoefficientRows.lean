import Mathlib.LinearAlgebra.LinearIndependent.Defs

/-!
Expanding a field-valued cubic evaluation into prime-field coordinate rows
preserves its vanishing when the coordinate map is injective. The concrete
P-192 field basis and cubic evaluations are replayed by Sage.
-/

namespace P192ExtensionCover

variable {K L : Type*} [Field K] [AddCommGroup L] [Module K L]

theorem cubic_evaluation_zero_iff_coordinate_rows_zero {n : ℕ}
    (coordinates : L →ₗ[K] (Fin n → K))
    (hinjective : Function.Injective coordinates) (evaluation : L) :
    evaluation = 0 ↔ ∀ i, coordinates evaluation i = 0 := by
  constructor
  · intro h i
    simp [h]
  · intro h
    apply hinjective
    funext i
    simpa using h i

end P192ExtensionCover

#print axioms P192ExtensionCover.cubic_evaluation_zero_iff_coordinate_rows_zero
