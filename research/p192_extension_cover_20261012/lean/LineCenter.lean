import Mathlib.Tactic.Abel

/-!
Four-point line equations with one common quotient-group constant can be
centered by its fourth part. Translation preserves distinct point labels;
it does not change the incidence pattern or remove private columns.
-/

namespace P192ExtensionCover

theorem four_point_line_centered_zero {G : Type*} [AddCommGroup G]
    (a b c d center : G)
    (hsum : a + b + c + d = center + center + center + center) :
    (a - center) + (b - center) + (c - center) + (d - center) = 0 := by
  calc
    (a - center) + (b - center) + (c - center) + (d - center)
        = (a + b + c + d) - (center + center + center + center) := by abel
    _ = 0 := by rw [hsum]; simp

theorem translate_point_labels_injective {G : Type*} [AddGroup G]
    (center : G) : Function.Injective (fun q : G => q - center) := by
  intro a b h
  have h' := congrArg (fun q : G => q + center) h
  simpa using h'

end P192ExtensionCover

#print axioms P192ExtensionCover.four_point_line_centered_zero
#print axioms P192ExtensionCover.translate_point_labels_injective
