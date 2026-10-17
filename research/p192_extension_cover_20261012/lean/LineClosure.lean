import Mathlib.LinearAlgebra.Span.Defs

/-!
When every new geometric relation row repeats an existing incidence row,
the row span cannot grow. The P-192 second-generation census provides the
concrete row-identification input for this abstract result.
-/

namespace P192ExtensionCover

theorem duplicate_relation_rows_span_le
    {R M I J : Type*} [Semiring R] [AddCommMonoid M] [Module R M]
    (old : I → M) (new : J → M)
    (hrepeat : ∀ j, ∃ i, new j = old i) :
    Submodule.span R (Set.range new) ≤ Submodule.span R (Set.range old) := by
  apply Submodule.span_mono
  rintro x ⟨j, rfl⟩
  obtain ⟨i, hi⟩ := hrepeat j
  exact ⟨i, hi.symm⟩

end P192ExtensionCover

#print axioms P192ExtensionCover.duplicate_relation_rows_span_le
