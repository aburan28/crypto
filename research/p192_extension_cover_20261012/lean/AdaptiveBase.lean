import ScalarAnchor

/-!
An adaptive collector can force a relation by defining one new base point
from each known target and old-base summands. Every such relation introduces
one private logarithm column. The theorem below shows that these equations
extend every possible assignment of old logarithms uniquely, so they do not
constrain the old columns until a second equation eliminates a private one.
-/

namespace P192ExtensionCover

variable {K : Type*} [AddCommGroup K]

/-- A family of known-scalar equations, each with one fresh private column,
    has a unique extension for every assignment of the old columns. The
    functions `row` can be arbitrary, including nonlinear expressions; the
    conclusion depends only on each new column occurring once with unit
    coefficient and nowhere in the other rows. -/
theorem private_columns_extend {n m : ℕ}
    (row : Fin m → (Fin n → K) → K) (rhs : Fin m → K) :
    ∀ old : Fin n → K,
      ∃! fresh : Fin m → K, ∀ j, row j old + fresh j = rhs j := by
  intro old
  refine ⟨fun j => rhs j - row j old, ?_, ?_⟩
  · intro j
    dsimp
    abel
  · intro fresh h
    funext j
    calc
      fresh j = (row j old + fresh j) - row j old := by abel
      _ = rhs j - row j old := by rw [h j]

/-- Existence of old solutions is exactly existence of augmented solutions
    when the added equations have private columns. -/
theorem private_columns_preserve_old_feasibility {n m : ℕ}
    (oldCondition : (Fin n → K) → Prop)
    (row : Fin m → (Fin n → K) → K) (rhs : Fin m → K) :
    (∃ old : Fin n → K, oldCondition old) ↔
      ∃ old : Fin n → K, ∃ fresh : Fin m → K, oldCondition old ∧
        ∀ j, row j old + fresh j = rhs j := by
  constructor
  · rintro ⟨old, hold⟩
    obtain ⟨fresh, hfresh, _⟩ := private_columns_extend row rhs old
    exact ⟨old, fresh, hold, hfresh⟩
  · rintro ⟨old, _, hold, _⟩
    exact ⟨old, hold⟩

end P192ExtensionCover

#print axioms P192ExtensionCover.private_columns_extend
#print axioms P192ExtensionCover.private_columns_preserve_old_feasibility
