import Mathlib.Algebra.BigOperators.Group.Finset.Basic
import Mathlib.Data.Fin.Basic
import Mathlib.Data.Fintype.Basic
import Mathlib.Tactic.NormNum

/-!
Each split-line equation has a private fresh-point column. If a linear
combination of line equations cancels all fresh columns, every coefficient
vanishes. This is the structural rank rule used in the P-192 line census;
the concrete incidence and geometric line identities are separate inputs.
-/

namespace P192ExtensionCover

theorem private_column_coefficients_zero {R : Type*} [Semiring R]
    {m : ℕ} (row : Fin m → Fin m → R) (coeff : Fin m → R)
    (hdiag : ∀ j, row j j = 1)
    (hoffdiag : ∀ i j, i ≠ j → row i j = 0)
    (hcancel : ∀ j, (∑ i : Fin m, coeff i * row i j) = 0) :
    ∀ j, coeff j = 0 := by
  classical
  intro j
  have hsum : (∑ i : Fin m, coeff i * row i j) = coeff j := by
    rw [Finset.sum_eq_single j]
    · simp [hdiag j]
    · intro i _ hij
      simp [hoffdiag i j hij]
    · simp
  have hj := hcancel j
  rw [hsum] at hj
  exact hj

theorem full_split_line_point_count : 141 + 2 * 4892 = 9925 := by
  norm_num

end P192ExtensionCover

#print axioms P192ExtensionCover.private_column_coefficients_zero
#print axioms P192ExtensionCover.full_split_line_point_count
