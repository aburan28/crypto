import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Ring

/-!
The algebraic core of the P-192 Legendre obstruction.  If Frobenius cycles
three distinct 2-division roots, their three cross-ratios have product -1.
The finite-field Frobenius and elliptic-curve instantiation is recorded as an
open formalization bridge in the adjoining report.
-/

namespace P192ExtensionCover

theorem cyclic_cross_ratio_product {K : Type*} [Field K]
    (a b c : K) (hab : a ≠ b) (hbc : b ≠ c) (hca : c ≠ a) :
    ((c - a) / (b - a)) * ((b - c) / (a - c)) * ((a - b) / (c - b)) = -1 := by
  have hba0 : b - a ≠ 0 := sub_ne_zero.mpr (Ne.symm hab)
  have hac0 : a - c ≠ 0 := sub_ne_zero.mpr (Ne.symm hca)
  have hcb0 : c - b ≠ 0 := sub_ne_zero.mpr (Ne.symm hbc)
  field_simp [hba0, hac0, hcb0]
  ring

theorem cyclic_cross_ratio_not_norm_one {K : Type*} [Field K]
    (a b c : K) (hab : a ≠ b) (hbc : b ≠ c) (hca : c ≠ a)
    (hchar : (-1 : K) ≠ 1) :
    ((c - a) / (b - a)) * ((b - c) / (a - c)) * ((a - b) / (c - b)) ≠ 1 := by
  rw [cyclic_cross_ratio_product a b c hab hbc hca]
  exact hchar

/-- The group-theoretic step for a degree-two cover: on any subgroup where
    doubling is injective, a pullback with pushforward-pullback equal to
    doubling is injective. -/
theorem pullback_injective_of_push_pull_two {A B : Type*} [AddCommGroup A]
    (pull : A → B) (push : B → A)
    (hpushpull : ∀ x, push (pull x) = x + x)
    (hdouble : Function.Injective (fun x : A => x + x)) :
    Function.Injective pull := by
  intro x y hxy
  apply hdouble
  calc
    x + x = push (pull x) := (hpushpull x).symm
    _ = push (pull y) := congrArg push hxy
    _ = y + y := hpushpull y

end P192ExtensionCover

#print axioms P192ExtensionCover.cyclic_cross_ratio_product
#print axioms P192ExtensionCover.cyclic_cross_ratio_not_norm_one
#print axioms P192ExtensionCover.pullback_injective_of_push_pull_two
