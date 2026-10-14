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

/-- Clearing denominators in equality of two Legendre j-invariants yields
    the six cross-ratio transformations. This is the algebraic part of the
    nonsplit-quadratic branch's exact-degree screen. -/
theorem legendre_j_difference_factor {K : Type*} [CommRing K] (x y : K) :
    (x ^ 2 - x + 1) ^ 3 * y ^ 2 * (y - 1) ^ 2 -
      (y ^ 2 - y + 1) ^ 3 * x ^ 2 * (x - 1) ^ 2 =
      (x - y) * (x + y - 1) * (x * y - x - y) *
        (x * y - x + 1) * (x * y - y + 1) * (x * y - 1) := by
  ring

/-- The polynomial identity that converts an elliptic cubic with a root to
    a palindromic Scholten cubic. Field and nonzero hypotheses needed to
    construct the scale are checked in the Sage instance. -/
theorem scholten_palindromic_identity {K : Type*} [CommRing K]
    (a b z s t : K) (hroot : z ^ 3 + a * z + b = 0)
    (hscale : 3 * z * s + 3 * z ^ 2 + a = 0) :
    (z + s * (t + 1)) ^ 3 + a * (z + s * (t + 1)) + b =
      s ^ 3 * t ^ 3 + (3 * s ^ 3 + 3 * z * s ^ 2) * t ^ 2 +
        (3 * s ^ 3 + 3 * z * s ^ 2) * t + s ^ 3 := by
  calc
    (z + s * (t + 1)) ^ 3 + a * (z + s * (t + 1)) + b =
        s ^ 3 * t ^ 3 + (3 * s ^ 3 + 3 * z * s ^ 2) * t ^ 2 +
          (3 * s ^ 3 + 3 * z * s ^ 2) * t + s ^ 3 +
          (z ^ 3 + a * z + b) + s * (t + 1) * (3 * z * s + 3 * z ^ 2 + a) := by
            ring
    _ = _ := by rw [hroot, hscale]; ring

/-- Algebraic verification of the degree-two projection from the plane
    quartic to the original elliptic cubic. -/
theorem quartic_projection_identity {K : Type*} [CommRing K]
    (a b x z : K) (h : z ^ 4 = x ^ 3 + a * x + b) :
    (z ^ 2) ^ 2 = x ^ 3 + a * x + b := by
  calc
    (z ^ 2) ^ 2 = z ^ 4 := by ring
    _ = x ^ 3 + a * x + b := h

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
#print axioms P192ExtensionCover.legendre_j_difference_factor
#print axioms P192ExtensionCover.scholten_palindromic_identity
#print axioms P192ExtensionCover.quartic_projection_identity
