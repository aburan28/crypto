import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Ring
import Mathlib.Data.Finset.Card
import Mathlib.Data.Finset.Fin
import Mathlib.Data.Finset.Prod
import Mathlib.Data.Fintype.Basic

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

/-- The coefficient-swapped Scholten cubic is the reciprocal of the
    original cubic. The field involution and genus-two descent are separate
    geometric obligations in the report. -/
theorem scholten_reciprocal_identity {K : Type*} [Field K]
    (A B c t : K) (ht : t ≠ 0) :
    A * t ^ 3 + B * c * t ^ 2 + B * c ^ 2 * t + A =
      t ^ 3 * (A * (t⁻¹) ^ 3 + B * c ^ 2 * (t⁻¹) ^ 2 + B * c * t⁻¹ + A) := by
  field_simp [ht]
  ring

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

/-- Projecting an exact three-point Jacobian relation by a cofactor also
    projects the known multiple of the generator. This is the group step
    used by the small-field relation-control implementation. -/
theorem projected_three_point_relation {A : Type*} [AddCommGroup A]
    (c k : ℤ) (T D P Q R : A) (h : T + k • D = P + Q + R) :
    c • T + (c * k) • D = c • P + c • Q + c • R := by
  calc
    c • T + (c * k) • D = c • (T + k • D) := by simp [zsmul_add, mul_zsmul]
    _ = c • (P + Q + R) := congrArg (fun X : A => c • X) h
    _ = c • P + c • Q + c • R := by simp [zsmul_add]

/-- Scalar recovery after cofactor projection must divide the sum of
    projected base logs by the cofactor, then remove the target offset. -/
theorem projected_target_scalar {K : Type*} [Field K]
    (c x k l₁ l₂ l₃ : K) (hc : c ≠ 0)
    (h : c * (x + k) = l₁ + l₂ + l₃) :
    x = (l₁ + l₂ + l₃) / c - k := by
  calc
    x = (c * (x + k)) / c - k := by field_simp [hc]; ring
    _ = (l₁ + l₂ + l₃) / c - k := by rw [h]

/-- An odd Frobenius power paired with its conjugate is a square. The
    finite-field bridge supplies the Frobenius-distributivity hypothesis
    in `type_i_subfield_rhs_square`. -/
theorem odd_power_pair_square {K : Type*} [CommRing K] (u : K) (n : ℕ) :
    u * u ^ (2 * n + 1) = (u ^ (n + 1)) ^ 2 := by
  calc
    u * u ^ (2 * n + 1) = u ^ ((2 * n + 1) + 1) := by rw [pow_succ]; ring
    _ = u ^ ((n + 1) * 2) := by congr 1; omega
    _ = (u ^ (n + 1)) ^ 2 := by rw [pow_mul]

/-- For an odd Frobenius power q, every subfield x-coordinate on the
    Type I quartic has a square right-hand side once Frobenius fixes x. -/
theorem type_i_subfield_rhs_square {K : Type*} [CommRing K]
    (x a b : K) (n : ℕ)
    (hconj : (x - a ^ (2 * n + 1)) * (x - b ^ (2 * n + 1)) =
      ((x - a) * (x - b)) ^ (2 * n + 1)) :
    (x - a) * (x - a ^ (2 * n + 1)) *
      (x - b) * (x - b ^ (2 * n + 1)) =
      (((x - a) * (x - b)) ^ (n + 1)) ^ 2 := by
  calc
    (x - a) * (x - a ^ (2 * n + 1)) *
        (x - b) * (x - b ^ (2 * n + 1)) =
        ((x - a) * (x - b)) *
          ((x - a ^ (2 * n + 1)) * (x - b ^ (2 * n + 1))) := by ring
    _ = ((x - a) * (x - b)) *
          (((x - a) * (x - b)) ^ (2 * n + 1)) := by rw [hconj]
    _ = (((x - a) * (x - b)) ^ (n + 1)) ^ 2 :=
      odd_power_pair_square ((x - a) * (x - b)) n

/-- A Frobenius-conjugate pair gives a fixed norm class whenever the
    induced action is an involution on the ambient divisor-class group. -/
theorem involution_norm_fixed {A : Type*} [AddCommGroup A]
    (σ : A →+ A) (hσ : ∀ x, σ (σ x) = x) (D : A) :
    σ (D + σ D) = D + σ D := by
  rw [map_add, hσ]
  exact add_comm (σ D) D

/-- The norm of an `r`-torsion divisor class remains `r`-torsion.  This is
    the group bridge used by the explicit Scholten norm computation. -/
theorem involution_norm_torsion {A : Type*} [AddCommGroup A]
    (σ : A →+ A) (r : ℤ) (D : A) (hD : r • D = 0) :
    r • (D + σ D) = 0 := by
  have hσD : r • σ D = 0 := by
    simpa only [map_zsmul, map_zero] using congrArg σ hD
  simp [zsmul_add, hD, hσD]

/-- Conorm of a torsion class remains torsion under an additive divisor
    class map; geometry must provide the particular conorm homomorphism. -/
theorem conorm_preserves_torsion {A B : Type*} [AddCommGroup A] [AddCommGroup B]
    (pull : A →+ B) (r : ℤ) (D : A) (hD : r • D = 0) :
    r • pull D = 0 := by
  simpa only [map_zsmul, map_zero] using congrArg pull hD

/-- Clearing the denominator in the quartic-to-cubic Möbius coordinate
    change. The exact P-192 coefficients and square twist are checked by
    the generator-map Sage receipt. -/
theorem mobius_quartic_clear {K : Type*} [Field K]
    (a b C D ρ X : K) (hX : X ≠ ρ) :
    ((C + D / (X - ρ)) ^ 3 + a * (C + D / (X - ρ)) + b) *
        (X - ρ) ^ 4 =
      (C * (X - ρ) + D) ^ 3 * (X - ρ) +
        a * (C * (X - ρ) + D) * (X - ρ) ^ 3 +
        b * (X - ρ) ^ 4 := by
  have h : X - ρ ≠ 0 := sub_ne_zero.mpr hX
  field_simp [h]

/-- Canonical unordered triples, with repetitions, from a 16-point base. -/
def canonicalTriples16 : Finset (Fin 16 × (Fin 16 × Fin 16)) :=
  ((Finset.univ : Finset (Fin 16)).product
    ((Finset.univ : Finset (Fin 16)).product
      (Finset.univ : Finset (Fin 16)))).filter
    (fun t => t.1 ≤ t.2.1 ∧ t.2.1 ≤ t.2.2)

set_option maxRecDepth 4096 in
theorem canonicalTriples16_card : canonicalTriples16.card = 816 := by
  decide

/-- A fixed 16-element base has at most 816 distinct three-sum images. -/
theorem three_sum_support_bound {A : Type*} [AddCommMonoid A] [DecidableEq A]
    (base : Fin 16 → A) :
    (canonicalTriples16.image
      (fun t => base t.1 + base t.2.1 + base t.2.2)).card ≤ 816 := by
  calc
    (canonicalTriples16.image
        (fun t => base t.1 + base t.2.1 + base t.2.2)).card ≤
        canonicalTriples16.card := Finset.card_image_le
    _ = 816 := canonicalTriples16_card

end P192ExtensionCover

#print axioms P192ExtensionCover.cyclic_cross_ratio_product
#print axioms P192ExtensionCover.cyclic_cross_ratio_not_norm_one
#print axioms P192ExtensionCover.pullback_injective_of_push_pull_two
#print axioms P192ExtensionCover.legendre_j_difference_factor
#print axioms P192ExtensionCover.scholten_palindromic_identity
#print axioms P192ExtensionCover.scholten_reciprocal_identity
#print axioms P192ExtensionCover.quartic_projection_identity
#print axioms P192ExtensionCover.projected_three_point_relation
#print axioms P192ExtensionCover.projected_target_scalar
#print axioms P192ExtensionCover.odd_power_pair_square
#print axioms P192ExtensionCover.type_i_subfield_rhs_square
#print axioms P192ExtensionCover.involution_norm_fixed
#print axioms P192ExtensionCover.involution_norm_torsion
#print axioms P192ExtensionCover.conorm_preserves_torsion
#print axioms P192ExtensionCover.mobius_quartic_clear
#print axioms P192ExtensionCover.canonicalTriples16_card
#print axioms P192ExtensionCover.three_sum_support_bound
