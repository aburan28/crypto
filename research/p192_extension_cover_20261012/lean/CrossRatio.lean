import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Abel
import Mathlib.Algebra.Module.LinearMap.Basic
import Mathlib.Data.Finset.Card
import Mathlib.Data.Finset.Fin
import Mathlib.Data.Finset.Prod
import Mathlib.Data.Fintype.Basic
import Mathlib.GroupTheory.OrderOfElement

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

/-- The Type I discriminant factor cannot vanish over a field in which
    minus one has no square root. In particular this covers prime fields
    of characteristic congruent to three modulo four. -/
theorem type_i_discriminant_factor_ne_zero {K : Type*} [Field K]
    (t : K) (htwo : (2 : K) ≠ 0)
    (hminus : ∀ u : K, u ^ 2 ≠ -1) :
    (t - 1) ^ 2 + 4 ≠ 0 := by
  intro hzero
  have hfour : (4 : K) ≠ 0 := by
    rw [show (4 : K) = 2 * 2 by ring]
    exact mul_ne_zero htwo htwo
  have hid : (t - 1) ^ 2 + 4 =
      4 * (((t - 1) / 2) ^ 2 + 1) := by
    field_simp [htwo]
    ring
  rw [hid] at hzero
  have hs : ((t - 1) / 2) ^ 2 + 1 = 0 :=
    (mul_eq_zero.mp hzero).resolve_left hfour
  apply hminus ((t - 1) / 2)
  calc
    ((t - 1) / 2) ^ 2 = ((t - 1) / 2) ^ 2 + 1 - 1 := by ring
    _ = -1 := by rw [hs]; ring

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

/-- A degree-two cover carries a nonzero `r`-torsion source class to a
    nonzero `r`-torsion conorm whenever its double is nonzero. The
    push-pull identity is the geometric hypothesis supplied by the cover. -/
theorem conorm_nonzero_torsion_of_push_pull_two
    {A B : Type*} [AddCommGroup A] [AddCommGroup B]
    (pull : A →+ B) (push : B →+ A) (r : ℤ) (D : A)
    (hpushpull : ∀ x, push (pull x) = x + x)
    (hrD : r • D = 0) (hdouble : D + D ≠ 0) :
    pull D ≠ 0 ∧ r • pull D = 0 := by
  constructor
  · intro hzero
    apply hdouble
    calc
      D + D = push (pull D) := (hpushpull D).symm
      _ = push 0 := by rw [hzero]
      _ = 0 := map_zero push
  · exact conorm_preserves_torsion pull r D hrD

set_option linter.style.haveILetI false

/-- For a prime-order source point, the nonzero conorm has that exact
    additive order. This specializes the preceding degree-two argument. -/
theorem conorm_exact_prime_order_of_push_pull_two
    {A B : Type*} [AddCommGroup A] [AddCommGroup B]
    (pull : A →+ B) (push : B →+ A) (r : ℕ) (D : A)
    (hpushpull : ∀ x, push (pull x) = x + x)
    (hr : r.Prime) (hrD : r • D = 0) (hdouble : D + D ≠ 0) :
    addOrderOf (pull D) = r := by
  haveI : Fact r.Prime := ⟨hr⟩
  have hne : pull D ≠ 0 := by
    intro hzero
    apply hdouble
    calc
      D + D = push (pull D) := (hpushpull D).symm
      _ = push 0 := by rw [hzero]
      _ = 0 := map_zero push
  have hmul : r • pull D = 0 := by
    simpa only [map_nsmul, map_zero] using congrArg pull hrD
  exact addOrderOf_eq_prime hmul hne

set_option linter.style.haveILetI true

/-- The three elliptic quotient coordinates determine an odd-order class
    when the Klein-four pull-push sum is doubling. The latter identity is
    the geometric hypothesis of the Type I cover. -/
theorem three_quotients_faithful_on_odd_torsion
    {J E₀ E₁ E₂ : Type*} [AddCommGroup J]
    [AddCommGroup E₀] [AddCommGroup E₁] [AddCommGroup E₂]
    (push₀ : J →+ E₀) (push₁ : J →+ E₁) (push₂ : J →+ E₂)
    (pull₀ : E₀ →+ J) (pull₁ : E₁ →+ J) (pull₂ : E₂ →+ J)
    (hdecomp : ∀ z, pull₀ (push₀ z) + pull₁ (push₁ z) +
      pull₂ (push₂ z) = z + z)
    (r k : ℕ) (hodd : r = 2 * k + 1)
    (x y : J) (hx : r • x = 0) (hy : r • y = 0)
    (h₀ : push₀ x = push₀ y) (h₁ : push₁ x = push₁ y)
    (h₂ : push₂ x = push₂ y) : x = y := by
  let z := x - y
  have hz₀ : push₀ z = 0 := by simp [z, map_sub, h₀]
  have hz₁ : push₁ z = 0 := by simp [z, map_sub, h₁]
  have hz₂ : push₂ z = 0 := by simp [z, map_sub, h₂]
  have hdouble : z + z = 0 := by
    calc
      z + z = pull₀ (push₀ z) + pull₁ (push₁ z) +
          pull₂ (push₂ z) := (hdecomp z).symm
      _ = 0 := by simp [hz₀, hz₁, hz₂]
  have hrz : r • z = 0 := by simp [z, nsmul_sub, hx, hy]
  have h2k : (2 * k) • z = 0 := by
    calc
      (2 * k) • z = (k + k) • z := by congr 1; omega
      _ = k • z + k • z := by rw [add_nsmul]
      _ = k • (z + z) := by rw [smul_add]
      _ = 0 := by rw [hdouble, nsmul_zero]
  have hrz' : (2 * k + 1) • z = 0 := by simpa [hodd] using hrz
  have hrepr : (2 * k + 1) • z = z := by
    rw [add_nsmul, h2k, one_nsmul, zero_add]
  rw [hrepr] at hrz'
  exact sub_eq_zero.mp hrz'

/-- On odd-order torsion, the three quotient pullbacks give an explicit
    inverse to the quotient-coordinate map: half their sum is the class. -/
theorem three_quotients_inverse_on_odd_torsion
    {J E₀ E₁ E₂ : Type*} [AddCommGroup J]
    [AddCommGroup E₀] [AddCommGroup E₁] [AddCommGroup E₂]
    (push₀ : J →+ E₀) (push₁ : J →+ E₁) (push₂ : J →+ E₂)
    (pull₀ : E₀ →+ J) (pull₁ : E₁ →+ J) (pull₂ : E₂ →+ J)
    (hdecomp : ∀ z, pull₀ (push₀ z) + pull₁ (push₁ z) +
      pull₂ (push₂ z) = z + z)
    (r k : ℕ) (hodd : r = 2 * k + 1)
    (D : J) (hrD : r • D = 0) :
    (k + 1) • (pull₀ (push₀ D) + pull₁ (push₁ D) +
      pull₂ (push₂ D)) = D := by
  rw [hdecomp]
  have hcoeff : (k + 1) * 2 = r + 1 := by omega
  calc
    (k + 1) • (D + D) = ((k + 1) * 2) • D := by
      simp only [mul_nsmul, two_nsmul, smul_add]
    _ = (r + 1) • D := by rw [hcoeff]
    _ = D := by rw [add_nsmul, hrD, one_nsmul, zero_add]

set_option linter.style.haveILetI false

/-- A three-step Frobenius norm is fixed and has exact prime order when
    one elliptic quotient sees a nonzero component in precisely one
    translate.  The concrete quotient coordinates are checked by Sage. -/
theorem cyclic_three_norm_exact_prime_order
    {J E : Type*} [AddCommGroup J] [AddCommGroup E]
    (σ : J →+ J) (f : J →+ E) (r : ℕ) (D : J) (U : E)
    (hr : r.Prime) (hrD : r • D = 0)
    (hσ3 : σ (σ (σ D)) = D)
    (h0 : f D = U) (h1 : f (σ D) = 0)
    (h2 : f (σ (σ D)) = 0) (hU : U ≠ 0) :
    let N := D + σ D + σ (σ D)
    σ N = N ∧ f N = U ∧ addOrderOf N = r := by
  let N := D + σ D + σ (σ D)
  have hfixed : σ N = N := by
    dsimp [N]
    rw [map_add, map_add, hσ3]
    ac_rfl
  have hproj : f N = U := by
    simp [N, map_add, h0, h1, h2]
  have hne : N ≠ 0 := by
    intro hzero
    apply hU
    calc
      U = f N := hproj.symm
      _ = f 0 := by rw [hzero]
      _ = 0 := map_zero f
  have hσD : r • σ D = 0 := by
    simpa only [map_nsmul, map_zero] using congrArg σ hrD
  have hσ2D : r • σ (σ D) = 0 := by
    simpa only [map_nsmul, map_zero] using congrArg σ hσD
  have hmul : r • N = 0 := by
    simp [N, nsmul_add, hrD, hσD, hσ2D]
  haveI : Fact r.Prime := ⟨hr⟩
  exact ⟨hfixed, hproj, addOrderOf_eq_prime hmul hne⟩

set_option linter.style.haveILetI true

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

/-- Conditioning queries on any geometric predicate cannot increase the
    number of covered classes in a fixed 16-point three-sum base. -/
theorem conditioned_three_sum_support_bound {A : Type*}
    [AddCommMonoid A] [DecidableEq A]
    (base : Fin 16 → A) (admissible : Finset A) :
    (admissible ∩ canonicalTriples16.image
      (fun t => base t.1 + base t.2.1 + base t.2.2)).card ≤ 816 := by
  calc
    (admissible ∩ canonicalTriples16.image
      (fun t => base t.1 + base t.2.1 + base t.2.2)).card ≤
      (canonicalTriples16.image
        (fun t => base t.1 + base t.2.1 + base t.2.2)).card :=
          Finset.card_le_card Finset.inter_subset_right
    _ ≤ 816 := three_sum_support_bound base

/-- Character-sum bounds from genus 7 and 31 double covers imply a
    fourfold lower bound on the accepted scalar count. The geometric
    Riemann--Hurwitz and Hasse--Weil hypotheses are explicit inputs here. -/
theorem character_filter_acceptance_lower_bound
    (r U T accepted s₁ s₁₂ identityAccepted : ℤ)
    (hcount : 4 * (accepted + identityAccepted) =
      r + 2 * s₁ + s₁₂)
    (h₁ : -(14 * U + T) ≤ s₁)
    (h₁₂ : -(62 * U + T) ≤ s₁₂) :
    r - 90 * U - 3 * T - 4 * identityAccepted ≤ 4 * accepted := by
  omega

/-- Exact integer specialization of the P-192 character-density bound.
    The geometric character-sum inequalities are supplied separately. -/
theorem p192_character_filter_integer_certificate :
    ((6277101735386680763835789423176059013767194773182842284081 : ℕ) -
      90 * 79228162514264337593543950336 -
      3 * 31607402316713927207482677199 + 3) / 4 - 1 =
      1569275433846670190958947353987675545133315652035359680560 ∧
    816 * (2 : ℕ)^180 <
      1569275433846670190958947353987675545133315652035359680560 := by
  constructor <;> norm_num

/-- A fixed shift of a 16-point three-sum base still has at most 816 images. -/
theorem shifted_three_sum_support_bound {A : Type*} [AddCommGroup A] [DecidableEq A]
    (base : Fin 16 → A) (shift : A) :
    (canonicalTriples16.image
      (fun t => base t.1 + base t.2.1 + base t.2.2 - shift)).card ≤ 816 := by
  calc
    (canonicalTriples16.image
        (fun t => base t.1 + base t.2.1 + base t.2.2 - shift)).card ≤
        canonicalTriples16.card := Finset.card_image_le
    _ = 816 := canonicalTriples16_card

/-- A known-scalar shifted query yields a usable base-log equation on a hit. -/
theorem shifted_three_point_relation {A : Type*} [AddCommGroup A]
    (a : ℤ) (N U V X Y Z : A)
    (h : a • N + U + V = X + Y + Z) :
    a • N = X + Y + Z - U - V := by
  calc
    a • N = (a • N + U + V) - U - V := by abel
    _ = X + Y + Z - U - V := by rw [h]

/-- A six-periodic eigenvector cannot have eigenvalue whose sixth power is
    different from one. Applied to elliptic r-torsion, this isolates the
    Frobenius-fixed line when the other eigenvalue is p and p^6 ≠ 1 mod r. -/
theorem six_periodic_eigenvector_zero {K A : Type*} [Field K]
    [AddCommGroup A] [Module K A] (F : A →ₗ[K] A) (p : K) (y : A)
    (heigen : F y = p • y)
    (hperiod : F (F (F (F (F (F y))))) = y)
    (hp : p ^ 6 ≠ 1) : y = 0 := by
  have h2 : F (F y) = p ^ 2 • y := by
    rw [heigen, map_smul, heigen]
    simp [pow_two, smul_smul]
  have h3 : F (F (F y)) = p ^ 3 • y := by
    rw [h2, map_smul, heigen]
    simp [pow_succ, smul_smul, pow_two, mul_assoc]
  have h4 : F (F (F (F y))) = p ^ 4 • y := by
    rw [h3, map_smul, heigen]
    simp [pow_succ, smul_smul, mul_assoc]
  have h5 : F (F (F (F (F y)))) = p ^ 5 • y := by
    rw [h4, map_smul, heigen]
    simp [pow_succ, smul_smul, mul_assoc]
  have h6 : F (F (F (F (F (F y))))) = p ^ 6 • y := by
    rw [h5, map_smul, heigen]
    simp [pow_succ, smul_smul, mul_assoc]
  have hsmul : (p ^ 6 - 1) • y = 0 := by
    rw [sub_smul, one_smul]
    rw [← h6, hperiod]
    simp
  have hnonzero : p ^ 6 - 1 ≠ 0 := sub_ne_zero.mpr hp
  calc
    y = (p ^ 6 - 1)⁻¹ • ((p ^ 6 - 1) • y) := by
      simp [smul_smul, hnonzero]
    _ = 0 := by rw [hsmul, smul_zero]

/- The primitive coordinate s = y₀+y₁+y₂ of the degree-four Type I cover
   satisfies this quartic when y₀y₁y₂ is the chosen norm component. -/
theorem type_i_primitive_quartic_identity {K : Type*} [CommRing K]
    (y₀ y₁ y₂ : K) :
    let s := y₀ + y₁ + y₂
    let A := y₀^2 + y₁^2 + y₂^2
    let B := y₀^2*y₁^2 + y₀^2*y₂^2 + y₁^2*y₂^2
    let H := y₀*y₁*y₂
    s^4 - 2*A*s^2 - 8*H*s + A^2 - 4*B = 0 := by
  dsimp
  ring

/-- Eliminating the base coordinate from the two quadratic quotient equations
    factors the Sylvester resultant as w^4 times a Ciani quartic. -/
theorem type_i_ciani_resultant_factor {K : Type*} [CommRing K]
    (s₀ s₁ s₂ n₀ n₁ n₂ u v w : K) :
    let a₂ := u^2-w^2
    let a₁ := -s₀*u^2+s₁*w^2
    let a₀ := n₀*u^2-n₁*w^2
    let b₂ := v^2-w^2
    let b₁ := -s₀*v^2+s₂*w^2
    let b₀ := n₀*v^2-n₂*w^2
    let L₀ := (n₀-n₂)*u^2+(n₁-n₀)*v^2+(n₂-n₁)*w^2
    let L₁ := (s₂-s₀)*u^2+(s₀-s₁)*v^2+(s₁-s₂)*w^2
    let L₂ := (s₀*n₂-n₀*s₂)*u^2+(s₁*n₀-n₁*s₀)*v^2+
      (n₁*s₂-s₁*n₂)*w^2
    (a₂*b₀-a₀*b₂)^2-(a₂*b₁-a₁*b₂)*(a₁*b₀-a₀*b₁) =
      w^4*(L₀^2-L₁*L₂) := by
  dsimp
  ring

/-- The exact neighboring capacities certify the half-coverage threshold for
    the fixed-base, unordered three-sum policy. Monotonicity of
    `B * (B+1) * (B+2)` supplies the general integer threshold. -/
theorem p192_conditioned_half_coverage_threshold_certificate :
    2 * ((16759979263520929693 : ℕ) * 16759979263520929694 *
      16759979263520929695 / 6) <
      1569275433846670190958947353987675545133315652035359680560 ∧
    1569275433846670190958947353987675545133315652035359680560 ≤
      2 * ((16759979263520929694 : ℕ) * 16759979263520929695 *
        16759979263520929696 / 6) := by
  constructor <;> norm_num

/-- A full unordered pair index at the half-coverage threshold has this many
    entries, before implementation overhead or collision handling. -/
theorem p192_threshold_pair_index_entries :
    (16759979263520929694 : ℕ) * 16759979263520929695 / 2 =
      140448452456825782461936609026805931665 := by
  norm_num

end P192ExtensionCover

#print axioms P192ExtensionCover.cyclic_cross_ratio_product
#print axioms P192ExtensionCover.cyclic_cross_ratio_not_norm_one
#print axioms P192ExtensionCover.type_i_discriminant_factor_ne_zero
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
#print axioms P192ExtensionCover.conorm_nonzero_torsion_of_push_pull_two
#print axioms P192ExtensionCover.conorm_exact_prime_order_of_push_pull_two
#print axioms P192ExtensionCover.cyclic_three_norm_exact_prime_order
#print axioms P192ExtensionCover.three_quotients_faithful_on_odd_torsion
#print axioms P192ExtensionCover.three_quotients_inverse_on_odd_torsion
#print axioms P192ExtensionCover.mobius_quartic_clear
#print axioms P192ExtensionCover.canonicalTriples16_card
#print axioms P192ExtensionCover.three_sum_support_bound
#print axioms P192ExtensionCover.conditioned_three_sum_support_bound
#print axioms P192ExtensionCover.character_filter_acceptance_lower_bound
#print axioms P192ExtensionCover.p192_character_filter_integer_certificate
#print axioms P192ExtensionCover.shifted_three_sum_support_bound
#print axioms P192ExtensionCover.shifted_three_point_relation
#print axioms P192ExtensionCover.six_periodic_eigenvector_zero
#print axioms P192ExtensionCover.type_i_primitive_quartic_identity
#print axioms P192ExtensionCover.type_i_ciani_resultant_factor
#print axioms P192ExtensionCover.p192_conditioned_half_coverage_threshold_certificate
#print axioms P192ExtensionCover.p192_threshold_pair_index_entries
