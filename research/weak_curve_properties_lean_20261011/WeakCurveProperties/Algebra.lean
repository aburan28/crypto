import Mathlib.Tactic.Ring
import Mathlib.Algebra.BigOperators.Group.Finset.Basic

set_option autoImplicit false

open scoped BigOperators

namespace WeakCurveProperties

variable {R : Type*} [CommRing R]

/-- The exact zero-sum-root discriminant identity retained in the catalogue. -/
theorem cubic_discriminant (r s : R) :
    (r - s)^2 * (2*r + s)^2 * (r + 2*s)^2 =
      -4 * (-(r^2 + r*s + s^2))^3 - 27 * (r*s*(r+s))^2 := by
  ring

theorem cleared_normalized_j (a b : R) :
    (6912*a^3) * 27*b^2 =
      4*a^3 * (1728*(4*a^3 + 27*b^2) - 6912*a^3) := by
  ring

theorem cleared_twist_a (a b : R) :
    (3*b)^2 * 4*a^3 = a*(2*a)^2 * 9*b^2 := by
  ring

theorem cleared_twist_b (a b : R) :
    (3*b)^3 * 8*a^3 = b*(2*a)^3 * 27*b^2 := by
  ring

theorem translated_frobenius_norm (x y c : R) :
    (x-c)*(y-c) = x*y - c*(x+y) + c^2 := by
  ring

def unitIndexPolynomial (q : R) (n : ℕ) : R :=
  ∑ i ∈ Finset.range n, q^i

/-- All eight saved unit-index identities follow from this universal theorem. -/
theorem unit_index_identity (q : R) (n : ℕ) :
    (q-1)*unitIndexPolynomial q n = q^n-1 := by
  induction n with
  | zero => simp [unitIndexPolynomial]
  | succ n ih =>
    simp only [unitIndexPolynomial, Finset.sum_range_succ] at *
    rw [mul_add, ih, pow_succ]
    ring

theorem unit_index_degree_1 (q : R) : (q-1)*unitIndexPolynomial q 1 = q^1-1 :=
  unit_index_identity q 1
theorem unit_index_degree_2 (q : R) : (q-1)*unitIndexPolynomial q 2 = q^2-1 :=
  unit_index_identity q 2
theorem unit_index_degree_3 (q : R) : (q-1)*unitIndexPolynomial q 3 = q^3-1 :=
  unit_index_identity q 3
theorem unit_index_degree_5 (q : R) : (q-1)*unitIndexPolynomial q 5 = q^5-1 :=
  unit_index_identity q 5
theorem unit_index_degree_6 (q : R) : (q-1)*unitIndexPolynomial q 6 = q^6-1 :=
  unit_index_identity q 6
theorem unit_index_degree_7 (q : R) : (q-1)*unitIndexPolynomial q 7 = q^7-1 :=
  unit_index_identity q 7
theorem unit_index_degree_10 (q : R) : (q-1)*unitIndexPolynomial q 10 = q^10-1 :=
  unit_index_identity q 10
theorem unit_index_degree_14 (q : R) : (q-1)*unitIndexPolynomial q 14 = q^14-1 :=
  unit_index_identity q 14

end WeakCurveProperties
