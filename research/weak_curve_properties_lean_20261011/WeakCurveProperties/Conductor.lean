import Mathlib.NumberTheory.Padics.PadicVal.Basic
import Mathlib.Tactic.Ring

set_option autoImplicit false

namespace WeakCurveProperties

/-- Norm of `a+b*omega` in an integral quadratic basis. -/
def quadraticNorm (traceOmega normOmega a b : ℤ) : ℤ :=
  a^2 + traceOmega*a*b + normOmega*b^2

/-- Coefficient divisibility in the lattice with basis `1,omega`. -/
def coefficientTorsion (n a b : ℤ) : Prop := n ∣ a-1 ∧ n ∣ b

theorem odd_norm_forces_odd_scalar {traceOmega normOmega a b : ℤ}
    (hodd : ¬(2 : ℤ) ∣ quadraticNorm traceOmega normOmega a b)
    (hb : (2 : ℤ) ∣ b) : ¬(2 : ℤ) ∣ a := by
  intro ha
  rcases ha with ⟨x, rfl⟩
  rcases hb with ⟨y, rfl⟩
  apply hodd
  refine ⟨2*x^2 + 2*traceOmega*x*y + 2*normOmega*y^2, ?_⟩
  unfold quadraticNorm
  ring

theorem two_divides_scalar_minus_one {a : ℤ} (ha : ¬(2 : ℤ) ∣ a) :
    (2 : ℤ) ∣ a-1 := by
  rw [Int.dvd_iff_emod_eq_zero] at *
  omega

/-- The odd norm makes the two coefficient conditions equivalent to `2 | b`. -/
theorem coefficient_two_iff {traceOmega normOmega a b : ℤ}
    (hodd : ¬(2 : ℤ) ∣ quadraticNorm traceOmega normOmega a b) :
    coefficientTorsion 2 a b ↔ (2 : ℤ) ∣ b := by
  constructor
  · exact And.right
  · intro hb
    exact ⟨two_divides_scalar_minus_one (odd_norm_forces_odd_scalar hodd hb), hb⟩

/-- The valuation gain is a consequence of the order index and actual divisibility. -/
theorem conductor_gain {fE fPi : ℕ} {b : ℤ} (hE : fE ≠ 0) (hb : b ≠ 0)
    (hindex : fPi = fE*b.natAbs) (k : ℕ) (hdiv : (2 : ℤ)^k ∣ b) :
    padicValNat 2 fE + k ≤ padicValNat 2 fPi := by
  let : Fact (Nat.Prime 2) := ⟨Nat.prime_two⟩
  have hlevel : k ≤ padicValInt 2 b :=
    ((padicValInt_dvd_iff k b).mp hdiv).resolve_left hb
  rw [hindex, padicValNat.mul hE (Int.natAbs_ne_zero.mpr hb)]
  exact Nat.add_le_add_left hlevel _

/-- Failure of full coefficient 2-torsion certifies the floor of this order model. -/
theorem conductor_floor {fE fPi : ℕ} {traceOmega normOmega a b : ℤ}
    (hE : fE ≠ 0) (hb : b ≠ 0) (hindex : fPi = fE*b.natAbs)
    (hodd : ¬(2 : ℤ) ∣ quadraticNorm traceOmega normOmega a b)
    (hnot : ¬coefficientTorsion 2 a b) :
    padicValNat 2 fE = padicValNat 2 fPi := by
  let : Fact (Nat.Prime 2) := ⟨Nat.prime_two⟩
  have hnotb : ¬(2 : ℤ) ∣ b := fun h => hnot ((coefficient_two_iff hodd).mpr h)
  have hz : padicValNat 2 b.natAbs = 0 :=
    padicValInt.eq_zero_of_not_dvd hnotb
  rw [hindex, padicValNat.mul hE (Int.natAbs_ne_zero.mpr hb), hz, Nat.add_zero]

theorem conductor_full_two {fE fPi : ℕ} {a b : ℤ}
    (hE : fE ≠ 0) (hb : b ≠ 0) (hindex : fPi = fE*b.natAbs)
    (hfull : coefficientTorsion 2 a b) :
    padicValNat 2 fE + 1 ≤ padicValNat 2 fPi := by
  exact conductor_gain hE hb hindex 1 (by simpa using hfull.2)

theorem conductor_full_four {fE fPi : ℕ} {a b : ℤ}
    (hE : fE ≠ 0) (hb : b ≠ 0) (hindex : fPi = fE*b.natAbs)
    (hfull : coefficientTorsion 4 a b) :
    padicValNat 2 fE + 2 ≤ padicValNat 2 fPi := by
  exact conductor_gain hE hb hindex 2 (by norm_num at *; exact hfull.2)

/-- The surface conclusion is exact when the proved bound collapses. -/
theorem conductor_surface_two {fE fPi : ℕ} {a b : ℤ}
    (hE : fE ≠ 0) (hb : b ≠ 0) (hindex : fPi = fE*b.natAbs)
    (hfull : coefficientTorsion 2 a b) (hdepth : padicValNat 2 fPi = 1) :
    padicValNat 2 fE = 0 := by
  have := conductor_full_two hE hb hindex hfull
  omega

theorem conductor_surface_four {fE fPi : ℕ} {a b : ℤ}
    (hE : fE ≠ 0) (hb : b ≠ 0) (hindex : fPi = fE*b.natAbs)
    (hfull : coefficientTorsion 4 a b) (hdepth : padicValNat 2 fPi = 2) :
    padicValNat 2 fE = 0 := by
  have := conductor_full_four hE hb hindex hfull
  omega

end WeakCurveProperties
