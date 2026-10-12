import WeakCurveProperties.Conductor
import WeakCurveProperties.Torsion

set_option autoImplicit false

namespace WeakCurveProperties

/-- Conditional application to a Frobenius action on geometric points.
The two coefficient equivalences and the CM order-index identity are exposed
as hypotheses. Their interpretation for the actual ordinary elliptic-curve
endomorphism ring is the remaining geometric dependency in the proof map. -/
theorem torsion_order_depth_bounds_of_bridges {A : Type*} [AddCommGroup A]
    (Frob : A →+ A) {fE fPi : ℕ} {traceOmega normOmega a b : ℤ}
    (hE : fE ≠ 0) (hb : b ≠ 0) (hindex : fPi = fE*b.natAbs)
    (hodd : ¬(2 : ℤ) ∣ quadraticNorm traceOmega normOmega a b)
    (hcoeff2 : fixesNTorsion Frob 2 ↔ coefficientTorsion 2 a b)
    (hcoeff4 : fixesNTorsion Frob 4 ↔ coefficientTorsion 4 a b) :
    (¬fixesNTorsion Frob 2 → padicValNat 2 fE = padicValNat 2 fPi) ∧
    (fixesNTorsion Frob 2 → padicValNat 2 fE+1 ≤ padicValNat 2 fPi) ∧
    (fixesNTorsion Frob 4 → padicValNat 2 fE+2 ≤ padicValNat 2 fPi) := by
  refine ⟨?_, ?_, ?_⟩
  · intro hnot
    exact conductor_floor hE hb hindex hodd (fun h => hnot (hcoeff2.mpr h))
  · intro hfull
    exact conductor_full_two hE hb hindex (hcoeff2.mp hfull)
  · intro hfull
    exact conductor_full_four hE hb hindex (hcoeff4.mp hfull)

end WeakCurveProperties
