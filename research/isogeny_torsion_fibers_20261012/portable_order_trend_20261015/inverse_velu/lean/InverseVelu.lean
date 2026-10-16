import Lean

/-!
The integer polynomial identity behind the selected degree-three inverse
Vélu x-coordinate. This module does not define the finite field, curve,
or Vélu morphism used by the Rust experiment; that formal bridge is open.
-/

namespace InverseVelu

def cubic (x u : Int) : Int :=
  x*x*x - (u+6)*x*x + (6*u+21)*x - (9*u+32)

theorem cubicIdentity (x u : Int) :
    cubic x u = -((u+x+12)*(x-3)*(x-3)-2*(x*x*x-21*x+38)) := by
  unfold cubic
  simp only [Int.add_mul, Int.mul_add, Int.sub_mul, Int.mul_sub]
  simp only [Int.mul_assoc]
  simp [Int.mul_comm, Int.mul_left_comm, Int.mul_assoc]
  have h3 : x * (x * x) = x*x*x := by ac_rfl
  have h12 : x * (x * 12) = 12 * (x*x) := by ac_rfl
  have h3x : x * (x * 3) = 3 * (x*x) := by ac_rfl
  have h6x : x * (x * 6) = 6 * (x*x) := by ac_rfl
  have hu3 : x * (u * 3) = 3 * (u*x) := by ac_rfl
  have hu6 : x * (u * 6) = 6 * (u*x) := by ac_rfl
  have h2xxx : x * (x * (x * 2)) = 2 * (x*x*x) := by ac_rfl
  simp only [h3, h12, h3x, h6x, hu3, hu6, h2xxx]
  omega

theorem cubicOfCrossMultiplication (x u : Int)
    (h : (u+x+12)*(x-3)*(x-3) = 2*(x*x*x-21*x+38)) :
    cubic x u = 0 := by
  rw [cubicIdentity]
  omega

theorem cubicOfSourceOrdinate (x y u : Int)
    (curve : y*y = x*x*x-21*x+37)
    (mapX : (u+x+12)*(x-3)*(x-3) = 2*(y*y+1)) :
    cubic x u = 0 := by
  have h : (u+x+12)*(x-3)*(x-3) = 2*(x*x*x-21*x+38) := by
    rw [curve] at mapX
    omega
  exact cubicOfCrossMultiplication x u h

end InverseVelu

#print axioms InverseVelu.cubicIdentity
#print axioms InverseVelu.cubicOfCrossMultiplication
#print axioms InverseVelu.cubicOfSourceOrdinate
