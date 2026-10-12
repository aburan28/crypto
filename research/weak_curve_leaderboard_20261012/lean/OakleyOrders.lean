import Mathlib.GroupTheory.OrderOfElement

namespace WeakCurveLeaderboard

def oakley3Prime : ℕ := 3805993847215893016155463826195386266397436443
def oakley4Prime : ℕ := 12259964326927110866866776214413170562013096250261263279

theorem oakley3_full_order_factor :
    45671926166590716193865565914344635196769237316 = 12 * oakley3Prime := by
  decide

theorem oakley4_full_order_factor :
    49039857307708443467467104857652682248052385001045053116 =
      4 * oakley4Prime := by
  decide

theorem oakley3_extension_factor : 155 = 5 * 31 := by decide
theorem oakley4_extension_factor : 185 = 5 * 37 := by decide

/-- A published generator of order c*r yields a point of exact order r after
    multiplying by c. This is the structural subgroup step used for Oakley. -/
theorem selected_order_of_generator {G : Type*} [AddGroup G]
    (P : G) (r c : ℕ) (hc : c ≠ 0) (hP : addOrderOf P = c * r) :
    addOrderOf (c • P) = r := by
  have hdiv : c ∣ addOrderOf P := by
    rw [hP]
    exact dvd_mul_right c r
  rw [addOrderOf_nsmul_of_dvd hc hdiv, hP]
  exact Nat.mul_div_cancel_left r (Nat.pos_of_ne_zero hc)

theorem oakley3_selected_order {G : Type*} [AddGroup G] (P : G)
    (hP : addOrderOf P = 4 * oakley3Prime) :
    addOrderOf (4 • P) = oakley3Prime := by
  exact selected_order_of_generator P oakley3Prime 4 (by decide) hP

theorem oakley4_selected_order {G : Type*} [AddGroup G] (P : G)
    (hP : addOrderOf P = 2 * oakley4Prime) :
    addOrderOf (2 • P) = oakley4Prime := by
  exact selected_order_of_generator P oakley4Prime 2 (by decide) hP

end WeakCurveLeaderboard
