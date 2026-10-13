/- Numeric capacity bounds for a proposed 75%-load packed-root table.
   A proposed request of ceil(2*states/3) passed to RootTable128's existing
   2x power-of-two allocator would reserve at least four slots per three
   states. This module proves capacity implications only; the proposed
   request rule is not implemented in the measured producer. -/

namespace N83Load75Capacity

theorem k2202_request_fits_eight_gib (states : Nat)
    (hStates : states ≤ 83 * 2202 ^ 2) :
    2 * ((2 * states + 2) / 3) ≤ 2 ^ 29 := by
  have hScaled : 2 * states + 2 ≤ 2 * (83 * 2202 ^ 2) + 2 :=
    Nat.add_le_add_right (Nat.mul_le_mul_left 2 hStates) 2
  calc
    2 * ((2 * states + 2) / 3)
        ≤ 2 * ((2 * (83 * 2202 ^ 2) + 2) / 3) :=
          Nat.mul_le_mul_left 2 (Nat.div_le_div_right hScaled)
    _ ≤ 2 ^ 29 := by decide

theorem k2203_maximum_crosses_eight_gib :
    2 ^ 29 < 2 * ((2 * (83 * 2203 ^ 2) + 2) / 3) := by decide

#print axioms k2202_request_fits_eight_gib
#print axioms k2203_maximum_crosses_eight_gib

end N83Load75Capacity
