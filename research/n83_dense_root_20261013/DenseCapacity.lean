/- Numeric capacity bounds for the opt-in dense packed root table.
   The Rust source requests one entry per regular state and RootTable128
   allocates the next power of two at or above twice that request.
   The source loop has at most 83*K*K states for N83; this module proves
   arithmetic implications under that bound. The loop-to-bound bridge
   remains a source inspection and measured-state obligation. -/

namespace N83DenseRootCapacity

theorem k1271_requested_slots_fit_four_gib (states : Nat)
    (hStates : states ≤ 83 * 1271 ^ 2) :
    2 * states ≤ 2 ^ 28 := by
  calc
    2 * states ≤ 2 * (83 * 1271 ^ 2) := Nat.mul_le_mul_left 2 hStates
    _ ≤ 2 ^ 28 := by decide

theorem k1272_maximum_crosses_four_gib :
    2 ^ 28 < 2 * (83 * 1272 ^ 2) := by decide

theorem k1798_requested_slots_fit_eight_gib (states : Nat)
    (hStates : states ≤ 83 * 1798 ^ 2) :
    2 * states ≤ 2 ^ 29 := by
  calc
    2 * states ≤ 2 * (83 * 1798 ^ 2) := Nat.mul_le_mul_left 2 hStates
    _ ≤ 2 ^ 29 := by decide

theorem k1799_maximum_crosses_eight_gib :
    2 ^ 29 < 2 * (83 * 1799 ^ 2) := by decide

theorem four_gib_packed_slot_bytes : 16 * 2 ^ 28 = 4294967296 := by decide

theorem eight_gib_packed_slot_bytes : 16 * 2 ^ 29 = 8589934592 := by decide

#print axioms k1271_requested_slots_fit_four_gib
#print axioms k1272_maximum_crosses_four_gib
#print axioms k1798_requested_slots_fit_eight_gib
#print axioms k1799_maximum_crosses_eight_gib
#print axioms four_gib_packed_slot_bytes
#print axioms eight_gib_packed_slot_bytes

end N83DenseRootCapacity
