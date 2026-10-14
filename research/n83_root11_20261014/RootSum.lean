import Std.Tactic

/-! Bit reconstruction and capacity obligations for one retained N83 root. -/

namespace N83RootSum

def stateBound (k : Nat) : Nat := 83 * (k * (k - 1) / 2) + 42 * k

theorem xor_reconstruct (first second : Nat) :
    first ^^^ (first ^^^ second) = second := by
  rw [← Nat.xor_assoc, Nat.xor_self, Nat.zero_xor]

theorem k3400_complete_key_bound_fits_fixed_slots :
    2 * stateBound 3400 < 2 ^ 30 := by decide

theorem k3400_state_id_fits_29_bits :
    stateBound 3400 + 1 < 2 ^ 29 := by decide

theorem k3400_one_root_and_slot_payload_below_10_gibibytes :
    11 * stateBound 3400 + 5 * (2 ^ 30) < 10 * (2 ^ 30) := by decide

#print axioms xor_reconstruct
#print axioms k3400_complete_key_bound_fits_fixed_slots
#print axioms k3400_state_id_fits_29_bits
#print axioms k3400_one_root_and_slot_payload_below_10_gibibytes

end N83RootSum
