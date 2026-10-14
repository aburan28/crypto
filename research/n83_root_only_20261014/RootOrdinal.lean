import Std.Tactic

/-! Arithmetic of the guarded N83 symmetric state-ordinal layout. -/

namespace N83RootOrdinal

def stateBound (k : Nat) : Nat := 83 * (k * (k - 1) / 2) + 42 * k
def requestedEntries (k : Nat) : Nat := (4 * stateBound k + 2) / 3

theorem diagonal_offset_roundtrip (pfx relative : Nat) (h : relative < 42) :
    pfx + relative < pfx + 42 ∧
    pfx + relative - pfx = relative := by omega

theorem off_diagonal_offset_roundtrip
    (rightOffset relative : Nat) (h : relative < 83) :
    (42 + 83 * rightOffset + relative - 42) / 83 = rightOffset ∧
    (42 + 83 * rightOffset + relative - 42) % 83 = relative := by omega

theorem k3114_request_fits_2pow30_slots :
    2 * requestedEntries 3114 ≤ 2 ^ 30 := by decide

theorem k3115_request_crosses_2pow30_slots :
    2 ^ 30 < 2 * requestedEntries 3115 := by decide

theorem k3114_state_id_fits_29_bits :
    stateBound 3114 + 1 < 2 ^ 29 := by decide

theorem k3114_root_and_slot_payload_below_13point3_gibibytes :
    10 * (22 * stateBound 3114 + 5 * (2 ^ 30)) < 133 * (2 ^ 30) := by decide

#print axioms diagonal_offset_roundtrip
#print axioms off_diagonal_offset_roundtrip
#print axioms k3114_request_fits_2pow30_slots
#print axioms k3115_request_crosses_2pow30_slots
#print axioms k3114_state_id_fits_29_bits
#print axioms k3114_root_and_slot_payload_below_13point3_gibibytes

end N83RootOrdinal
