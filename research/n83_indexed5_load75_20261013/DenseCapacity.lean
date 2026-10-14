/-!
Exact source-loop bounds for the optional 75-percent request on the N83
symmetric S3 five-byte index. The table constructor rounds twice the request
to a power of two, so the 2^29-slot bucket ends at K=3114.
-/

namespace N83Indexed5Load75

def stateBound (k : Nat) : Nat :=
  83 * (k * (k - 1) / 2) + 42 * k

def requestEntries (k : Nat) : Nat :=
  (2 * stateBound k + 2) / 3

theorem k3114_request_fits_2pow29_slots :
    2 * requestEntries 3114 ≤ 2 ^ 29 := by decide

theorem k3115_request_crosses_2pow29_slots :
    2 ^ 29 < 2 * requestEntries 3115 := by decide

theorem k3114_state_id_fits_29_bits :
    stateBound 3114 + 1 < 2 ^ 29 := by decide

theorem k3114_state_and_slot_payload_under_13_gibibytes :
    28 * stateBound 3114 + 5 * (2 ^ 29) <
      13 * (2 ^ 30) := by decide

#print axioms k3114_request_fits_2pow29_slots
#print axioms k3115_request_crosses_2pow29_slots
#print axioms k3114_state_id_fits_29_bits
#print axioms k3114_state_and_slot_payload_under_13_gibibytes

end N83Indexed5Load75
