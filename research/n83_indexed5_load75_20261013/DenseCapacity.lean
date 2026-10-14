/-!
Exact source-loop bounds for the optional 75-percent request on the N83
symmetric S3 five-byte index. The symmetric quotient has almost two distinct
root keys per retained state, so the request is ceil(4*S/3). The table
constructor rounds twice that request to a power of two.
-/

namespace N83Indexed5Load75

def stateBound (k : Nat) : Nat :=
  83 * (k * (k - 1) / 2) + 42 * k

def requestEntries (k : Nat) : Nat :=
  (4 * stateBound k + 2) / 3

theorem k2830_request_fits_2pow30_slots :
    2 * requestEntries 2830 ≤ 2 ^ 30 := by decide

theorem k3114_request_fits_2pow30_slots :
    2 * requestEntries 3114 ≤ 2 ^ 30 := by decide

theorem k3115_request_crosses_2pow30_slots :
    2 ^ 30 < 2 * requestEntries 3115 := by decide

theorem k2830_state_id_fits_29_bits :
    stateBound 2830 + 1 < 2 ^ 29 := by decide

theorem k2830_state_and_slot_payload_under_13point7_gibibytes :
    10 * (28 * stateBound 2830 + 5 * (2 ^ 30)) <
      137 * (2 ^ 30) := by decide

#print axioms k2830_request_fits_2pow30_slots
#print axioms k3114_request_fits_2pow30_slots
#print axioms k3115_request_crosses_2pow30_slots
#print axioms k2830_state_id_fits_29_bits
#print axioms k2830_state_and_slot_payload_under_13point7_gibibytes

end N83Indexed5Load75
