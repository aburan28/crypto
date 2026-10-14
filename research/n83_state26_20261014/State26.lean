import Std.Tactic

/-! The 31-bit label and 26-byte payload used by the optional N83 state store. -/

namespace N83State26

def encodeLabels (left right relative : Nat) : Nat :=
  left + 2 ^ 12 * right + 2 ^ 24 * relative

theorem labels_roundtrip
    (left right relative : Nat)
    (hl : left < 2 ^ 12) (hr : right < 2 ^ 12)
    (ht : relative < 2 ^ 7) :
    encodeLabels left right relative < 2 ^ 31 ∧
    encodeLabels left right relative % (2 ^ 12) = left ∧
    encodeLabels left right relative / (2 ^ 12) % (2 ^ 12) = right ∧
    encodeLabels left right relative / (2 ^ 24) % (2 ^ 7) = relative := by
  unfold encodeLabels
  omega

def stateBound (k : Nat) : Nat := 83 * (k * (k - 1) / 2) + 42 * k
def requestedEntries (k : Nat) : Nat := (4 * stateBound k + 2) / 3

theorem k2940_slot_request_fits :
    2 * requestedEntries 2940 ≤ 2 ^ 30 := by decide

theorem k2940_state_indices_fit :
    stateBound 2940 + 1 < 2 ^ 29 := by decide

theorem k2940_payload_under_13point7_gibibytes :
    10 * (26 * stateBound 2940 + 5 * (2 ^ 30)) < 137 * (2 ^ 30) := by decide

#print axioms labels_roundtrip
#print axioms k2940_slot_request_fits
#print axioms k2940_state_indices_fit
#print axioms k2940_payload_under_13point7_gibibytes

end N83State26
