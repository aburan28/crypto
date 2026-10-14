import Std.Tactic

/-! Exact arithmetic obligations for the N83 four-byte indexed root slot. -/

namespace N83RootFour

def stateBound (k : Nat) : Nat := 83 * (k * (k - 1) / 2) + 42 * k

/-- The source stores `state + 1` in bits 1..29, a root selector in bit 0,
    and a two-bit fingerprint in bits 30..31. Zero is the empty slot. -/
def encode (state selector fingerprint : Nat) : Nat :=
  2 ^ 30 * fingerprint + 2 * (state + 1) + selector

theorem decode_selector (state selector fingerprint : Nat)
    (hs : selector < 2) : encode state selector fingerprint % 2 = selector := by
  unfold encode
  omega

theorem decode_state (state selector fingerprint : Nat)
    (hstate : state + 1 < 2 ^ 29) (hs : selector < 2) :
    (encode state selector fingerprint % (2 ^ 30)) / 2 - 1 = state := by
  unfold encode
  omega

theorem decode_fingerprint (state selector fingerprint : Nat)
    (hstate : state + 1 < 2 ^ 29) (hs : selector < 2) :
    encode state selector fingerprint / (2 ^ 30) = fingerprint := by
  unfold encode
  omega

theorem encode_nonzero (state selector fingerprint : Nat) :
    encode state selector fingerprint ≠ 0 := by
  unfold encode
  omega

theorem encode_fits_u32 (state selector fingerprint : Nat)
    (hstate : state + 1 < 2 ^ 29) (hs : selector < 2)
    (hfp : fingerprint < 4) :
    encode state selector fingerprint < 2 ^ 32 := by
  unfold encode
  omega

theorem k3550_complete_key_bound_fits_fixed_slots :
    2 * stateBound 3550 < 2 ^ 31 := by decide

theorem k3550_state_id_fits_29_bits :
    stateBound 3550 + 1 < 2 ^ 29 := by decide

theorem k3550_one_root_and_slot_payload_below_14_gibibytes :
    11 * stateBound 3550 + 4 * (2 ^ 31) < 14 * (2 ^ 30) := by decide

theorem k3596_state_id_fits_29_bits :
    stateBound 3596 + 1 < 2 ^ 29 := by decide

theorem k3597_state_id_exceeds_29_bits :
    2 ^ 29 ≤ stateBound 3597 + 1 := by decide

theorem k3596_one_root_and_slot_payload_below_14_gibibytes :
    11 * stateBound 3596 + 4 * (2 ^ 31) < 14 * (2 ^ 30) := by decide

#print axioms decode_selector
#print axioms decode_state
#print axioms decode_fingerprint
#print axioms encode_nonzero
#print axioms encode_fits_u32
#print axioms k3550_complete_key_bound_fits_fixed_slots
#print axioms k3550_state_id_fits_29_bits
#print axioms k3550_one_root_and_slot_payload_below_14_gibibytes
#print axioms k3596_state_id_fits_29_bits
#print axioms k3597_state_id_exceeds_29_bits
#print axioms k3596_one_root_and_slot_payload_below_14_gibibytes

end N83RootFour
