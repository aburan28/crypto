import Lean

/-!
Exact state-count ceiling for the opt-in N83 five-byte root slot with a
30-bit state ID, one root-selector bit, and nine fingerprint bits. The
implementation bridge is the source's `state_index+1 < 2^30` assertion.
-/

namespace N83StateId30

def stateCount (k : Nat) : Nat :=
  83 * (k * (k - 1) / 2) + 42 * k

def encodeSlot (state selector fingerprint : Nat) : Nat :=
  fingerprint * 2 ^ 31 + (state + 1) * 2 + selector

def decodeState (word : Nat) : Nat :=
  ((word % 2 ^ 31) / 2) - 1

def decodeSelector (word : Nat) : Nat :=
  word % 2

def decodeFingerprint (word : Nat) : Nat :=
  word / 2 ^ 31

theorem slot_roundtrip (state selector fingerprint : Nat)
    (hstate : state + 1 < 2 ^ 30)
    (hselector : selector < 2) :
    decodeState (encodeSlot state selector fingerprint) = state ∧
      decodeSelector (encodeSlot state selector fingerprint) = selector ∧
      decodeFingerprint (encodeSlot state selector fingerprint) = fingerprint := by
  unfold encodeSlot decodeState decodeSelector decodeFingerprint
  omega

theorem slot_fits_five_bytes (state selector fingerprint : Nat)
    (hstate : state + 1 < 2 ^ 30)
    (hselector : selector < 2)
    (hfp : fingerprint < 2 ^ 9) :
    encodeSlot state selector fingerprint < 2 ^ 40 := by
  unfold encodeSlot
  omega

theorem k5086_fits : stateCount 5086 < 2 ^ 30 := by decide
theorem k5087_exceeds : 2 ^ 30 ≤ stateCount 5087 := by decide

theorem excluded_above_5086 (k : Nat) (hk : 5087 ≤ k) :
    2 ^ 30 ≤ stateCount k := by
  have hsub : 5086 ≤ k - 1 := by omega
  have hmul : 5087 * 5086 ≤ k * (k - 1) := Nat.mul_le_mul hk hsub
  unfold stateCount
  omega

#print axioms k5086_fits
#print axioms k5087_exceeds
#print axioms excluded_above_5086
#print axioms slot_roundtrip
#print axioms slot_fits_five_bytes

end N83StateId30
