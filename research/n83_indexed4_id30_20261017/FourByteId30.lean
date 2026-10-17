import Init.Tactics
import Lean.Elab.Tactic.Omega
import Lean

/-!
The four-byte root slot has a one-bit fingerprint, a nonzero 30-bit
state ordinal, and a one-bit root selector. Zero is the empty marker.
These arithmetic statements correspond to the explicit Rust encode/decode
functions; the table's group-law and runtime behavior are checked separately.
-/

namespace N83FourByteId30

def encode (state selector fingerprint : Nat) : Nat :=
  fingerprint * 2 ^ 31 + (state + 1) * 2 + selector

theorem encode_nonzero_and_fits
    (state selector fingerprint : Nat)
    (hstate : state + 1 < 2 ^ 30)
    (hselector : selector < 2)
    (hfingerprint : fingerprint < 2) :
    0 < encode state selector fingerprint ∧
      encode state selector fingerprint < 2 ^ 32 := by
  unfold encode
  omega

theorem decode_roundtrip
    (state selector fingerprint : Nat)
    (hstate : state + 1 < 2 ^ 30)
    (hselector : selector < 2) :
    ((encode state selector fingerprint % (2 ^ 31)) / 2 - 1 = state) ∧
      ((encode state selector fingerprint % (2 ^ 31)) % 2 = selector) ∧
      (encode state selector fingerprint / (2 ^ 31) = fingerprint) := by
  unfold encode
  omega

theorem load75_bound (states slots : Nat)
    (hslots : 2 * ((4 * states + 2) / 3) ≤ slots) :
    8 * states ≤ 3 * slots := by
  omega

#print axioms encode_nonzero_and_fits
#print axioms decode_roundtrip
#print axioms load75_bound

end N83FourByteId30
