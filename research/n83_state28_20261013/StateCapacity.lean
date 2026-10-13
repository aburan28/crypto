/- Numeric and encoding bounds for 28-byte N83 root states.
   The Rust field-root-to-byte-array bridge is checked by the source test and
   exact target witness replay; these theorems state the arithmetic bounds. -/

namespace N83State28

theorem root83_fits_eleven_bytes (x : Nat) (h : x < 2 ^ 83) :
    x % (256 ^ 11) = x := by
  exact Nat.mod_eq_of_lt (Nat.lt_trans h (by decide))

theorem k2202_packed_state_bytes_bounded (states : Nat)
    (h : states ≤ 83 * 2202 ^ 2) :
    28 * states ≤ 28 * (83 * 2202 ^ 2) :=
  Nat.mul_le_mul_left 28 h

theorem k2202_maximum_saves_more_than_seven_gib :
    7 * 2 ^ 30 < (48 - 28) * (83 * 2202 ^ 2) := by decide

#print axioms root83_fits_eleven_bytes
#print axioms k2202_packed_state_bytes_bounded
#print axioms k2202_maximum_saves_more_than_seven_gib

end N83State28
