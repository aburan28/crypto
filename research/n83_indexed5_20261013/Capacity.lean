/-!
Numeric capacity gate for the opt-in five-byte N83 root index. The Rust
source-loop bridge is kept separate: these statements assume its bound of
`83 * K^2` states and at most two root keys per state.
-/

namespace N83Indexed5

theorem k2543_state_id_fits : 83 * 2543 ^ 2 < 2 ^ 29 := by decide

theorem k2544_state_id_crosses : 2 ^ 29 ≤ 83 * 2544 ^ 2 := by decide

theorem k2543_regular_request_fits_bucket :
    4 * (83 * 2543 ^ 2) ≤ 2 ^ 31 := by decide

theorem k2544_regular_request_crosses_bucket :
    2 ^ 31 < 4 * (83 * 2544 ^ 2) := by decide

theorem ten_gibibytes_of_slots : 5 * (2 ^ 31) = 10 * (2 ^ 30) := by decide

#print axioms k2543_state_id_fits
#print axioms k2544_state_id_crosses
#print axioms k2543_regular_request_fits_bucket
#print axioms k2544_regular_request_crosses_bucket
#print axioms ten_gibibytes_of_slots

end N83Indexed5
