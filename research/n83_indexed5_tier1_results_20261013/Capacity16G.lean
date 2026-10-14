/-!
The regular five-byte N83 table has a 5 GiB bucket through K1798 and
crosses to a 10 GiB bucket at K1799, conditional on the source-loop bound
of `83 * K^2` states and two candidate keys per state.
-/

namespace N83Indexed5Tier1

def stateBound (k : Nat) : Nat := 83 * (k * k)

theorem state_bound_monotone {a b : Nat} (h : a ≤ b) :
    stateBound a ≤ stateBound b := by
  unfold stateBound
  have hsq : a * a ≤ b * b := Nat.mul_le_mul h h
  exact Nat.mul_le_mul_left 83 hsq

theorem k1798_fits_five_gibibyte_bucket :
    4 * stateBound 1798 ≤ 2 ^ 30 := by decide

theorem all_smaller_k_fit_bucket {k : Nat} (h : k ≤ 1798) :
    4 * stateBound k ≤ 2 ^ 30 := by
  have hbound := state_bound_monotone h
  exact Nat.le_trans (Nat.mul_le_mul_left 4 hbound)
    k1798_fits_five_gibibyte_bucket

theorem k1799_crosses_five_gibibyte_bucket :
    2 ^ 30 < 4 * stateBound 1799 := by decide

theorem all_larger_k_cross_bucket {k : Nat} (h : 1799 ≤ k) :
    2 ^ 30 < 4 * stateBound k := by
  have hbound := state_bound_monotone h
  have hscaled := Nat.mul_le_mul_left 4 hbound
  exact Nat.lt_of_lt_of_le k1799_crosses_five_gibibyte_bucket hscaled

theorem five_gibibytes_of_slots :
    5 * (2 ^ 30) = 5368709120 := by decide

theorem ten_gibibytes_of_slots :
    5 * (2 ^ 31) = 10 * (2 ^ 30) := by decide

#print axioms state_bound_monotone
#print axioms k1798_fits_five_gibibyte_bucket
#print axioms all_smaller_k_fit_bucket
#print axioms k1799_crosses_five_gibibyte_bucket
#print axioms all_larger_k_cross_bucket
#print axioms five_gibibytes_of_slots
#print axioms ten_gibibytes_of_slots

end N83Indexed5Tier1
