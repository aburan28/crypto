import Lean

/-!
Exact N83 root-state and table-slot capacity arithmetic for the bounded
online range. The implementation bridge is the pinned source's state-count
and dynamic/fixed slot request policy; the extra GiB is a frozen working-
memory reserve, not an observed RSS value.
-/

namespace N83OnlineRange

def stateCount (k : Nat) : Nat :=
  83 * (k * (k - 1) / 2) + 42 * k

def dynamicRequest (k : Nat) : Nat :=
  (4 * stateCount k + 2) / 3

def rssCapBytes : Nat := 61 * 2 ^ 30 / 4
def workingReserveBytes : Nat := 2 ^ 30
def oneRootStateBytes (k : Nat) : Nat := 11 * stateCount k
def tableBytes (slots : Nat) : Nat := 5 * slots

theorem k3450_state_count : stateCount 3450 = 493955475 := by decide
theorem k3584_state_count : stateCount 3584 = 533071616 := by decide
theorem k4096_state_count : stateCount 4096 = 696256512 := by decide

theorem k3450_dynamic_bucket :
    2 ^ 30 < 2 * dynamicRequest 3450 ∧
      2 * dynamicRequest 3450 ≤ 2 ^ 31 := by decide

theorem k4096_dynamic_bucket :
    2 ^ 30 < 2 * dynamicRequest 4096 ∧
      2 * dynamicRequest 4096 ≤ 2 ^ 31 := by decide

theorem k3450_dynamic_exceeds_guard :
    tableBytes (2 ^ 31) + oneRootStateBytes 3450 +
      workingReserveBytes > rssCapBytes := by decide

theorem k4096_dynamic_exceeds_guard :
    tableBytes (2 ^ 31) + oneRootStateBytes 4096 +
      workingReserveBytes > rssCapBytes := by decide

theorem k3450_fixed_fits_capacity_model :
    tableBytes (2 ^ 30) + oneRootStateBytes 3450 +
      workingReserveBytes ≤ rssCapBytes := by decide

theorem k3450_fixed_root_count_fits :
    2 * stateCount 3450 < 2 ^ 30 := by decide

theorem k3584_fixed_fits_capacity_model :
    tableBytes (2 ^ 30) + oneRootStateBytes 3584 +
      workingReserveBytes ≤ rssCapBytes := by decide

theorem k3584_fixed_root_count_fits :
    2 * stateCount 3584 < 2 ^ 30 := by decide

theorem k3596_state_id_fits : stateCount 3596 < 2 ^ 29 := by decide
theorem k3597_state_id_exceeds : 2 ^ 29 ≤ stateCount 3597 := by decide

theorem packed_state_id_excluded_above_3596 (k : Nat) (hk : 3597 ≤ k) :
    2 ^ 29 ≤ stateCount k := by
  have hsub : 3596 ≤ k - 1 := by omega
  have hmul : 3597 * 3596 ≤ k * (k - 1) := Nat.mul_le_mul hk hsub
  unfold stateCount
  omega

#print axioms k3450_state_count
#print axioms k3584_state_count
#print axioms k4096_state_count
#print axioms k3450_dynamic_bucket
#print axioms k4096_dynamic_bucket
#print axioms k3450_dynamic_exceeds_guard
#print axioms k4096_dynamic_exceeds_guard
#print axioms k3450_fixed_fits_capacity_model
#print axioms k3450_fixed_root_count_fits
#print axioms k3584_fixed_fits_capacity_model
#print axioms k3584_fixed_root_count_fits
#print axioms k3596_state_id_fits
#print axioms k3597_state_id_exceeds
#print axioms packed_state_id_excluded_above_3596

end N83OnlineRange
