import Lean

/-!
Exact dynamic root-table bucket boundary for the N83 symmetric one-root
factor-base family. The source bridge is the pinned producer's retained-state
formula and its load-75 request, followed by a power-of-two table allocation
at least twice the request. This is a capacity calculation, not a timing claim.
-/

namespace N83OnlineFinegrid

def stateCount (k : Nat) : Nat :=
  83 * (k * (k - 1) / 2) + 42 * k

def dynamicRequest (k : Nat) : Nat :=
  (4 * stateCount k + 2) / 3

theorem k3114_states : stateCount 3114 = 402426891 := by decide
theorem k3115_states : stateCount 3115 = 402685395 := by decide
theorem k3114_request : dynamicRequest 3114 = 536569188 := by decide
theorem k3115_request : dynamicRequest 3115 = 536913860 := by decide

theorem k3114_uses_two_to_30_slots :
    2 ^ 29 < 2 * dynamicRequest 3114 ∧
      2 * dynamicRequest 3114 ≤ 2 ^ 30 := by decide

theorem k3115_uses_two_to_31_slots :
    2 ^ 30 < 2 * dynamicRequest 3115 ∧
      2 * dynamicRequest 3115 ≤ 2 ^ 31 := by decide

theorem k3114_state_id_fits : stateCount 3114 < 2 ^ 29 := by decide

theorem every_k_above_3114_crosses_bucket (k : Nat) (hk : 3115 ≤ k) :
    2 ^ 30 < 2 * dynamicRequest k := by
  have hsub : 3114 ≤ k - 1 := by omega
  have hmul : 3115 * 3114 ≤ k * (k - 1) := Nat.mul_le_mul hk hsub
  have hs : stateCount 3115 ≤ stateCount k := by
    unfold stateCount
    omega
  have hr : dynamicRequest 3115 ≤ dynamicRequest k := by
    unfold dynamicRequest
    omega
  have hb := k3115_uses_two_to_31_slots
  omega

#print axioms k3114_states
#print axioms k3115_states
#print axioms k3114_request
#print axioms k3115_request
#print axioms k3114_uses_two_to_30_slots
#print axioms k3115_uses_two_to_31_slots
#print axioms k3114_state_id_fits
#print axioms every_k_above_3114_crosses_bucket

end N83OnlineFinegrid
