/-!
Exact symmetric-state counts and five-byte table buckets for the three
one-root dynamic candidates in the frozen N83 screen. The source bridge is
the checked traversal assertion that every symmetric state is retained.
-/

namespace N83Root11Dynamic

def stateCount (k : Nat) : Nat :=
  83 * (k * (k - 1) / 2) + 42 * k

def requestEntries (k : Nat) : Nat :=
  (4 * stateCount k + 2) / 3

theorem k1024_state_count : stateCount 1024 = 43516416 := by decide
theorem k1271_state_count : stateCount 1271 = 67041437 := by decide
theorem k1536_state_count : stateCount 1536 = 97911552 := by decide

theorem k1024_request : requestEntries 1024 = 58021888 := by decide
theorem k1271_request : requestEntries 1271 = 89388583 := by decide
theorem k1536_request : requestEntries 1536 = 130548736 := by decide

theorem k1024_slot_bucket :
    2 ^ 26 < 2 * requestEntries 1024 ∧
      2 * requestEntries 1024 ≤ 2 ^ 27 := by decide

theorem k1271_slot_bucket :
    2 ^ 27 < 2 * requestEntries 1271 ∧
      2 * requestEntries 1271 ≤ 2 ^ 28 := by decide

theorem k1536_slot_bucket :
    2 ^ 27 < 2 * requestEntries 1536 ∧
      2 * requestEntries 1536 ≤ 2 ^ 28 := by decide

theorem k1536_state_ids_fit : stateCount 1536 + 1 < 2 ^ 29 := by decide

#print axioms k1024_state_count
#print axioms k1271_state_count
#print axioms k1536_state_count
#print axioms k1024_request
#print axioms k1271_request
#print axioms k1536_request
#print axioms k1024_slot_bucket
#print axioms k1271_slot_bucket
#print axioms k1536_slot_bucket
#print axioms k1536_state_ids_fit

end N83Root11Dynamic
