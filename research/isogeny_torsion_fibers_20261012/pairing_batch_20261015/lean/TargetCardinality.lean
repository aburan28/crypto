import Init

/-!
The arithmetic implication of the distinct-target cardinality bound.
The finite-group counting lemma that supplies `n ≤ r - 1` is recorded in
PROOF_MAP.md as a missing bridge, not hidden as an axiom of this theorem.
-/

theorem order19_batch_bound (n : Nat) (h : n ≤ 19 - 1) : n < 64 := by
  omega

theorem order61_batch_bound (n : Nat) (h : n ≤ 61 - 1) : n < 64 := by
  omega

theorem order_for_64_nonidentity (r : Nat) (h : 64 ≤ r - 1) : 65 ≤ r := by
  omega

#print axioms order19_batch_bound
#print axioms order61_batch_bound
#print axioms order_for_64_nonidentity
