import Std.Tactic

/- The symmetric N83 root-only construction allocates one state for each
   off-diagonal orbit pair and 42 diagonal states per factor-base column.
   The four-byte root slot encodes a state index in 29 bits, and the Rust
   producer requires the state count to be strictly below 2^29. -/
def stateBound (k : Nat) : Nat := 83 * k * (k - 1) / 2 + 42 * k

theorem stateBound_3596 : stateBound 3596 = 536647262 := by decide
theorem stateBound_3597 : stateBound 3597 = 536945772 := by decide
theorem fits_3596 : stateBound 3596 < 2^29 := by decide
theorem excludes_3597 : 2^29 ≤ stateBound 3597 := by decide

#print axioms stateBound_3596
#print axioms stateBound_3597
#print axioms fits_3596
#print axioms excludes_3597

theorem excludes_all_larger (k : Nat) (h : 3597 ≤ k) : 2^29 ≤ stateBound k := by
  have hprev : 3596 ≤ k - 1 := by omega
  have hp : 3597 * 3596 ≤ k * (k - 1) := Nat.mul_le_mul h hprev
  have hbig : 83 * (3597 * 3596) ≤ 83 * (k * (k - 1)) := Nat.mul_le_mul_left 83 hp
  have hdiv : 83 * (3597 * 3596) / 2 ≤ 83 * (k * (k - 1)) / 2 := Nat.div_le_div_right hbig
  have hdiag : 42 * 3597 ≤ 42 * k := Nat.mul_le_mul_left 42 h
  have hmono : stateBound 3597 ≤ stateBound k := by
    simpa only [stateBound, Nat.mul_assoc] using Nat.add_le_add hdiv hdiag
  exact Nat.le_trans excludes_3597 hmono

#print axioms excludes_all_larger

theorem max_k_current_layout :
    stateBound 3596 < 2^29 ∧
    ∀ k : Nat, 3596 < k → 2^29 ≤ stateBound k := by
  constructor
  · exact fits_3596
  · intro k hk
    exact excludes_all_larger k (by omega)

#print axioms max_k_current_layout
