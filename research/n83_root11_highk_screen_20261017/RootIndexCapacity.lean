import Init.Tactics
import Lean.Elab.Tactic.Omega
import Lean

/-!
The source-side one-root state-count guard simplifies to
`83 * k * k + k < 2 ^ idBoundBits` for the N83 symmetric state
layout. These theorems give the exact integer ceiling for the two
index-ID modes used by the experiment. They do not formalize the Rust
implementation or its memory and timing behavior.
-/

namespace N83RootIndexCapacity

theorem id29_exact_ceiling (k : Nat) :
    (83 * k * k + k < 2 ^ 30) ↔ k ≤ 3596 := by
  constructor
  · intro h
    by_cases hk : k ≤ 3596
    · exact hk
    have hge : 3597 ≤ k := by omega
    have hs : 3597 * 3597 ≤ k * k := Nat.mul_self_le_mul_self hge
    have hscale : 83 * (3597 * 3597) ≤ 83 * (k * k) :=
      Nat.mul_le_mul_left 83 hs
    simp only [Nat.mul_assoc] at h
    omega
  · intro hk
    have hs : k * k ≤ 3596 * 3596 := Nat.mul_self_le_mul_self hk
    have hscale : 83 * (k * k) ≤ 83 * (3596 * 3596) :=
      Nat.mul_le_mul_left 83 hs
    simp only [Nat.mul_assoc]
    omega

theorem id30_exact_ceiling (k : Nat) :
    (83 * k * k + k < 2 ^ 31) ↔ k ≤ 5086 := by
  constructor
  · intro h
    by_cases hk : k ≤ 5086
    · exact hk
    have hge : 5087 ≤ k := by omega
    have hs : 5087 * 5087 ≤ k * k := Nat.mul_self_le_mul_self hge
    have hscale : 83 * (5087 * 5087) ≤ 83 * (k * k) :=
      Nat.mul_le_mul_left 83 hs
    simp only [Nat.mul_assoc] at h
    omega
  · intro hk
    have hs : k * k ≤ 5086 * 5086 := Nat.mul_self_le_mul_self hk
    have hscale : 83 * (k * k) ≤ 83 * (5086 * 5086) :=
      Nat.mul_le_mul_left 83 hs
    simp only [Nat.mul_assoc]
    omega

#print axioms id29_exact_ceiling
#print axioms id30_exact_ceiling

end N83RootIndexCapacity
