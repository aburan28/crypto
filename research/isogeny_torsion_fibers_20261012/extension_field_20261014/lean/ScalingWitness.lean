import Init

/-!
Exact coefficient witnesses for the four selected quadratic-extension
isomorphisms. These theorems check the public scaling coefficients only;
they do not formalize elliptic-curve points or the group law.
-/

theorem scaling_sixth_power_factorization (u : Int) :
    (u ^ 2 + 3) * (u ^ 4 - 3 * u ^ 2 + 9) = u ^ 6 + 27 := by
  simp [Int.pow_succ, Int.pow_zero, Int.mul_add, Int.add_mul,
    Int.sub_eq_add_neg, Int.mul_neg, Int.mul_assoc,
    Int.mul_left_comm, Int.mul_comm]
  omega

theorem scaling_witness_p61 :
    27 ^ 2 % 61 = 58 ∧ 27 ^ 6 % 61 = 34 := by decide

theorem scaling_witness_p271 :
    57 ^ 2 % 271 = 268 ∧ 57 ^ 6 % 271 = 244 := by decide

theorem scaling_witness_p787 :
    28 ^ 2 % 787 = 784 ∧ 28 ^ 6 % 787 = 760 := by decide

theorem scaling_witness_p3319 :
    264 ^ 2 % 3319 = 3316 ∧ 264 ^ 6 % 3319 = 3292 := by decide

#print axioms scaling_witness_p61
#print axioms scaling_witness_p271
#print axioms scaling_witness_p787
#print axioms scaling_witness_p3319
#print axioms scaling_sixth_power_factorization
