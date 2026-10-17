import Mathlib.GroupTheory.OrderOfElement

/-!
The finite-group order calculation behind the minimal pairing output field.
Instantiate `x` with the class of `p` in `(ZMod r)ˣ` when `gcd p r = 1`.
Then `orderOf x` is the minimal degree `d` containing the r-th roots of
unity, while `orderOf (x ^ m)` is the embedding degree over `q = p ^ m`.
The finite-field roots-of-unity identification is a separate bridge.
-/

namespace WeakCurveLeaderboard

theorem embedding_degree_formula {G : Type*} [Group G] (x : G) (m : ℕ)
    (hm : m ≠ 0) :
    orderOf (x ^ m) = orderOf x / Nat.gcd (orderOf x) m := by
  exact orderOf_pow' x hm

theorem minimal_degree_divides_pairing_extension {G : Type*} [Group G]
    (x : G) (m : ℕ) :
    orderOf x ∣ m * orderOf (x ^ m) := by
  apply orderOf_dvd_of_pow_eq_one
  rw [pow_mul, pow_orderOf_eq_one]

theorem minimal_degree_is_minimal {G : Type*} [Group G]
    (x : G) (n : ℕ) :
    x ^ n = 1 ↔ orderOf x ∣ n := by
  exact (orderOf_dvd_iff_pow_eq_one).symm

end WeakCurveLeaderboard
