import Std.Tactic

set_option autoImplicit false

namespace CensusCounting

structure Flags where
  cubic : Bool
  quadratic : Bool

def countCubic : List Flags → Nat
  | [] => 0
  | x :: xs => (if x.cubic then 1 else 0) + countCubic xs
def countQuadratic : List Flags → Nat
  | [] => 0
  | x :: xs => (if x.quadratic then 1 else 0) + countQuadratic xs
def countBoth : List Flags → Nat
  | [] => 0
  | x :: xs => (if x.cubic && x.quadratic then 1 else 0) + countBoth xs
def countUnion : List Flags → Nat
  | [] => 0
  | x :: xs => (if x.cubic || x.quadratic then 1 else 0) + countUnion xs
def countZero : List Flags → Nat
  | [] => 0
  | x :: xs => (if x.cubic || x.quadratic then 0 else 1) + countZero xs

/-- Exact inclusion-exclusion for arbitrary finite lists of two family flags.
This theorem verifies aggregation, not the truth of geometric input labels. -/
theorem inclusion_exclusion (xs : List Flags) :
    countCubic xs + countQuadratic xs = countUnion xs + countBoth xs := by
  induction xs with
  | nil => rfl
  | cons x xs ih =>
    cases hc : x.cubic <;> cases hq : x.quadratic <;>
      simp [countCubic, countQuadratic, countUnion, countBoth, hc, hq] <;> omega

theorem union_zero_partition (xs : List Flags) :
    countUnion xs + countZero xs = xs.length := by
  induction xs with
  | nil => rfl
  | cons x xs ih =>
    cases hc : x.cubic <;> cases hq : x.quadratic <;>
      simp [countUnion, countZero, hc, hq] <;> omega

theorem union_bounded_by_population (xs : List Flags) :
    countUnion xs ≤ xs.length := by
  have := union_zero_partition xs
  omega

end CensusCounting
