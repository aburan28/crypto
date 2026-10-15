import RelationSupport
import Mathlib.Algebra.Order.BigOperators.Group.Finset

/-!
An integer weight can count how many sampler seeds reach each target class.
No uniformity assumption is needed: the maximum target weight controls all
finite support-hit mass, including geometry-guided or conditioned query laws.
-/

namespace WeakCurveLeaderboard

theorem weighted_relation_support_bound {A : Type*} [DecidableEq A]
    (B m : ℕ) (target : Sym (Fin B) m → A)
    (weight : A → ℕ) (maxWeight : ℕ)
    (hmax : ∀ a, weight a ≤ maxWeight) :
    ∑ a ∈ ((Finset.univ : Finset (Sym (Fin B) m)).image target), weight a ≤
      (B + m - 1).choose m * maxWeight := by
  classical
  let support := (Finset.univ : Finset (Sym (Fin B) m)).image target
  calc
    ∑ a ∈ support, weight a ≤ ∑ _a ∈ support, maxWeight := by
      exact Finset.sum_le_sum (fun a _ => hmax a)
    _ = support.card * maxWeight := by simp
    _ ≤ (B + m - 1).choose m * maxWeight := by
      exact Nat.mul_le_mul_right maxWeight
        (unordered_relation_support_bound B m target)

theorem conditioned_weighted_relation_support_bound {A : Type*} [DecidableEq A]
    (B m : ℕ) (target : Sym (Fin B) m → A)
    (accepted : Finset A) (weight : A → ℕ) (maxWeight : ℕ)
    (hmax : ∀ a ∈ accepted, weight a ≤ maxWeight) :
    ∑ a ∈ accepted ∩
      ((Finset.univ : Finset (Sym (Fin B) m)).image target), weight a ≤
      (B + m - 1).choose m * maxWeight := by
  classical
  let support := (Finset.univ : Finset (Sym (Fin B) m)).image target
  calc
    ∑ a ∈ accepted ∩ support, weight a ≤
        ∑ _a ∈ accepted ∩ support, maxWeight := by
      exact Finset.sum_le_sum (fun a ha => hmax a (Finset.mem_inter.mp ha).1)
    _ = (accepted ∩ support).card * maxWeight := by simp
    _ ≤ (B + m - 1).choose m * maxWeight := by
      exact Nat.mul_le_mul_right maxWeight
        (conditioned_relation_support_bound B m target accepted)

end WeakCurveLeaderboard

#print axioms WeakCurveLeaderboard.weighted_relation_support_bound
#print axioms WeakCurveLeaderboard.conditioned_weighted_relation_support_bound
