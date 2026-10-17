import Mathlib.Data.Sym.Card
import Mathlib.Tactic.NormNum

/-!
Any fixed-base relation policy that chooses an unordered multiset of `m`
base labels has at most `multichoose B m` distinct outputs, regardless of
its curve, group law, or query predicate. The concrete P-192 base-size
certificate is supplied by the finite-field geometry receipt.
-/

namespace P192ExtensionCover

theorem unordered_policy_support_card_le_multichoose
    {Base Output : Type*} [Fintype Base] [DecidableEq Base]
    [DecidableEq Output] (m : ℕ)
    (target : Sym Base m → Output) :
    ((Finset.univ : Finset (Sym Base m)).image target).card ≤
      Nat.multichoose (Fintype.card Base) m := by
  classical
  calc
    ((Finset.univ : Finset (Sym Base m)).image target).card ≤
        (Finset.univ : Finset (Sym Base m)).card := Finset.card_image_le
    _ = Nat.multichoose (Fintype.card Base) m := by
      simp [Sym.card_sym_eq_multichoose]

/-- Arithmetic certificate for the optimistic 1,411,637-label P-192 base. -/
theorem all_pivot_three_sum_capacity :
    Nat.multichoose 1411637 3 = 468833646178255739 := by
  rw [Nat.multichoose_eq, Nat.choose_eq_descFactorial_div_factorial]
  norm_num [Nat.descFactorial, Nat.factorial]

/-- Under the certified accepted-scalar lower bound, even that generous
    three-sum capacity has ideal-uniform hit probability below `2^-131`. -/
theorem all_pivot_three_sum_hit_bound :
    468833646178255739 * (2 ^ 131 : ℕ) <
      1569275433846670190958947353987675545133315652035359680560 := by
  norm_num

/-- The same optimistic base first has enough unordered tuple labels to
    match the accepted-set lower bound at arity eleven. This is a cardinality
    threshold, independent of relation yield. -/
theorem all_pivot_arity_capacity_threshold :
    Nat.multichoose 1411637 10 <
      1569275433846670190958947353987675545133315652035359680560 ∧
    1569275433846670190958947353987675545133315652035359680560 ≤
      Nat.multichoose 1411637 11 := by
  rw [Nat.multichoose_eq, Nat.choose_eq_descFactorial_div_factorial,
      Nat.multichoose_eq, Nat.choose_eq_descFactorial_div_factorial]
  norm_num [Nat.descFactorial, Nat.factorial]

end P192ExtensionCover

#print axioms P192ExtensionCover.unordered_policy_support_card_le_multichoose
#print axioms P192ExtensionCover.all_pivot_three_sum_capacity
#print axioms P192ExtensionCover.all_pivot_three_sum_hit_bound
#print axioms P192ExtensionCover.all_pivot_arity_capacity_threshold
