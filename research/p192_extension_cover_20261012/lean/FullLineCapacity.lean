import Mathlib.Tactic.GCongr
import Mathlib.Tactic.NormNum

/-!
Even if every rational outside residual from the complete P-192 old-anchored
line census supplied one distinct usable subgroup class, the three-sum
support remains bounded by a small fraction of the certified four-cell
known-scalar population.  The actual usable-base and character-population
hypotheses stay explicit.
-/

namespace P192ExtensionCover

/-- Existing 9,925 selected classes plus at most 1,360,948 rational outside
    point presentations from the complete old-anchored screen. -/
theorem p192_old_anchor_maximal_base_size :
    9925 + 1360948 = 1370873 := by
  norm_num

/-- Exact multiset three-sum capacity at the maximal hypothetical base size. -/
theorem p192_old_anchor_maximal_three_sum_capacity :
    1370873 * 1370874 * 1370875 / 6 = 429379561999436125 := by
  norm_num

/-- Conditional ideal-uniform hit bound for any actual usable base drawn
    from the selected and old-anchored rational residual points. `M` is its
    distinct three-sum support and `D` its accepted-scalar cell population. -/
theorem p192_old_anchor_three_sum_hit_bound
    (B M D : ℕ)
    (hB : B ≤ 1370873)
    (hM : M ≤ B * (B + 1) * (B + 2) / 6)
    (hD : 1569275433846670190958947353987675545133315652035359680560 ≤ D) :
    M * 2 ^ 131 < D := by
  have hprod : B * (B + 1) * (B + 2) ≤
      1370873 * 1370874 * 1370875 := by
    gcongr <;> omega
  have hcap : M ≤ 429379561999436125 := by
    calc
      M ≤ B * (B + 1) * (B + 2) / 6 := hM
      _ ≤ (1370873 * 1370874 * 1370875) / 6 := Nat.div_le_div_right hprod
      _ = 429379561999436125 := p192_old_anchor_maximal_three_sum_capacity
  calc
    M * 2 ^ 131 ≤ 429379561999436125 * 2 ^ 131 :=
      Nat.mul_le_mul_right _ hcap
    _ < 1569275433846670190958947353987675545133315652035359680560 := by
      norm_num
    _ ≤ D := hD

end P192ExtensionCover

#print axioms P192ExtensionCover.p192_old_anchor_maximal_base_size
#print axioms P192ExtensionCover.p192_old_anchor_maximal_three_sum_capacity
#print axioms P192ExtensionCover.p192_old_anchor_three_sum_hit_bound
