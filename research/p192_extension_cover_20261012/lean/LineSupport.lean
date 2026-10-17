import AdaptiveBase
import Mathlib.Tactic.NormNum

/-!
Two fresh residual points on each split line section cannot by themselves
constrain old factor-base logarithms. The P-192 incidence and subgroup-point
distinctness instances are checked by PARI, Python, and Sage separately.
-/

namespace P192ExtensionCover

variable {K : Type*} [AddCommGroup K]

/-- Fixed-origin point representatives preserve the source divisor class.
    The concrete plane-quartic non-speciality test is checked separately. -/
theorem fixed_origin_point_divisor_identity (P A B C : K) :
    (P + B + C) - (A + B + C) = P - A := by
  abel

/-- One private point column per line is enough to extend any assignment of
    the old and companion point columns to the line equations. -/
theorem line_private_pairs_extend {n m : ℕ}
    (row : Fin m → (Fin n → K) → K) (rhs : Fin m → K) :
    ∀ old : Fin n → K, ∀ companion : Fin m → K,
      ∃! fresh : Fin m → K,
        ∀ j, row j old + fresh j + companion j = rhs j := by
  intro old companion
  refine ⟨fun j => rhs j - row j old - companion j, ?_, ?_⟩
  · intro j
    dsimp
    abel
  · intro other hother
    funext j
    calc
      other j = (row j old + other j + companion j) -
          row j old - companion j := by abel
      _ = rhs j - row j old - companion j := by rw [hother j]

/-- Exact unordered three-sum support cap for 9,925 distinct projected
    rational point classes, counting repetition. -/
theorem p192_line_point_three_sum_capacity :
    9925 * 9926 * 9927 / 6 = 162993977475 := by
  norm_num

/-- Relative to the certified accepted scalar-set lower bound, even the
    9,925-point three-sum support has ideal-uniform hit probability below
    two to the minus 152. -/
theorem p192_line_point_hit_bound_certificate :
    162993977475 * (2 ^ 152 : ℕ) <
      1569275433846670190958947353987675545133315652035359680560 := by
  norm_num

end P192ExtensionCover

#print axioms P192ExtensionCover.line_private_pairs_extend
#print axioms P192ExtensionCover.fixed_origin_point_divisor_identity
#print axioms P192ExtensionCover.p192_line_point_three_sum_capacity
#print axioms P192ExtensionCover.p192_line_point_hit_bound_certificate
