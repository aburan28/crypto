import Mathlib.Tactic.Abel
import Mathlib.Tactic.NormNum

/-!
The six-point inverse fiber divisor is Frobenius-fixed when the third
Frobenius iterate either fixes the two points or swaps them. This does not
assume that either individual point is rational over the cubic field.
The concrete cover-point orbits and cubic-section uniqueness are replayed
by Sage and PARI, respectively.
-/

namespace P192ExtensionCover

variable {G : Type*} [AddCommGroup G]

def sixFiber (frobenius : G →+ G) (plus minus : G) : G :=
  plus + minus + frobenius plus + frobenius minus +
    frobenius (frobenius plus) + frobenius (frobenius minus)

theorem six_fiber_fixed_of_three_sum (frobenius : G →+ G)
    (plus minus : G)
    (hthree : frobenius (frobenius (frobenius plus)) +
      frobenius (frobenius (frobenius minus)) = plus + minus) :
    frobenius (sixFiber frobenius plus minus) =
      sixFiber frobenius plus minus := by
  calc
    frobenius (sixFiber frobenius plus minus) =
        frobenius plus + frobenius minus +
          frobenius (frobenius plus) + frobenius (frobenius minus) +
          (frobenius (frobenius (frobenius plus)) +
            frobenius (frobenius (frobenius minus))) := by
            simp only [sixFiber, map_add]
            abel
    _ = frobenius plus + frobenius minus +
          frobenius (frobenius plus) + frobenius (frobenius minus) +
          (plus + minus) := by rw [hthree]
    _ = sixFiber frobenius plus minus := by unfold sixFiber; abel

theorem six_fiber_fixed_of_three_swap (frobenius : G →+ G)
    (plus minus : G)
    (hplus : frobenius (frobenius (frobenius plus)) = minus)
    (hminus : frobenius (frobenius (frobenius minus)) = plus) :
    frobenius (sixFiber frobenius plus minus) =
      sixFiber frobenius plus minus := by
  apply six_fiber_fixed_of_three_sum frobenius plus minus
  rw [hplus, hminus]
  abel

theorem six_fiber_fixed_of_three_fixed (frobenius : G →+ G)
    (plus minus : G)
    (hplus : frobenius (frobenius (frobenius plus)) = plus)
    (hminus : frobenius (frobenius (frobenius minus)) = minus) :
    frobenius (sixFiber frobenius plus minus) =
      sixFiber frobenius plus minus := by
  apply six_fiber_fixed_of_three_sum frobenius plus minus
  rw [hplus, hminus]

theorem p192_line_quotient_character_cells :
    30 + 26 + 40 + 45 = 141 := by norm_num

end P192ExtensionCover

#print axioms P192ExtensionCover.six_fiber_fixed_of_three_sum
#print axioms P192ExtensionCover.six_fiber_fixed_of_three_swap
#print axioms P192ExtensionCover.six_fiber_fixed_of_three_fixed
#print axioms P192ExtensionCover.p192_line_quotient_character_cells
