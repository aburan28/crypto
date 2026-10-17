import Mathlib.Algebra.Polynomial.Eval.Defs
import Mathlib.Tactic.NormNum

/-!
An exact Bézout identity for two representative x-polynomials excludes a
common affine support point. The P-192 coefficient instances and gcd results
are checked separately by Python and PARI.
-/

namespace P192ExtensionCover

open Polynomial

variable {F : Type*} [Field F]

/-- A polynomial Bézout identity forbids a common root over the field. -/
theorem no_common_x_of_bezout (f g a b : F[X])
    (hbezout : a * f + b * g = 1) (x : F) :
    ¬(f.eval x = 0 ∧ g.eval x = 0) := by
  rintro ⟨hf, hg⟩
  have heval := congrArg (fun h : F[X] => h.eval x) hbezout
  have hbad : (0 : F) = 1 := by
    simpa [hf, hg] using heval
  exact zero_ne_one hbad

/-- If two affine divisor supports are cut out by coprime x-polynomials,
    no affine point belongs to both supports. -/
theorem affine_support_disjoint_of_bezout (f g a b : F[X])
    (hbezout : a * f + b * g = 1)
    (left right : Set (F × F))
    (hleft : ∀ q ∈ left, f.eval q.1 = 0)
    (hright : ∀ q ∈ right, g.eval q.1 = 0) :
    ∀ q, q ∈ left → q ∈ right → False := by
  intro q hqleft hqright
  exact no_common_x_of_bezout f g a b hbezout q.1
    ⟨hleft q hqleft, hright q hqright⟩

/-- Exact integer boundary for the fully split 20-class support policy:
    at most 143 eligible rank-bearing targets among the certified accepted
    scalar set. -/
theorem p192_split_support_hit_bound_certificate :
    143 * (2 ^ 182 : ℕ) <
      1569275433846670190958947353987675545133315652035359680560 := by
  norm_num

end P192ExtensionCover

#print axioms P192ExtensionCover.no_common_x_of_bezout
#print axioms P192ExtensionCover.affine_support_disjoint_of_bezout
#print axioms P192ExtensionCover.p192_split_support_hit_bound_certificate
