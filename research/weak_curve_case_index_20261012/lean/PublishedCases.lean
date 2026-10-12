import Mathlib.Data.ZMod.Basic
import Mathlib.Algebra.Polynomial.Basic
import Mathlib.Tactic.Ring

set_option autoImplicit false

namespace PublishedCases

private abbrev P := Polynomial (ZMod 2)

/-- The factorization used by the cubic binary ramification-mask screen.
This is the polynomial step of IMC 2015, Section 3.3, not a covering theorem. -/
theorem cubic_frobenius_factorization :
    (Polynomial.X ^ 2 + Polynomial.X + 1 : P) * (Polynomial.X + 1) =
      Polynomial.X ^ 3 + 1 := by
  have htwo : (2 : P) = 0 := by
    change Polynomial.C (2 : ZMod 2) = 0
    have hz : (2 : ZMod 2) = 0 := by decide
    simp [hz]
  ring_nf
  simp [htwo]

/-- The factorization underlying the degree-seven representation choices.
The factors have degrees one, three, and three over the binary field. -/
theorem septic_frobenius_factorization :
    (Polynomial.X + 1 : P) *
      (Polynomial.X ^ 3 + Polynomial.X + 1) *
      (Polynomial.X ^ 3 + Polynomial.X ^ 2 + 1) =
      Polynomial.X ^ 7 + 1 := by
  have htwo : (2 : P) = 0 := by
    change Polynomial.C (2 : ZMod 2) = 0
    have hz : (2 : ZMod 2) = 0 := by decide
    simp [hz]
  have hfour : (4 : P) = 0 := by
    calc
      (4 : P) = 2 + 2 := by ring
      _ = 0 := by rw [htwo]; simp
  ring_nf
  simp [htwo, hfour]

end PublishedCases
