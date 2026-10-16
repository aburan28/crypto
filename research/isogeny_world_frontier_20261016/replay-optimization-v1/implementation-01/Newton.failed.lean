import Lean.Elab.Tactic.Omega

namespace FixedModulusNewton

-- This statement is over all integer inputs, rather than evaluation points.
theorem newton_scalar (b g : Int) :
    (b * g) * (2 - b * g) = 1 - (1 - b * g) * (1 - b * g) := by
  generalize b * g = t
  simp only [Int.mul_sub, Int.sub_mul, Int.mul_one, Int.one_mul, Int.mul_two]
  omega

-- Coefficient n of a square is the finite convolution of coefficients 0..n.
def convolutionPrefix (e : Nat → Int) (n : Nat) : Nat → Int
  | 0 => 0
  | r + 1 => convolutionPrefix e n r + e r * e (n - r)

-- A residual with its first k coefficients zero has a square with its first
-- 2k coefficients zero. This proves the precision-doubling structural step
-- for integer coefficient sequences, independently of the sampled IDC1.
theorem residual_square_precision
    (e : Nat → Int) (k : Nat) (h : ∀ i, i < k → e i = 0)
    (n : Nat) (hn : n < 2 * k) :
    convolutionPrefix e n (n + 1) = 0 := by
  have terms : ∀ i, i ≤ n → e i * e (n - i) = 0 := by
    intro i hi
    by_cases hik : i < k
    · simp only [h i hik, Int.zero_mul]
    · have hj : n - i < k := by omega
      simp only [h (n - i) hj, Int.mul_zero]
  have sums : ∀ r, r ≤ n + 1 → convolutionPrefix e n r = 0 := by
    intro r
    induction r with
    | zero => intro _; rfl
    | succ r ih =>
      intro hr
      have hi : r ≤ n := by omega
      simp only [convolutionPrefix, ih (by omega), terms r hi, Int.add_zero]
  exact sums (n + 1) (by omega)

#print axioms newton_scalar
#print axioms residual_square_precision
end FixedModulusNewton

