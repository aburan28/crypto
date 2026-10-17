import LineSupport

/-!
Two plane-quartic line sections through the same point have equal residual
three-point classes after cancelling that common point. The concrete quartic
line cuts and quotient projections are checked by the coordinate receipts.
-/

namespace P192ExtensionCover

variable {K : Type*} [AddCommGroup K]

/-- Cancelling the common point of two line-section divisor identities
    produces an equal three-sum target. -/
theorem pencil_three_sum_eq_of_line_sections
    (P A B X Y U V H : K)
    (h₁ : P + A + X + Y = H)
    (h₂ : P + B + U + V = H) :
    A + X + Y = B + U + V := by
  have left : A + X + Y = H - P := by
    calc
      A + X + Y = (P + A + X + Y) - P := by abel
      _ = H - P := by rw [h₁]
  have right : B + U + V = H - P := by
    calc
      B + U + V = (P + B + U + V) - P := by abel
      _ = H - P := by rw [h₂]
  exact left.trans right.symm

end P192ExtensionCover

#print axioms P192ExtensionCover.pencil_three_sum_eq_of_line_sections
#print axioms P192ExtensionCover.line_private_pairs_extend
