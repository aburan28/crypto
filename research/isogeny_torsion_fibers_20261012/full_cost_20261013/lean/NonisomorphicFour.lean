import Init

namespace NonisomorphicFour

set_option maxRecDepth 65536
set_option maxHeartbeats 0

/-- The four concrete destination/source coefficient ratios are not sixth
powers among canonical residues. Linking this finite statement to the
classification of base-field isomorphisms of the two Weierstrass models
remains a mathematical bridge outside this Lean module. -/
theorem ratio617 : ∀ u : Fin 617, u.val ^ 6 % 617 ≠ 590 := by decide
theorem ratio2633 : ∀ u : Fin 2633, u.val ^ 6 % 2633 ≠ 2606 := by decide
-- The analogous 9677 and 38561 statements remain open in this Lean module.
-- Kernel evaluation of their direct `decide` proofs exceeded the bounded
-- build attempt recorded in build-v1.log and build-v2.log.

theorem coprime103 : Nat.Coprime 103 3 := by decide
theorem coprime439 : Nat.Coprime 439 3 := by decide
theorem coprime1613 : Nat.Coprime 1613 3 := by decide
theorem coprime6427 : Nat.Coprime 6427 3 := by decide

#print axioms NonisomorphicFour.ratio617
#print axioms NonisomorphicFour.ratio2633
#print axioms NonisomorphicFour.coprime103
#print axioms NonisomorphicFour.coprime439
#print axioms NonisomorphicFour.coprime1613
#print axioms NonisomorphicFour.coprime6427

end NonisomorphicFour
