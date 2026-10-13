import Init

namespace Nonisomorphic617

set_option maxRecDepth 8192
set_option maxHeartbeats 0

/-- The destination/source coefficient ratio 590 is not a sixth power
of any field element represented by a canonical residue modulo 617. -/
theorem ratioNotSixthPower : ∀ u : Fin 617, u.val ^ 6 % 617 ≠ 590 := by
  decide

/-- The target prime is coprime to the rational degree-three kernel. -/
theorem targetCoprimeKernel : Nat.Coprime 103 3 := by decide

/-- On the order-103 scalar ring, a point killed by the kernel order 3
has zero scalar. Applying this to an isogeny kernel still requires the
group-law bridge from the concrete curve. -/
theorem kernelScalarVanishes :
    ∀ k : Fin 103, (3 * k : Fin 103) = 0 → k = 0 := by
  decide

#print axioms Nonisomorphic617.ratioNotSixthPower
#print axioms Nonisomorphic617.targetCoprimeKernel
#print axioms Nonisomorphic617.kernelScalarVanishes

end Nonisomorphic617
