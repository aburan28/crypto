import Init

namespace FiberCofactor

/-- The divisibility step in the corrected relation-projection statement.
`k` is the number of rational kernel points. The hypotheses `k ∣ h*r`
and `k ∣ d` come from Lagrange and the isogeny degree, respectively. -/
theorem rationalKernelOrderDvdCofactor
    (h r d k : Nat)
    (kernelSubgroup : k ∣ h * r)
    (kernelOfIsogeny : k ∣ d)
    (coprimeTargetDegree : Nat.Coprime r d) : k ∣ h := by
  have coprimeKernelTarget : Nat.Coprime k r :=
    coprimeTargetDegree.symm.coprime_dvd_left kernelOfIsogeny
  exact coprimeKernelTarget.dvd_of_dvd_mul_right kernelSubgroup

/-- The concrete witness's full degree is not a divisor of its cofactor. -/
theorem degreeThreeNotDvdCofactorOne : ¬ (3 : Nat) ∣ 1 := by
  decide

#print axioms FiberCofactor.rationalKernelOrderDvdCofactor
#print axioms FiberCofactor.degreeThreeNotDvdCofactorOne

end FiberCofactor
