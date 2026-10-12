import Mathlib.Algebra.Group.Subgroup.Basic
import Mathlib.Algebra.Group.Hom.Instances

set_option autoImplicit false

namespace WeakCurveProperties

variable {A : Type*} [AddCommGroup A]

def mulBy (n : ℕ) : A →+ A where
  toFun := fun x => n • x
  map_zero' := nsmul_zero n
  map_add' := fun x y => nsmul_add x y n

def fixesNTorsion (F : A →+ A) (n : ℕ) : Prop :=
  ∀ x : A, n • x = 0 → F x = x

/-- At the point-group level, fixing n-torsion is exactly divisibility of F-id
through the surjective n-multiplication map. Algebraicity of the resulting
homomorphism is a separate geometric obligation in the proof map. -/
theorem fixes_torsion_iff_factors_difference (F : A →+ A) (n : ℕ)
    (hsurj : Function.Surjective (mulBy (A := A) n)) :
    fixesNTorsion F n ↔
      ∃ psi : A →+ A, psi.comp (mulBy n) = F - AddMonoidHom.id A := by
  constructor
  · intro hfix
    have hker : (mulBy (A := A) n).ker ≤ (F - AddMonoidHom.id A).ker := by
      intro x hx
      change F x - x = 0
      exact sub_eq_zero.mpr (hfix x hx)
    refine ⟨(mulBy n).liftOfSurjective hsurj ⟨F - AddMonoidHom.id A, hker⟩, ?_⟩
    ext x
    exact (mulBy (A := A) n).liftOfRightInverse_comp_apply
      (Function.surjInv hsurj) (Function.rightInverse_surjInv hsurj)
      ⟨F - AddMonoidHom.id A, hker⟩ x
  · rintro ⟨psi, hpsi⟩ x hx
    have heq := DFunLike.congr_fun hpsi x
    change psi (n • x) = F x - x at heq
    rw [hx, map_zero] at heq
    exact sub_eq_zero.mp heq.symm

end WeakCurveProperties
