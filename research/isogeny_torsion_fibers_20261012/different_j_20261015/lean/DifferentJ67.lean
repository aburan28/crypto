import Init

namespace DifferentJ67

/-- The two measured j residues differ in the frozen base field. -/
theorem distinctJ : (7 : Fin 67) ≠ 11 := by decide

/-- The degree-three kernel order is coprime to the target subgroup order. -/
theorem targetCoprimeKernel : Nat.Coprime 61 3 := by decide

/-- Conditional scalar version of kernel injectivity on the order-61 subgroup. -/
theorem kernelScalarVanishes :
    ∀ k : Fin 61, (3 * k : Fin 61) = 0 → k = 0 := by
  decide

/-- An image map is injective on the order-61 scalar subgroup if every
    collision has difference killed by the degree-three kernel. The bridge
    from this premise to the concrete Vélu map is still to be formalized. -/
theorem injectiveOfThreeKernel {β : Type} (image : Fin 61 → β)
    (kernel : ∀ a b : Fin 61, image a = image b → (3 * (a - b) : Fin 61) = 0) :
    Function.Injective image := by
  intro a b equalImage
  have differenceZero := kernelScalarVanishes (a - b) (kernel a b equalImage)
  omega

/-- Arithmetic check of the recorded rational point and image counts. -/
theorem recordedFiberCardinality : 4575 = 3 * 1525 := by decide

/-- Arithmetic check of the reported subgroup and cofactor. -/
theorem recordedSubgroupOrder : 4575 = 75 * 61 := by decide

#print axioms DifferentJ67.distinctJ
#print axioms DifferentJ67.targetCoprimeKernel
#print axioms DifferentJ67.kernelScalarVanishes
#print axioms DifferentJ67.injectiveOfThreeKernel
#print axioms DifferentJ67.recordedFiberCardinality
#print axioms DifferentJ67.recordedSubgroupOrder

end DifferentJ67
