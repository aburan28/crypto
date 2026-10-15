import Init

namespace SemaevRootTransfer

/-- Once an additive image map is injective on the working subgroup,
    an image relation is equivalent to the source relation. This is the
    conditional algebraic bridge used by the quotient-root control. -/
theorem imageRelationIff {A B : Type} [Add A] [Add B]
    (image : A → B)
    (mapAdd : ∀ a b : A, image (a + b) = image a + image b)
    (injective : Function.Injective image) (a b t : A) :
    image a + image b = image t ↔ a + b = t := by
  constructor
  · intro imageRelation
    apply injective
    simpa only [mapAdd] using imageRelation
  · intro sourceRelation
    rw [← mapAdd a b, sourceRelation]

/-- Arithmetic kernel tests for the four frozen prime subgroup orders. -/
theorem kernelVanishes19 : ∀ k : Fin 19, (3 * k : Fin 19) = 0 → k = 0 := by decide
theorem kernelVanishes61 : ∀ k : Fin 61, (3 * k : Fin 61) = 0 → k = 0 := by decide
theorem kernelVanishes211 : ∀ k : Fin 211, (3 * k : Fin 211) = 0 → k = 0 := by
  intro k killed
  omega
theorem kernelVanishes811 : ∀ k : Fin 811, (3 * k : Fin 811) = 0 → k = 0 := by
  intro k killed
  omega

/-- If a frozen order-19 image collision differs by a point killed by
    the degree-three kernel, the image is injective on that subgroup. -/
theorem injectiveOfKernel19 {B : Type} (image : Fin 19 → B)
    (kernel : ∀ a b : Fin 19, image a = image b → (3 * (a - b) : Fin 19) = 0) :
    Function.Injective image := by
  intro a b sameImage
  have zero := kernelVanishes19 (a - b) (kernel a b sameImage)
  omega

theorem injectiveOfKernel61 {B : Type} (image : Fin 61 → B)
    (kernel : ∀ a b : Fin 61, image a = image b → (3 * (a - b) : Fin 61) = 0) :
    Function.Injective image := by
  intro a b sameImage
  have zero := kernelVanishes61 (a - b) (kernel a b sameImage)
  omega

theorem injectiveOfKernel211 {B : Type} (image : Fin 211 → B)
    (kernel : ∀ a b : Fin 211, image a = image b → (3 * (a - b) : Fin 211) = 0) :
    Function.Injective image := by
  intro a b sameImage
  have zero := kernelVanishes211 (a - b) (kernel a b sameImage)
  omega

theorem injectiveOfKernel811 {B : Type} (image : Fin 811 → B)
    (kernel : ∀ a b : Fin 811, image a = image b → (3 * (a - b) : Fin 811) = 0) :
    Function.Injective image := by
  intro a b sameImage
  have zero := kernelVanishes811 (a - b) (kernel a b sameImage)
  omega

#print axioms SemaevRootTransfer.imageRelationIff
#print axioms SemaevRootTransfer.kernelVanishes19
#print axioms SemaevRootTransfer.kernelVanishes61
#print axioms SemaevRootTransfer.kernelVanishes211
#print axioms SemaevRootTransfer.kernelVanishes811
#print axioms SemaevRootTransfer.injectiveOfKernel19
#print axioms SemaevRootTransfer.injectiveOfKernel61
#print axioms SemaevRootTransfer.injectiveOfKernel211
#print axioms SemaevRootTransfer.injectiveOfKernel811
