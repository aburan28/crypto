import SignedBase

/-!
The finite P-192 support census is checked by native arithmetic and an
independent signed-multiset count. These abstract equations identify exactly
what an anchor or a repeated target would contribute to the logarithm system.
-/

namespace P192ExtensionCover

variable {G K : Type*} [AddCommGroup G] [AddCommGroup K]

/-- A collision between a rank-bearing support point and a known-scalar
    anchor gives a directly usable equation in the unknown contribution. -/
theorem scalar_anchor_collision_equation (log : G →+ K)
    (unknown known anchor : G) (knownLog anchorLog : K)
    (hknown : log known = knownLog) (hanchor : log anchor = anchorLog)
    (hcollision : unknown + known = anchor) :
    log unknown = anchorLog - knownLog := by
  have h := congrArg log hcollision
  simp only [map_add, hknown, hanchor] at h
  calc
    log unknown = (log unknown + knownLog) - knownLog := by abel
    _ = anchorLog - knownLog := by rw [h]

/-- Two distinct formal representations of one target give a relation
    through the difference of their unknown and known contributions. -/
theorem repeated_target_relation_equation (log : G →+ K)
    (unknown₁ unknown₂ known₁ known₂ : G)
    (knownLog₁ knownLog₂ : K)
    (hknown₁ : log known₁ = knownLog₁)
    (hknown₂ : log known₂ = knownLog₂)
    (hcollision : unknown₁ + known₁ = unknown₂ + known₂) :
    log unknown₁ - log unknown₂ = knownLog₂ - knownLog₁ := by
  have h := congrArg log hcollision
  simp only [map_add, hknown₁, hknown₂] at h
  calc
    log unknown₁ - log unknown₂ =
        (log unknown₁ + knownLog₁) - (log unknown₂ + knownLog₂) +
          (knownLog₂ - knownLog₁) := by abel
    _ = knownLog₂ - knownLog₁ := by rw [h]; abel

end P192ExtensionCover

#print axioms P192ExtensionCover.scalar_anchor_collision_equation
#print axioms P192ExtensionCover.repeated_target_relation_equation
