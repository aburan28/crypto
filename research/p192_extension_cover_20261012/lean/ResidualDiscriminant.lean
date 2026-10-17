import Mathlib.Tactic.Ring

/-!
Algebraic obligations used by the complete P-192 quartic-line residual
screen. The concrete finite-field quartic and its 49,247,850 pair inputs
are checked by the separate coordinate replays.
-/

namespace P192ExtensionCover

theorem quartic_two_known_roots_factor_identity
    {R : Type*} [CommRing R] (t d A B C : R) :
    t * (t - d) * (A * t ^ 2 + B * t + C) =
      A * t ^ 4 + (B - d * A) * t ^ 3 +
      (C - d * B) * t ^ 2 - d * C * t := by
  ring

theorem residual_quadratic_root_is_simple
    {R : Type*} [CommRing R] (A B C x : R)
    (hdisc : B ^ 2 - (4 : R) * A * C ≠ 0)
    (hroot : A * x ^ 2 + B * x + C = 0) :
    (2 : R) * A * x + B ≠ 0 := by
  intro hderiv
  apply hdisc
  calc
    B ^ 2 - (4 : R) * A * C =
        ((2 : R) * A * x + B) ^ 2 -
          (4 : R) * A * (A * x ^ 2 + B * x + C) := by ring
    _ = 0 := by rw [hderiv, hroot]; ring

theorem residual_root_avoids_known_points
    {R : Type*} [CommRing R] (A B C d x : R)
    (hconstant : C ≠ 0)
    (hsecond : A * d ^ 2 + B * d + C ≠ 0)
    (hdisc : B ^ 2 - (4 : R) * A * C ≠ 0)
    (hroot : A * x ^ 2 + B * x + C = 0) :
    x ≠ 0 ∧ x ≠ d ∧ (2 : R) * A * x + B ≠ 0 := by
  refine ⟨?_, ?_, residual_quadratic_root_is_simple A B C x hdisc hroot⟩
  · intro hx
    exact hconstant (by simpa [hx] using hroot)
  · intro hx
    exact hsecond (by simpa [hx] using hroot)

end P192ExtensionCover

#print axioms P192ExtensionCover.quartic_two_known_roots_factor_identity
#print axioms P192ExtensionCover.residual_quadratic_root_is_simple
#print axioms P192ExtensionCover.residual_root_avoids_known_points
