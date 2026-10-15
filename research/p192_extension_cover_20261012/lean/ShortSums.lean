import AdaptiveBase

/-!
Short-sum collision implications for the adaptive Type I collector. These
are abstract group identities. The concrete P-192 absence of collisions is
an independently replayed finite computation, not a Lean finite-field proof.
-/

namespace P192ExtensionCover

variable {K : Type*} [AddCommGroup K]

/-- A second two-summand representation eliminates the fresh point from a
    known-target adaptive equation. -/
theorem two_sum_collision_eliminates_fresh
    (a b c d e q : K) (hadaptive : a + b + c = q)
    (hsecond : c = d + e) : a + b + d + e = q := by
  rw [hsecond] at hadaptive
  convert hadaptive using 1 <;> abel

/-- The analogous four-summand representation gives a six-summand equation
    entirely in the original classes. -/
theorem four_sum_collision_eliminates_fresh
    (a b c d e f g q : K) (hadaptive : a + b + c = q)
    (hsecond : c = d + e + f + g) :
    a + b + d + e + f + g = q := by
  rw [hsecond] at hadaptive
  convert hadaptive using 1 <;> abel

end P192ExtensionCover

#print axioms P192ExtensionCover.two_sum_collision_eliminates_fresh
#print axioms P192ExtensionCover.four_sum_collision_eliminates_fresh
