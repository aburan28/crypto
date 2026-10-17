import Init.Grind.Ring.Basic

namespace CubicAutomorphism

/-- The explicit field element used for the F13 curve pair has order three. -/
theorem zeta13 : 3 ^ 3 % 13 = 1 := by decide

/-- The explicit field element used for the F199 curve pair has order three. -/
theorem zeta199 : 92 ^ 3 % 199 = 1 := by decide

/-- The F13 destination parameter is a sixth-power scaling of the source. -/
theorem scaling13 : 2 ^ 6 * 2 % 13 = 11 := by decide

/-- The F199 destination parameter is a sixth-power scaling of the source. -/
theorem scaling199 : 14 ^ 6 * 3 % 199 = 118 := by decide

/-- A cubic root of unity preserves the affine equation of a j=0 curve.
This statement does not itself construct the elliptic-curve group law. -/
theorem cubicPreservesEquation {R : Type} [Lean.Grind.CommRing R]
    (zeta x y b : R) (hzeta : zeta ^ 3 = 1)
    (onCurve : y ^ 2 = x ^ 3 + b) : y ^ 2 = (zeta * x) ^ 3 + b := by
  calc
    y ^ 2 = x ^ 3 + b := onCurve
    _ = (zeta * x) ^ 3 + b := by grind

/-- Sixth-power scaling preserves the short Weierstrass equation while
changing `b` to `u^6*b`. Group-law preservation is a separate obligation. -/
theorem scalingPreservesEquation {R : Type} [Lean.Grind.CommRing R]
    (u x y b : R) (onCurve : y ^ 2 = x ^ 3 + b) :
    (u ^ 3 * y) ^ 2 = (u ^ 2 * x) ^ 3 + u ^ 6 * b := by
  grind

#print axioms CubicAutomorphism.zeta13
#print axioms CubicAutomorphism.zeta199
#print axioms CubicAutomorphism.scaling13
#print axioms CubicAutomorphism.scaling199
#print axioms CubicAutomorphism.cubicPreservesEquation
#print axioms CubicAutomorphism.scalingPreservesEquation

end CubicAutomorphism
