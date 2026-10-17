import Init.Grind.Ring.Basic

namespace SixClassCoefficients

/-- If `sigma` is additive, respects scalar multiplication, and acts by
`lambda` on both public generators, then it updates both rho coefficients
by the same multiplier. This is conditional on the actual curve map having
those properties; it does not formalize the curve group law. -/
theorem endomorphismWalkCoefficients {R P : Type} [Mul R]
    (add : P → P → P) (smul : R → P → P) (sigma : P → P)
    (lambda a b : R) (G Q : P)
    (hadd : ∀ x y, sigma (add x y) = add (sigma x) (sigma y))
    (hsmul : ∀ c x, sigma (smul c x) = smul c (sigma x))
    (hassoc : ∀ c d x, smul c (smul d x) = smul (c * d) x)
    (hG : sigma G = smul lambda G)
    (hQ : sigma Q = smul lambda Q) :
    sigma (add (smul a G) (smul b Q)) =
      add (smul (a * lambda) G) (smul (b * lambda) Q) := by
  rw [hadd, hsmul, hsmul, hG, hQ, hassoc, hassoc]

/-- The scalar value represented by rho's coefficient pair is preserved
when the point and both coefficients are multiplied by the same scalar. -/
theorem scaledCoefficientValue {R : Type} [Lean.Grind.CommRing R]
    (lambda a b d : R) :
    a * lambda + (b * lambda) * d = lambda * (a + b * d) := by
  grind

/-- Negating a walk point negates both coefficients. -/
theorem negatedCoefficientValue {R : Type} [Lean.Grind.CommRing R]
    (a b d : R) :
    -a + (-b) * d = -(a + b * d) := by
  grind

#print axioms SixClassCoefficients.endomorphismWalkCoefficients
#print axioms SixClassCoefficients.scaledCoefficientValue
#print axioms SixClassCoefficients.negatedCoefficientValue

end SixClassCoefficients
