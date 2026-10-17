import TraceProjection
import Mathlib.Algebra.Field.Basic

/-!
For a quotient image fixed by the third Frobenius power, its six-step
trace is twice its three-step trace. The field lemma below records the
actual rational expression used to recover the first cover coordinate
and map it to the elliptic quotient. The P-192 receipt checks that the
model constants, input plane coordinates, and inverse abscissa lie in the
fixed subfield. The elliptic-curve coordinate-to-point bridge remains
an explicit hypothesis of the additive-group result.
-/

namespace P192ExtensionCover

variable {K : Type*} [AddCommGroup K]

def threeTrace (f : K →+ K) (P : K) : K := P + f P + f (f P)

theorem sixTrace_eq_two_threeTrace (f : K →+ K) (P : K)
    (h3 : f (f (f P)) = P) :
    sixTrace f P = (2 : ℕ) • threeTrace f P := by
  simp only [sixTrace, threeTrace, h3]
  abel

theorem threeTrace_cofactor_projection
    (f : K →+ K) (P : K) (h k : ℕ)
    (h3 : f (f (f P)) = P)
    (hprojection : k • sixTrace f P = h • P) :
    (2 * k) • threeTrace f P = h • P := by
  rw [sixTrace_eq_two_threeTrace f P h3] at hprojection
  simpa only [mul_nsmul, nsmul_left_comm] using hprojection

/-- The three-step projector follows directly from the checked six-step
cofactor hypotheses and the subfield condition, without assuming a
projection equation as a separate input. -/
theorem threeTrace_cofactor_projection_from_conditions
    (f : K →+ K) (P : K) (r h k u m v : ℕ)
    (h3 : f (f (f P)) = P)
    (hfix : f (h • P) = h • P)
    (hrTrace : r • sixTrace f P = 0)
    (hrQ : r • (h • P) = 0)
    (hk : h = 6 * k + r * m)
    (hu : 6 * u = 1 + r * v) :
    (2 * k) • threeTrace f P = h • P := by
  exact threeTrace_cofactor_projection f P h k h3
    (sixTrace_cofactor_projection f P r h k u m v
      hfix hrTrace hrQ hk hu)

variable {F : Type*} [Field F]

/-- Every field operation in the recovered first cover coordinate and
the quotient image preserves the fixed subfield of the cube Frobenius.
The five nonzero hypotheses are the inverse's open-locus conditions. -/
theorem recovered_quotient_image_cube_fixed
    (σ : F →+* F) (A0 A1 A2 f0 norm X rho c0 d0 kappa : F)
    (hA0 : σ A0 = A0) (hA1 : σ A1 = A1) (hA2 : σ A2 = A2)
    (hf0 : σ f0 = f0) (hnorm : σ norm = norm)
    (hX : σ X = X) (hrho : σ rho = rho)
    (hc0 : σ c0 = c0) (hd0 : σ d0 = d0)
    (hkappa : σ kappa = kappa)
    (_hA0ne : A0 ≠ 0) (_hA1ne : A1 ≠ 0)
    (_hA2ne : A2 ≠ 0) (_hnormne : norm ≠ 0)
    (_hXne : X ≠ rho) :
    let lam := A0 * f0 / (A1 * A2 * norm)
    let y0 := 1 / (lam * A0)
    σ (c0 + d0 / (X - rho)) = c0 + d0 / (X - rho) ∧
      σ (kappa * y0 / (X - rho) ^ 2) =
        kappa * y0 / (X - rho) ^ 2 := by
  dsimp
  constructor <;> simp [hA0, hA1, hA2, hf0, hnorm, hX, hrho,
    hc0, hd0, hkappa]

end P192ExtensionCover

#print axioms P192ExtensionCover.sixTrace_eq_two_threeTrace
#print axioms P192ExtensionCover.threeTrace_cofactor_projection
#print axioms P192ExtensionCover.threeTrace_cofactor_projection_from_conditions
#print axioms P192ExtensionCover.recovered_quotient_image_cube_fixed
