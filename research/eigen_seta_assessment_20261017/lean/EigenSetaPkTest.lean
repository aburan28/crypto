/-!
# Arithmetic behind the public-key degree test for Eigen-Séta v2

Report: `research/eigen_seta_assessment_20261017/REPORT.md`, Lemma 1 and
Lemma 2. Core Lean only (no Mathlib).

Lemma 1 (`pkTest_kernel_isotropic`). Let `M = 2^e`, let `s` satisfy
`s^2 ≡ -7 (mod M)`, and let `c^2 + c'^2 + 7 u^2 f^2 = M`. With
`λ_S = u f s + c` and `λ_T = c - u f s` we have
`c'^2 + λ_S λ_T ≡ 0 (mod M)`. This is the Weil-pairing condition that makes
`⟨(c'P, -λ_S P), (c'Q, -λ_T Q)⟩` isotropic, for every candidate `f`
(true or wrong): the test chain always runs to its end.

Lemma 2 (`pkTest_kani_degree`). In `ℤ[f√-7]`, written as pairs `(c, w)` for
`c + w·f√-7`, the norm of `c + u·f√-7` is `c^2 + 7 u^2 f^2`, so
`deg γ + c'^2 = 2^e` exactly when `c^2 + c'^2 + 7 u^2 f^2 = 2^e`. The bridge
from this norm to the degree of the endomorphism `u ψ_A + c` of `E_A`
(`deg = nrd`, `ψ_A^2 = -7 f^2`) is not formalised: Lean core has no
elliptic curves or quaternion orders.
-/

theorem pkTest_kernel_isotropic (M s u f c c' k : Int)
    (hs : s * s + 7 = M * k)
    (hk : c * c + c' * c' + 7 * (u * u * (f * f)) = M) :
    ∃ q : Int, c' * c' + (u * f * s + c) * (c - u * f * s) = M * q := by
  refine ⟨1 - u * u * (f * f) * k, ?_⟩
  have e1 : (u * f * s + c) * (c - u * f * s)
      = c * c - u * u * (f * f) * (s * s) := by
    simp only [Int.add_mul, Int.mul_sub]
    simp only [Int.mul_comm, Int.mul_left_comm]
    omega
  have e2 : u * u * (f * f) * (s * s) = u * u * (f * f) * (M * k) - 7 * (u * u * (f * f)) := by
    have : s * s = M * k - 7 := by omega
    rw [this, Int.mul_sub]
    simp only [Int.mul_comm]
  rw [e1, e2, Int.mul_sub, Int.mul_one]
  have e3 : u * u * (f * f) * (M * k) = M * (u * u * (f * f) * k) := by
    simp only [Int.mul_comm, Int.mul_left_comm]
  rw [e3]
  omega

/-- Norm on pairs `(c, w)` representing `c + w·f√-7`. -/
def normF (f : Int) (x : Int × Int) : Int := x.1 * x.1 + 7 * (f * f) * (x.2 * x.2)

theorem pkTest_kani_degree (e u f c c' : Int) :
    normF f (c, u) + c' * c' = e ↔ c * c + c' * c' + 7 * (u * u * (f * f)) = e := by
  unfold normF
  have : 7 * (f * f) * (u * u) = 7 * (u * u * (f * f)) := by
    simp only [Int.mul_comm, Int.mul_assoc, Int.mul_left_comm]
  simp only [this]
  omega
