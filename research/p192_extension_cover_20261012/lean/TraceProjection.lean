import Mathlib.Tactic.Abel
import Mathlib.Tactic.NormNum

/-!
The six-step Frobenius trace replaces a large extension-field cofactor
multiplication by a prime-field scalar multiplication.  The group-theoretic
statement keeps the actual elliptic curve and Frobenius hypotheses explicit:
the group order and rational subgroup order are checked in the P-192 receipt.
-/

namespace P192ExtensionCover

variable {K : Type*} [AddCommGroup K]

/-- Sum of the six Galois conjugates under one additive Frobenius map. -/
def sixTrace (f : K →+ K) (P : K) : K :=
  P + f P + f (f P) + f (f (f P)) +
    f (f (f (f P))) + f (f (f (f (f P))))

theorem sixTrace_nsmul (f : K →+ K) (n : ℕ) (P : K) :
    sixTrace f (n • P) = n • sixTrace f P := by
  simp [sixTrace, nsmul_add]

theorem sixTrace_fixed (f : K →+ K) (Q : K) (hfix : f Q = Q) :
    sixTrace f Q = (6 : ℕ) • Q := by
  simp only [sixTrace, hfix]
  abel

/-- Abstract trace projection: if the cofactor component is Frobenius-fixed,
    the six-step trace and a scalar congruent to h/6 modulo r recover it. -/
theorem sixTrace_cofactor_projection
    (f : K →+ K) (P : K) (r h k u m v : ℕ)
    (hfix : f (h • P) = h • P)
    (hrTrace : r • sixTrace f P = 0)
    (hrQ : r • (h • P) = 0)
    (hk : h = 6 * k + r * m)
    (hu : 6 * u = 1 + r * v) :
    k • sixTrace f P = h • P := by
  let T := sixTrace f P
  let Q := h • P
  have hTrace : h • T = (6 : ℕ) • Q := by
    change h • sixTrace f P = (6 : ℕ) • (h • P)
    rw [← sixTrace_nsmul, sixTrace_fixed f (h • P) hfix]
  have hSix : (6 : ℕ) • (k • T - Q) = 0 := by
    rw [nsmul_sub]
    have hleft : (6 : ℕ) • (k • T) = h • T := by
      have hrhs : h • T = (6 : ℕ) • (k • T) := by
        rw [hk, add_nsmul, mul_nsmul', mul_nsmul, hrTrace]
        simp
      exact hrhs.symm
    rw [hleft, hTrace, sub_self]
  have hR : r • (k • T - Q) = 0 := by
    rw [nsmul_sub, nsmul_left_comm, hrTrace, nsmul_zero, hrQ, sub_self]
  have hUnit : (6 * u) • (k • T - Q) = k • T - Q := by
    rw [hu, add_nsmul, mul_nsmul, hR]
    simp only [one_nsmul, nsmul_zero, add_zero]
  have hZero : (6 * u) • (k • T - Q) = 0 := by
    rw [mul_nsmul, hSix, nsmul_zero]
  have : k • T - Q = 0 := by rw [← hUnit, hZero]
  exact sub_eq_zero.mp this

private def p192R : ℕ :=
  6277101735386680763835789423176059013767194773182842284081

private def p192H : ℕ :=
  9745314011399999080353382387924259374680358052274105047343790102210948785313621242842357105607587396653353823515292585640417051742014748800772927851507971775934426802951755139147189236582634345051735952092818586369948498241481992249433898803432562525942393045320537614777930915594708057584

private def p192K : ℕ :=
  1740964300588857170778051573990638825050470570621855798409

private def p192U : ℕ :=
  1046183622564446793972631570529343168961199128863807047347

private def p192M : ℕ :=
  1552518092300708935148979488478137498247665834575054412699184799862574496911786977738779413886835400471016139625245869296159123915085802963714115974752757105509600886564046380812406380982252977281006291348943313700818729561470048730

/-- Exact arithmetic certificate for the scalar used in the checked Sage
    trace projector; in particular, 6 is invertible modulo the P-192 order. -/
theorem p192_trace_scalar_certificate :
    p192H = 6 * p192K + p192R * p192M ∧
    6 * p192U = 1 + p192R := by
  norm_num [p192H, p192K, p192R, p192M, p192U]

/-- A line-section identity determines the second residual column from the
    first and the two selected endpoint columns. -/
theorem line_partner_from_constant (S A B X Y : K)
    (hline : A + B + X + Y = S) : Y = S - A - B - X := by
  calc
    Y = (A + B + X + Y) - A - B - X := by abel
    _ = S - A - B - X := by rw [hline]

end P192ExtensionCover

#print axioms P192ExtensionCover.sixTrace_nsmul
#print axioms P192ExtensionCover.sixTrace_fixed
#print axioms P192ExtensionCover.sixTrace_cofactor_projection
#print axioms P192ExtensionCover.p192_trace_scalar_certificate
#print axioms P192ExtensionCover.line_partner_from_constant
