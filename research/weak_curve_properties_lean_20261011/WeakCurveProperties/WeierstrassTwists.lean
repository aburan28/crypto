import Mathlib.AlgebraicGeometry.EllipticCurve.NormalForms
import Mathlib.Tactic.FieldSimp

set_option autoImplicit false

namespace WeakCurveProperties

variable {F : Type*} [Field F] [NeZero (2 : F)] [NeZero (3 : F)]

def shortModel (a b : F) : WeierstrassCurve F := ⟨0, 0, 0, a, b⟩

instance shortModel_isShortNF (a b : F) : (shortModel a b).IsShortNF :=
  ⟨rfl, rfl, rfl⟩

/-- Every admissible change between short models has zero translation/shear. -/
theorem short_change_has_zero_rst (a b a' b' : F)
    (C : WeierstrassCurve.VariableChange F)
    (heq : C • shortModel a b = shortModel a' b') :
    C.r = 0 ∧ C.s = 0 ∧ C.t = 0 := by
  have hu : ((C.u : F)⁻¹) ≠ 0 := inv_ne_zero (Units.ne_zero C.u)
  have hs := congrArg WeierstrassCurve.a₁ heq
  simp only [WeierstrassCurve.variableChange_a₁, shortModel, zero_add,
    Units.val_inv_eq_inv_val] at hs
  have hs0 : C.s = 0 := (mul_eq_zero.mp ((mul_eq_zero.mp hs).resolve_left hu)).resolve_left
    (NeZero.ne (2 : F))
  have hr := congrArg WeierstrassCurve.a₂ heq
  simp only [WeierstrassCurve.variableChange_a₂, shortModel, hs0, zero_mul,
    zero_add, zero_pow (by decide : 2 ≠ 0), sub_zero,
    Units.val_inv_eq_inv_val] at hr
  have hr0 : C.r = 0 :=
    (mul_eq_zero.mp ((mul_eq_zero.mp hr).resolve_left (pow_ne_zero 2 hu))).resolve_left
      (NeZero.ne (3 : F))
  have ht := congrArg WeierstrassCurve.a₃ heq
  simp only [WeierstrassCurve.variableChange_a₃, shortModel, mul_zero,
    zero_add, Units.val_inv_eq_inv_val] at ht
  have ht0 : C.t = 0 :=
    (mul_eq_zero.mp ((mul_eq_zero.mp ht).resolve_left (pow_ne_zero 3 hu))).resolve_left
      (NeZero.ne (2 : F))
  exact ⟨hr0, hs0, ht0⟩

/-- Full admissible-change classification of short models in characteristic >3. -/
theorem short_isomorphism_iff (a b a' b' : F) :
    (∃ C : WeierstrassCurve.VariableChange F,
      C • shortModel a b = shortModel a' b') ↔
    ∃ u : Fˣ, a = (u : F)^4*a' ∧ b = (u : F)^6*b' := by
  constructor
  · rintro ⟨C, hC⟩
    obtain ⟨hr, hs, ht⟩ := short_change_has_zero_rst a b a' b' C hC
    have ha := congrArg WeierstrassCurve.a₄ hC
    have hb := congrArg WeierstrassCurve.a₆ hC
    simp [WeierstrassCurve.variableChange_a₄, shortModel, hr, hs, ht] at ha
    simp [WeierstrassCurve.variableChange_a₆, shortModel, hr, ht] at hb
    have ha' : ((C.u : F)^4)⁻¹*a = a' := by simpa using ha
    have hb' : ((C.u : F)^6)⁻¹*b = b' := by simpa using hb
    refine ⟨C.u, ?_, ?_⟩
    · exact (inv_mul_eq_iff_eq_mul₀ (pow_ne_zero 4 (Units.ne_zero C.u))).mp ha'
    · exact (inv_mul_eq_iff_eq_mul₀ (pow_ne_zero 6 (Units.ne_zero C.u))).mp hb'
  · rintro ⟨u,ha,hb⟩
    refine ⟨⟨u,0,0,0⟩, ?_⟩
    ext <;> simp [shortModel, WeierstrassCurve.variableChange_a₁,
      WeierstrassCurve.variableChange_a₂, WeierstrassCurve.variableChange_a₃,
      WeierstrassCurve.variableChange_a₄, WeierstrassCurve.variableChange_a₆,
      ha, hb]

def genericTwistModel (a0 b0 : F) (D : Fˣ) : WeierstrassCurve F :=
  shortModel (a0*(D : F)^2) (b0*(D : F)^3)

/-- The quadratic power class is equivalent to actual admissible changes of
the normalized generic model, with no assumed isomorphism conclusion. -/
theorem generic_twist_isomorphism_iff (a0 b0 : F)
    (ha0 : a0 ≠ 0) (hb0 : b0 ≠ 0) (D D' : Fˣ) :
    (∃ C : WeierstrassCurve.VariableChange F,
      C • genericTwistModel a0 b0 D = genericTwistModel a0 b0 D') ↔
    ∃ u : Fˣ, (D : F) = (u : F)^2*(D' : F) := by
  rw [genericTwistModel, genericTwistModel, short_isomorphism_iff]
  constructor
  · rintro ⟨u, ha, hb⟩
    have ha' : (D : F)^2 = (u : F)^4*(D' : F)^2 :=
      mul_left_cancel₀ ha0 (by simpa [mul_assoc, mul_left_comm, mul_comm] using ha)
    have hb' : (D : F)^3 = (u : F)^6*(D' : F)^3 :=
      mul_left_cancel₀ hb0 (by simpa [mul_assoc, mul_left_comm, mul_comm] using hb)
    refine ⟨u, ?_⟩
    calc
      (D : F) = (D : F)^3/(D : F)^2 := by
        field_simp
      _ = ((u : F)^6*(D' : F)^3)/((u : F)^4*(D' : F)^2) := by rw [ha', hb']
      _ = (u : F)^2*(D' : F) := by
        field_simp
  · rintro ⟨u,hu⟩
    refine ⟨u, ?_, ?_⟩ <;> rw [hu] <;> ring

theorem sextic_twist_isomorphism_iff (b b' : F) :
    (∃ C : WeierstrassCurve.VariableChange F,
      C • shortModel 0 b = shortModel 0 b') ↔
    ∃ u : Fˣ, b = (u : F)^6*b' := by
  simp only [short_isomorphism_iff, mul_zero, true_and]

theorem quartic_twist_isomorphism_iff (a a' : F) :
    (∃ C : WeierstrassCurve.VariableChange F,
      C • shortModel a 0 = shortModel a' 0) ↔
    ∃ u : Fˣ, a = (u : F)^4*a' := by
  simp only [short_isomorphism_iff, mul_zero, and_true]

def shortJ (a b : F) : F := 6912*a^3/(4*a^3+27*b^2)

omit [NeZero (2 : F)] [NeZero (3 : F)] in
theorem short_j_agrees_with_weierstrass (a b : F) [(shortModel a b).IsElliptic] :
    (shortModel a b).j = shortJ a b := by
  simpa only [shortJ, shortModel] using
    WeierstrassCurve.j_of_isShortNF (shortModel a b)

/-- The normalizing coefficients in the report follow over every field of
characteristic different from 2 and 3, with the generic denominators exposed. -/
theorem generic_normalization (a b : F) (ha : a ≠ 0) (hb : b ≠ 0)
    (hdisc : 4*a^3+27*b^2 ≠ 0) :
    let j := shortJ a b
    let D := 3*b/(2*a)
    j ≠ 0 ∧ 1728-j ≠ 0 ∧ D ≠ 0 ∧
      a = D^2*(3*j/(1728-j)) ∧ b = D^3*(2*j/(1728-j)) := by
  dsimp only
  have h2 : (2 : F) ≠ 0 := NeZero.ne (2 : F)
  have h3 : (3 : F) ≠ 0 := NeZero.ne (3 : F)
  have h6912 : (6912 : F) ≠ 0 := by
    convert mul_ne_zero (pow_ne_zero 8 h2) (pow_ne_zero 3 h3) using 1
    norm_num
  have h46656 : (46656 : F) ≠ 0 := by
    convert mul_ne_zero (pow_ne_zero 6 h2) (pow_ne_zero 6 h3) using 1
    norm_num
  have hj : shortJ a b ≠ 0 := div_ne_zero
    (mul_ne_zero h6912 (pow_ne_zero 3 ha)) hdisc
  have hjdisc : shortJ a b*(4*a^3+27*b^2) = 6912*a^3 :=
    div_mul_cancel₀ _ hdisc
  have hdiff : (1728-shortJ a b)*(4*a^3+27*b^2) = 46656*b^2 := by
    rw [sub_mul, hjdisc]
    ring
  have hj1728 : 1728-shortJ a b ≠ 0 := by
    intro h
    have hz : (46656 : F)*b^2 = 0 := by rw [← hdiff, h, zero_mul]
    exact mul_ne_zero h46656 (pow_ne_zero 2 hb) hz
  have hD : 3*b/(2*a) ≠ 0 := div_ne_zero
    (mul_ne_zero h3 hb) (mul_ne_zero h2 ha)
  refine ⟨hj, hj1728, hD, ?_, ?_⟩
  · rw [← mul_div_assoc, eq_div_iff hj1728]
    apply mul_right_cancel₀ hdisc
    simp only [mul_assoc, hdiff, hjdisc]
    field_simp [h2, ha]
    ring
  · rw [← mul_div_assoc, eq_div_iff hj1728]
    apply mul_right_cancel₀ hdisc
    simp only [mul_assoc, hdiff, hjdisc]
    field_simp [h2, ha]
    ring

end WeakCurveProperties
