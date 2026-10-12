import WeakCurveProperties.CyclicTwists
import WeakCurveProperties.WeierstrassTwists
import Mathlib.RingTheory.IntegralDomain
import Mathlib.Algebra.GroupWithZero.Units.Fintype

set_option autoImplicit false

namespace WeakCurveProperties

variable {F : Type*} [Field F] [Finite F]

/-- The power test is specialized to actual finite-field units. -/
theorem finite_field_unit_twist_test (S : Subgroup Fˣ) (e : ℕ) (x : Fˣ) :
    x ∈ S ⊔ (powMonoidHom e : Fˣ →* Fˣ).range ↔
      x^((Nat.card F-1) / Nat.gcd e ((Nat.card F-1) / Nat.card S)) = 1 := by
  have hc : Nat.card Fˣ = Nat.card F-1 := Nat.card_units F
  simpa only [hc] using finite_cyclic_twist_descent S e x

def subfieldUnitsSubgroup (L : Subfield F) : Subgroup Fˣ :=
  (Units.map L.subtype.toMonoidHom).range

omit [Finite F] in
theorem card_subfield_units (L : Subfield F) :
    Nat.card (subfieldUnitsSubgroup L) = Nat.card L-1 := by
  have hi : Function.Injective (Units.map L.subtype.toMonoidHom) := by
    intro x y h
    apply Units.ext
    apply Subtype.ext
    have he := congrArg (fun z : Fˣ => (z : F)) h
    exact he
  change Nat.card (Set.range (Units.map L.subtype.toMonoidHom)) = Nat.card L-1
  rw [Nat.card_range_of_injective hi]
  exact Nat.card_units L

/-- The report's exact subfield unit-image criterion, with M=|F|-1 and m=|L|-1.
No generator or assumed image membership is supplied as a hypothesis. -/
theorem finite_field_subfield_twist_test (L : Subfield F) (e : ℕ) (x : Fˣ) :
    x ∈ subfieldUnitsSubgroup L ⊔ (powMonoidHom e : Fˣ →* Fˣ).range ↔
      x^((Nat.card F-1) / Nat.gcd e ((Nat.card F-1)/(Nat.card L-1))) = 1 := by
  simpa only [card_subfield_units] using finite_field_unit_twist_test
    (subfieldUnitsSubgroup L) e x

omit [Finite F] in
theorem subfield_twist_witness_iff (L : Subfield F) (e : ℕ) (x : Fˣ) :
    (∃ d : Lˣ, ∃ u : Fˣ, x = Units.map L.subtype.toMonoidHom d * u^e) ↔
      x ∈ subfieldUnitsSubgroup L ⊔ (powMonoidHom e : Fˣ →* Fˣ).range := by
  constructor
  · rintro ⟨d,u,hx⟩
    exact Subgroup.mem_sup.mpr ⟨_, ⟨d,rfl⟩, _, ⟨u,rfl⟩, hx.symm⟩
  · intro hx
    obtain ⟨s,hs,v,hv,heq⟩ := Subgroup.mem_sup.mp hx
    rcases hs with ⟨d,rfl⟩
    rcases hv with ⟨u,rfl⟩
    exact ⟨d,u,heq.symm⟩

/-- A construction witness exists exactly when the subfield power test passes. -/
theorem subfield_twist_witness_power_test (L : Subfield F) (e : ℕ) (x : Fˣ) :
    (∃ d : Lˣ, ∃ u : Fˣ, x = Units.map L.subtype.toMonoidHom d * u^e) ↔
      x^((Nat.card F-1) / Nat.gcd e ((Nat.card F-1)/(Nat.card L-1))) = 1 := by
  rw [subfield_twist_witness_iff, finite_field_subfield_twist_test]

omit [Finite F] in
theorem mapped_generic_model (L : Subfield F) (a0 b0 : L) (d : Lˣ) :
    (genericTwistModel a0 b0 d).map L.subtype =
      genericTwistModel (a0 : F) (b0 : F) (Units.map L.subtype.toMonoidHom d) := by
  ext <;> simp [genericTwistModel, shortModel, WeierstrassCurve.map]

def normalizedGenericDescends (L : Subfield F) (a0 b0 : L) (D : Fˣ) : Prop :=
  ∃ d : Lˣ, ∃ C : WeierstrassCurve.VariableChange F,
    C • genericTwistModel (a0 : F) (b0 : F) D =
      (genericTwistModel a0 b0 d).map L.subtype

variable [NeZero (2 : F)] [NeZero (3 : F)]

/-- Descent here means an actual admissible change to a model over L in the
normalized generic family. The power test supplies both the subfield parameter
and the change of variables. -/
theorem normalized_generic_descent_iff (L : Subfield F) (a0 b0 : L)
    (ha0 : a0 ≠ 0) (hb0 : b0 ≠ 0) (D : Fˣ) :
    normalizedGenericDescends L a0 b0 D ↔
      D^((Nat.card F-1) / Nat.gcd 2 ((Nat.card F-1)/(Nat.card L-1))) = 1 := by
  have haF : (a0 : F) ≠ 0 := fun h => ha0 (Subtype.ext h)
  have hbF : (b0 : F) ≠ 0 := fun h => hb0 (Subtype.ext h)
  rw [← subfield_twist_witness_power_test L 2 D]
  constructor
  · rintro ⟨d,C,hC⟩
    rw [mapped_generic_model] at hC
    obtain ⟨u,hu⟩ := (generic_twist_isomorphism_iff _ _ haF hbF D _).mp ⟨C,hC⟩
    refine ⟨d,u,?_⟩
    apply Units.ext
    simpa only [Units.val_mul, Units.val_pow_eq_pow_val, mul_comm] using hu
  · rintro ⟨d,u,hu⟩
    refine ⟨d,?_⟩
    rw [mapped_generic_model]
    apply (generic_twist_isomorphism_iff _ _ haF hbF D _).mpr
    refine ⟨u,?_⟩
    have hc := congrArg (fun x : Fˣ => (x : F)) hu
    simpa only [Units.val_mul, Units.val_pow_eq_pow_val, mul_comm] using hc

end WeakCurveProperties
