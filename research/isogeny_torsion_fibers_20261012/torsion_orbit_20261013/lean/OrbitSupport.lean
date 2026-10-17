import Init

/-!
An abstract quotient-support theorem. `addE` may be elliptic-curve addition,
`phi` either the isogeny or the subgroup projection, and `kernel` the rational
three-torsion kernel. The curve-level homomorphism and kernel facts remain
external assumptions; this theorem proves their exact support consequence.
-/

universe u v

theorem torsion_orbit_pair_support
    {E : Type u} {D : Type v}
    (addE : E → E → E) (addD : D → D → D) (zeroD : D)
    (phi : E → D) (base kernel : E → Prop)
    (hhom : ∀ a b, phi (addE a b) = addD (phi a) (phi b))
    (hzero : ∀ a, addD a zeroD = a)
    (hkernel : ∀ t, kernel t → phi t = zeroD)
    (hidentity : ∃ t, kernel t ∧ ∀ a, addE a t = a)
    (q : D) :
    (∃ a b t s, base a ∧ base b ∧ kernel t ∧ kernel s ∧
      phi (addE (addE a t) (addE b s)) = q) ↔
    (∃ a b, base a ∧ base b ∧ phi (addE a b) = q) := by
  constructor
  · rintro ⟨a, b, t, s, ha, hb, ht, hs, heq⟩
    refine ⟨a, b, ha, hb, ?_⟩
    have hcollapse :
        phi (addE (addE a t) (addE b s)) = phi (addE a b) := by
      rw [hhom, hhom, hhom, hhom, hkernel t ht, hkernel s hs]
      rw [hzero, hzero]
    rw [hcollapse] at heq
    exact heq
  · rintro ⟨a, b, ha, hb, heq⟩
    obtain ⟨t, ht, hid⟩ := hidentity
    refine ⟨a, b, t, t, ha, hb, ht, ht, ?_⟩
    rw [hid a, hid b]
    exact heq

#print axioms torsion_orbit_pair_support
