/-!
# Gauss-period windows are scaled polynomial windows (K1)

Let `γ` be a unit of a commutative ring with inverse `γ'` and put
`c k = γ^k + γ'^k`. For a type-II optimal normal basis, `γ` is a primitive
`(2n+1)`-th root of unity and the `c k` are the Gauss periods.

* `dickson_mul`: `c i · c j = c (i + j) + c (i − j)` for `j ≤ i`, in any
  commutative ring.
* `window_eq` (characteristic 2): the additive span of `c 1, …, c ℓ` equals
  the additive span of `c₁, c₁², …, c₁^ℓ`. So a Gauss-period window is
  `c₁ · ⟨1, c₁, …, c₁^{ℓ−1}⟩`: the coordinate window `C1` with generator
  `z = c₁`, scaled by `c₁`.

What this does and does not dispose of (README, register entry K1): it shows
the *set* is a scaled `C1` window, so it inherits that window's
product-dimension behaviour. It says nothing about using that set in the
2-torsion-invariant coordinate `s = x + 1/x`. That combination is untested
and stays open as experiment E-K1.

The additive span is over `ℤ`. In characteristic 2 this equals the
`F₂`-span.
-/

namespace FactorBaseDesigns

variable {α : Type} [Lean.Grind.CommRing α]

/-- `c k = γ^k + γ'^k`. -/
def gp (γ γ' : α) (k : Nat) : α := γ ^ k + γ' ^ k

theorem pow_mul_inv_pow (γ γ' : α) (h : γ * γ' = 1) (a j : Nat) :
    γ ^ (a + j) * γ' ^ j = γ ^ a := by
  induction j with
  | zero => simp [Lean.Grind.Semiring.pow_zero, Lean.Grind.Semiring.mul_one]
  | succ j ih =>
    have e1 : γ ^ (a + (j + 1)) = γ ^ (a + j) * γ := by
      rw [← Nat.add_assoc, Lean.Grind.Semiring.pow_succ]
    have e2 : γ' ^ (j + 1) = γ' ^ j * γ' := Lean.Grind.Semiring.pow_succ _ _
    rw [e1, e2]
    calc γ ^ (a + j) * γ * (γ' ^ j * γ') = (γ ^ (a + j) * γ' ^ j) * (γ * γ') := by grind
      _ = γ ^ a := by rw [ih, h]; grind

/-- **Dickson product identity.** `c i · c j = c (i+j) + c (i−j)` for `j ≤ i`. -/
theorem dickson_mul (γ γ' : α) (h : γ * γ' = 1) {i j : Nat} (hij : j ≤ i) :
    gp γ γ' i * gp γ γ' j = gp γ γ' (i + j) + gp γ γ' (i - j) := by
  have h' : γ' * γ = 1 := by grind
  have A : γ ^ i * γ' ^ j = γ ^ (i - j) := by
    have := pow_mul_inv_pow γ γ' h (i - j) j
    rwa [Nat.sub_add_cancel hij] at this
  have B : γ' ^ i * γ ^ j = γ' ^ (i - j) := by
    have := pow_mul_inv_pow γ' γ h' (i - j) j
    rwa [Nat.sub_add_cancel hij] at this
  have C : γ ^ (i + j) = γ ^ i * γ ^ j := Lean.Grind.Semiring.pow_add _ _ _
  have D : γ' ^ (i + j) = γ' ^ i * γ' ^ j := Lean.Grind.Semiring.pow_add _ _ _
  unfold gp
  rw [C, D]
  grind

/-- The additive (ℤ-) span of the elements satisfying `S`. -/
inductive ZSpan (S : α → Prop) : α → Prop
  | zero : ZSpan S 0
  | gen {x} : S x → ZSpan S x
  | add {x y} : ZSpan S x → ZSpan S y → ZSpan S (x + y)
  | neg {x} : ZSpan S x → ZSpan S (-x)

theorem ZSpan.mono {S T : α → Prop} (hST : ∀ x, S x → T x) {x : α}
    (h : ZSpan S x) : ZSpan T x := by
  induction h with
  | zero => exact .zero
  | gen hx => exact .gen (hST _ hx)
  | add _ _ ih1 ih2 => exact .add ih1 ih2
  | neg _ ih => exact .neg ih

/-- Multiplying a span element by `a` lands in the span of `a · S`. -/
theorem ZSpan.mul_left {S T : α → Prop} (a : α) (hST : ∀ x, S x → ZSpan T (a * x))
    {x : α} (h : ZSpan S x) : ZSpan T (a * x) := by
  induction h with
  | zero => have : a * (0 : α) = 0 := by grind
            rw [this]; exact .zero
  | gen hx => exact hST _ hx
  | add _ _ ih1 ih2 =>
      rename_i y z _ _
      have : a * (y + z) = a * y + a * z := Lean.Grind.Semiring.left_distrib _ _ _
      rw [this]; exact .add ih1 ih2
  | neg _ ih =>
      rename_i y _
      have : a * -y = -(a * y) := by grind
      rw [this]; exact .neg ih

/-- Generators of the Gauss-period window of length `ℓ`. -/
def GPGen (γ γ' : α) (ℓ : Nat) (x : α) : Prop := ∃ k, 1 ≤ k ∧ k ≤ ℓ ∧ x = gp γ γ' k

/-- Generators of the power window `c₁, …, c₁^ℓ`. -/
def PowGen (γ γ' : α) (ℓ : Nat) (x : α) : Prop :=
  ∃ k, 1 ≤ k ∧ k ≤ ℓ ∧ x = gp γ γ' 1 ^ k

section Char2

variable (γ γ' : α) (h : γ * γ' = 1) (two : ∀ x : α, x + x = 0)
include h two

omit h in
theorem gp_zero : gp γ γ' 0 = 0 := by
  unfold gp; simp only [Lean.Grind.Semiring.pow_zero]; exact two 1

/-- In characteristic 2: `c (k+1) = c₁ · c k + c (k−1)`. -/
theorem gp_succ (k : Nat) (hk : 1 ≤ k) :
    gp γ γ' (k + 1) = gp γ γ' 1 * gp γ γ' k + gp γ γ' (k - 1) := by
  have := dickson_mul γ γ' h (i := k) (j := 1) hk
  have e : gp γ γ' 1 * gp γ γ' k = gp γ γ' k * gp γ γ' 1 :=
    Lean.Grind.CommSemiring.mul_comm _ _
  rw [e, this]
  have := two (gp γ γ' (k - 1))
  grind

/-- Every `c k` (`1 ≤ k ≤ ℓ`) lies in the span of `c₁, …, c₁^ℓ`. -/
theorem gp_mem_pow_window : ∀ k, 1 ≤ k → ∀ ℓ, k ≤ ℓ →
    ZSpan (PowGen γ γ' ℓ) (gp γ γ' k) := by
  intro k
  induction k using Nat.strongRecOn with
  | _ k ih =>
    intro hk ℓ hkl
    match k, hk, ih with
    | 1, _, _ => exact .gen ⟨1, Nat.le_refl _, hkl, (Lean.Grind.Semiring.pow_one _).symm⟩
    | k + 2, _, ih =>
      rw [gp_succ γ γ' h two (k + 1) (by omega)]
      refine .add ?_ ?_
      · -- c₁ · c (k+1): c (k+1) ∈ span(c₁..c₁^(k+1)), times c₁ ∈ span(c₁..c₁^(k+2))
        have hmem := ih (k + 1) (by omega) (by omega) (k + 1) (Nat.le_refl _)
        refine ZSpan.mul_left _ ?_ hmem
        rintro x ⟨j, hj1, hj2, rfl⟩
        refine .gen ⟨j + 1, by omega, by omega, ?_⟩
        rw [Lean.Grind.Semiring.pow_succ]; exact Lean.Grind.CommSemiring.mul_comm _ _
      · show ZSpan _ (gp γ γ' (k + 1 - 1))
        by_cases hk0 : k = 0
        · subst hk0; rw [show (0 + 1 - 1) = 0 from rfl, gp_zero γ γ' two]; exact .zero
        · exact ih (k + 1 - 1) (by omega) (by omega) ℓ (by omega)

/-- Every `c₁^k` (`1 ≤ k ≤ ℓ`) lies in the span of `c 1, …, c ℓ`. -/
theorem pow_mem_gp_window : ∀ k, 1 ≤ k → ∀ ℓ, k ≤ ℓ →
    ZSpan (GPGen γ γ' ℓ) (gp γ γ' 1 ^ k) := by
  intro k hk
  induction k with
  | zero => omega
  | succ k ih =>
    intro ℓ hkl
    by_cases hk0 : k = 0
    · subst hk0; exact .gen ⟨1, Nat.le_refl _, hkl, Lean.Grind.Semiring.pow_one _⟩
    · have hmem := ih (by omega) k (by omega)
      have e : gp γ γ' 1 ^ (k + 1) = gp γ γ' 1 * gp γ γ' 1 ^ k := by
        rw [Lean.Grind.Semiring.pow_succ]; exact Lean.Grind.CommSemiring.mul_comm _ _
      rw [e]
      refine ZSpan.mul_left _ ?_ hmem
      rintro x ⟨j, hj1, hj2, rfl⟩
      -- c₁ · c j = c (j+1) + c (j−1)
      have e2 : gp γ γ' 1 * gp γ γ' j = gp γ γ' (j + 1) + gp γ γ' (j - 1) := by
        rw [Lean.Grind.CommSemiring.mul_comm]; exact dickson_mul γ γ' h hj1
      rw [e2]
      refine .add (.gen ⟨j + 1, by omega, by omega, rfl⟩) ?_
      by_cases hj : j = 1
      · subst hj; rw [show 1 - 1 = 0 from rfl, gp_zero γ γ' two]; exact .zero
      · exact .gen ⟨j - 1, by omega, by omega, rfl⟩

/-- **K1, set-level statement.** In characteristic 2 the Gauss-period window
`span(c 1, …, c ℓ)` equals the power window `span(c₁, c₁², …, c₁^ℓ)`. -/
theorem window_eq (ℓ : Nat) (x : α) :
    ZSpan (GPGen γ γ' ℓ) x ↔ ZSpan (PowGen γ γ' ℓ) x := by
  constructor
  · intro hx
    induction hx with
    | zero => exact .zero
    | gen hg => obtain ⟨k, h1, h2, rfl⟩ := hg; exact gp_mem_pow_window γ γ' h two k h1 ℓ h2
    | add _ _ ih1 ih2 => exact .add ih1 ih2
    | neg _ ih => exact .neg ih
  · intro hx
    induction hx with
    | zero => exact .zero
    | gen hg => obtain ⟨k, h1, h2, rfl⟩ := hg; exact pow_mem_gp_window γ γ' h two k h1 ℓ h2
    | add _ _ ih1 ih2 => exact .add ih1 ih2
    | neg _ ih => exact .neg ih

end Char2

end FactorBaseDesigns
