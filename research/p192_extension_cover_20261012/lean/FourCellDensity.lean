import Mathlib.Tactic.NormNum

/-!
The same character-sum bounds control all four two-fiber sign cells.
The concrete genus-seven and genus-thirty-one covers, Hasse--Weil bounds,
and finite-field character identity are supplied by the P-192 density
certificate and the written geometric proof.
-/

namespace P192ExtensionCover

/-- If doubling permutes the elliptic group, the second one-character sum
    equals the first. Every choice of two signs therefore has the same
    worst-case lower bound, allowing one exceptional identity point. -/
theorem four_cell_character_lower_bound
    (A r S₁ S₁₂ B₁ B₁₂ δ ε η : ℤ)
    (hε : ε = 1 ∨ ε = -1) (hη : η = 1 ∨ η = -1)
    (hδ : δ ≤ 1)
    (hcount : 4 * (A + δ) = r + (ε + η) * S₁ + ε * η * S₁₂)
    (hS₁ : -B₁ ≤ S₁ ∧ S₁ ≤ B₁)
    (hS₁₂ : -B₁₂ ≤ S₁₂ ∧ S₁₂ ≤ B₁₂) :
    r - 2 * B₁ - B₁₂ - 4 ≤ 4 * A := by
  rcases hε with hε | hε <;> rcases hη with hη | hη
  all_goals
    subst ε
    subst η
    omega

/-- Exact P-192 integer consequence of the four-cell character estimate.
    The hypotheses are the concrete character-sum identity and its two
    Hasse--Weil bounds, with the possible identity exception explicit. -/
theorem p192_four_cell_nonzero_scalar_lower_bound
    (A S₁ S₁₂ δ ε η : ℤ)
    (hε : ε = 1 ∨ ε = -1) (hη : η = 1 ∨ η = -1)
    (hδ : δ ≤ 1)
    (hcount : 4 * (A + δ) =
      6277101735386680763835789423176059013767194773182842284081 +
        (ε + η) * S₁ + ε * η * S₁₂)
    (hS₁ : -(14 * 79228162514264337593543950336 +
        31607402316713927207482677199) ≤ S₁ ∧
        S₁ ≤ 14 * 79228162514264337593543950336 +
          31607402316713927207482677199)
    (hS₁₂ : -(62 * 79228162514264337593543950336 +
        31607402316713927207482677199) ≤ S₁₂ ∧
        S₁₂ ≤ 62 * 79228162514264337593543950336 +
          31607402316713927207482677199) :
    1569275433846670190958947353987675545133315652035359680560 ≤ A := by
  have h := four_cell_character_lower_bound A
    6277101735386680763835789423176059013767194773182842284081
    S₁ S₁₂
    (14 * 79228162514264337593543950336 + 31607402316713927207482677199)
    (62 * 79228162514264337593543950336 + 31607402316713927207482677199)
    δ ε η hε hη hδ hcount hS₁ hS₁₂
  omega

/-- Every cell of the exact 205-class support has fewer than 364,825
    targets, giving the stated ideal-uniform conditional hit inequality
    whenever its accepted nonzero scalar count meets the theorem above. -/
theorem p192_205_class_four_cell_hit_integer_certificate :
    364824 * (2 ^ 171 : ℕ) <
      1569275433846670190958947353987675545133315652035359680560 := by
  norm_num

end P192ExtensionCover

#print axioms P192ExtensionCover.four_cell_character_lower_bound
#print axioms P192ExtensionCover.p192_four_cell_nonzero_scalar_lower_bound
#print axioms P192ExtensionCover.p192_205_class_four_cell_hit_integer_certificate
