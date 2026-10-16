/-!
# Mean decomposition yield is fixed by the domain size (K2, K5)

A decomposition *pattern* is a map `φ : τ → G` from a finite domain `D` of
candidate tuples (base points, Frobenius-rotated slots, partial sums, …) to
the group `G`. The yield at a target `g` is `count D φ g`, the number of
tuples landing on `g`.

* `sum_count`: summed over all targets, the yield is exactly `|D|`, whatever
  `φ` is. So the **mean** yield per target is `|D| / |G|` for every pattern.
* **K2** (`patterns_neutral`): two patterns on domains of equal size have
  equal total yield. Solving `k` rotated-slot patterns costs `k` solves and
  gives `k` times the mean yield of one, so per solve it is neutral **in the
  mean over targets**. This theorem does not speak to a single fixed target,
  nor to a merged solve whose cost differs from `k` separate ones. Those are
  the open parts of register entry K2.
* **K5** (`restrict_le`): restricting the domain (for example, requiring the
  intermediate points of a 3-term chain to lie in a subspace `W`) can only
  lower the yield at every target, and so it lowers the total too.
-/

namespace FactorBaseDesigns

variable {τ G : Type} [DecidableEq G]

/-- Yield of pattern `φ` on domain `D` at target `g`. -/
def count (D : List τ) (φ : τ → G) (g : G) : Nat :=
  (D.filter (fun x => decide (φ x = g))).length

omit [DecidableEq G] in
theorem sum_map_zero (Gs : List G) : (Gs.map (fun _ => (0 : Nat))).sum = 0 := by
  induction Gs with
  | nil => rfl
  | cons g Gs ih => simp [List.map_cons, List.sum_cons, ih]

/-- Sum over a duplicate-free target list of the indicator of `a`. -/
theorem sum_indicator {a : G} : ∀ {Gs : List G}, Gs.Nodup → a ∈ Gs →
    (Gs.map (fun g => if a = g then 1 else 0)).sum = 1
  | [], _, h => absurd h (List.not_mem_nil)
  | g :: Gs, hnd, hmem => by
    rw [List.nodup_cons] at hnd
    simp only [List.map_cons, List.sum_cons]
    by_cases hag : a = g
    · subst hag
      have : (Gs.map (fun g => if a = g then 1 else 0)).sum = 0 := by
        apply Nat.eq_zero_of_le_zero
        have hz : ∀ x ∈ Gs, (if a = x then 1 else 0) = 0 := by
          intro x hx
          have : a ≠ x := fun e => hnd.1 (e ▸ hx)
          simp [this]
        rw [List.map_congr_left hz]; exact Nat.le_of_eq (sum_map_zero Gs)
      simp [this]
    · have hmem' : a ∈ Gs := by
        rcases List.mem_cons.mp hmem with h | h
        · exact absurd h hag
        · exact h
      simp [hag, sum_indicator hnd.2 hmem']

theorem count_cons (x : τ) (D : List τ) (φ : τ → G) (g : G) :
    count (x :: D) φ g = count D φ g + (if φ x = g then 1 else 0) := by
  unfold count
  by_cases h : φ x = g <;> simp [h]

omit [DecidableEq G] in
theorem sum_add_map (Gs : List G) (f h : G → Nat) :
    (Gs.map (fun g => f g + h g)).sum = (Gs.map f).sum + (Gs.map h).sum := by
  induction Gs with
  | nil => simp
  | cons g Gs ih => simp [List.map_cons, List.sum_cons, ih]; omega

/-- **Total yield equals the domain size.** -/
theorem sum_count (D : List τ) (φ : τ → G) {Gs : List G} (hnd : Gs.Nodup)
    (hcov : ∀ x ∈ D, φ x ∈ Gs) : (Gs.map (count D φ)).sum = D.length := by
  induction D with
  | nil =>
    show (Gs.map (fun _ => 0)).sum = 0
    exact sum_map_zero Gs
  | cons x D ih =>
    have hcov' : ∀ y ∈ D, φ y ∈ Gs := fun y hy => hcov y (List.mem_cons_of_mem _ hy)
    have hx : φ x ∈ Gs := hcov x (List.mem_cons_self)
    have : (Gs.map (count (x :: D) φ)) =
        Gs.map (fun g => count D φ g + (if φ x = g then 1 else 0)) := by
      apply List.map_congr_left; intro g _; exact count_cons x D φ g
    rw [this, sum_add_map, ih hcov', sum_indicator hnd hx]
    simp

/-- **K2.** Two decomposition patterns over domains of the same size have the
same total yield over all targets, hence the same mean yield per target. -/
theorem patterns_neutral {σ : Type} (D : List τ) (E : List σ) (φ : τ → G) (ψ : σ → G)
    {Gs : List G} (hnd : Gs.Nodup) (hφ : ∀ x ∈ D, φ x ∈ Gs) (hψ : ∀ y ∈ E, ψ y ∈ Gs)
    (hlen : D.length = E.length) :
    (Gs.map (count D φ)).sum = (Gs.map (count E ψ)).sum := by
  rw [sum_count D φ hnd hφ, sum_count E ψ hnd hψ, hlen]

/-- **K5.** Restricting the domain by any predicate (e.g. "the chain's
intermediate points lie in `W`") never increases the yield at any target. -/
theorem restrict_le (D : List τ) (φ : τ → G) (keep : τ → Bool) (g : G) :
    count (D.filter keep) φ g ≤ count D φ g := by
  induction D with
  | nil => simp [count]
  | cons x D ih =>
    by_cases hk : keep x = true
    · simp only [List.filter_cons, hk, ite_true]
      rw [count_cons, count_cons]; omega
    · simp only [List.filter_cons, hk]
      rw [count_cons]; simp at ih ⊢; omega

end FactorBaseDesigns
