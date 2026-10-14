/-!
# A point moved by `τ³` has a `τ`-orbit of nine points

In the cover census, `τ` is an `F_q`-structure on `P¹` over `K = F_{q³}` and
`τ³` is the Frobenius of `K`. A root `θ` of an irreducible cubic over `K`
satisfies `τ⁹ θ = θ` and `τ³ θ ≠ θ`. This file proves that then no shorter
iterate returns to `θ`: `τ^k θ ≠ θ` for `1 ≤ k ≤ 8`. With `τ` injective the
nine iterates are therefore distinct, which gives `r = 1 + 9 = 10` branch
points and genus `2·(10 − 4) + 1 = 13` for a curve of odd order.

No Mathlib. Not formalized: the identification of `τ`, the branch set and
the genus formula with the geometric objects.
-/

namespace OrbitNine

variable {α : Type} (f : α → α)

def iter : Nat → α → α
  | 0, x => x
  | n + 1, x => iter n (f x)

theorem iter_add (a b : Nat) (x : α) : iter f (a + b) x = iter f b (iter f a x) := by
  induction a generalizing x with
  | zero => simp [iter]
  | succ a ih => rw [Nat.succ_add]; exact ih (f x)

/-- `x` returns after `a` steps. -/
def Ret (a : Nat) (x : α) : Prop := iter f a x = x

theorem ret_add {a b : Nat} {x : α} (ha : Ret f a x) (hb : Ret f b x) : Ret f (a + b) x := by
  unfold Ret at *; rw [iter_add, ha, hb]

theorem ret_sub {a b : Nat} {x : α} (ha : Ret f a x) (hab : Ret f (a + b) x) : Ret f b x := by
  unfold Ret at *; rw [iter_add, ha] at hab; exact hab

/-- If `τ⁹ x = x` and `τ³ x ≠ x`, then `τ^k x ≠ x` for every `1 ≤ k ≤ 8`. -/
theorem no_short_return {x : α} (h9 : Ret f 9 x) (h3 : ¬ Ret f 3 x) :
    ∀ k, 1 ≤ k → k ≤ 8 → ¬ Ret f k x := by
  have one : Ret f 1 x → False := fun h => h3 (ret_add f (ret_add f h h) h)
  intro k hk1 hk8 hk
  match k, hk1, hk8 with
  | 1, _, _ => exact one hk
  | 2, _, _ =>
    have h4 : Ret f (2 + 2) x := ret_add f hk hk
    have h8 : Ret f (2 + 2 + (2 + 2)) x := ret_add f h4 h4
    exact one (ret_sub f (a := 8) (b := 1) h8 h9)
  | 3, _, _ => exact h3 hk
  | 4, _, _ =>
    have h8 : Ret f (4 + 4) x := ret_add f hk hk
    exact one (ret_sub f (a := 8) (b := 1) h8 h9)
  | 5, _, _ =>
    have h10 : Ret f (5 + 5) x := ret_add f hk hk
    exact one (ret_sub f (a := 9) (b := 1) h9 h10)
  | 6, _, _ => exact h3 (ret_sub f (a := 6) (b := 3) hk h9)
  | 7, _, _ =>
    have h14 : Ret f (7 + 7) x := ret_add f hk hk
    have h5 : Ret f 5 x := ret_sub f (a := 9) (b := 5) h9 h14
    have h10 : Ret f (5 + 5) x := ret_add f h5 h5
    exact one (ret_sub f (a := 9) (b := 1) h9 h10)
  | 8, _, _ => exact one (ret_sub f (a := 8) (b := 1) hk h9)

end OrbitNine

#print axioms OrbitNine.no_short_return
