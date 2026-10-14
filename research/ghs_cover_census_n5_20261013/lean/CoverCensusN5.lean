/-!
# Lemmas behind the `n = 5` cover census (PROTOCOL.md, "Derived statements")

No Mathlib. What is formalized here, and what is not:

* `Orbit.ret_gcd`, `Orbit.period_three_or_fifteen`: if `τ¹⁵ θ = θ` and
  `τ⁵ θ ≠ θ` then every return of `θ` before 15 steps is a multiple of 3 and
  `τ θ ≠ θ` (orbit size 3 or 15). Used in Proposition O5 with `τ⁵` the
  Frobenius of `K = F_{q⁵}` and `θ` a root of an irreducible cubic over `K`.
  `Orbit.no_short_return_nine` re-derives the merged `OrbitNine` lemma from
  `ret_gcd`.
* `Rot5.span_sizes`: in `F_2^5` with the cyclic shift, the span of the five
  rotations of any vector has 1, 2, 16 or 32 elements, 16 exactly when the
  five rotations sum to 0 and the vector is not shift-invariant (Lemma M for
  a single free 5-orbit). Checked by `decide` over all 32 vectors.
* `Parity.even_one_sub_neg_pow`: `1 − (−q)^k` is even for odd `q`
  (the exponent in Proposition T).
* `Spectrum.n5_gap`: with `m = 4, r ∈ {5, 10}` or `m = 5, 6 ≤ r ≤ 20`, the
  genus `2^(m−2)(r−4)+1` is at least 5 and never in `6..16`.

Not formalized: Lang's theorem and the reduction to `σ`; Lemma M for an
arbitrary `σ`-set (the cyclic-module argument with the factorisation of
`x⁵ − 1`); Lemma E's configuration analysis; the simple transitivity in
Proposition O5; the cross-ratio invariant and the 2-descent step of
Proposition T; the genus formula and Lemma H.
-/

namespace Orbit

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

theorem ret_zero (x : α) : Ret f 0 x := rfl

theorem ret_add {a b : Nat} {x : α} (ha : Ret f a x) (hb : Ret f b x) : Ret f (a + b) x := by
  unfold Ret at *; rw [iter_add, ha, hb]

theorem ret_sub {a b : Nat} {x : α} (ha : Ret f a x) (hab : Ret f (a + b) x) : Ret f b x := by
  unfold Ret at *; rw [iter_add, ha] at hab; exact hab

theorem ret_mul {a : Nat} {x : α} (ha : Ret f a x) : ∀ k, Ret f (k * a) x
  | 0 => by simp [Ret, iter]
  | k + 1 => by rw [Nat.succ_mul]; exact ret_add f (ret_mul ha k) ha

theorem ret_mod {a b : Nat} {x : α} (ha : Ret f a x) (hb : Ret f b x) : Ret f (b % a) x := by
  have h := Nat.div_add_mod b a
  have hm : Ret f (a * (b / a)) x := by rw [Nat.mul_comm]; exact ret_mul f ha _
  exact ret_sub f hm (by rw [h]; exact hb)

/-- Returns are closed under `gcd`. -/
theorem ret_gcd {x : α} : ∀ a b, Ret f a x → Ret f b x → Ret f (Nat.gcd a b) x := by
  intro a
  induction a using Nat.strongRecOn with
  | _ a ih =>
    intro b ha hb
    cases a with
    | zero => rw [Nat.gcd_zero_left]; exact hb
    | succ a =>
      rw [Nat.gcd_rec]
      exact ih (b % (a + 1)) (Nat.mod_lt _ (Nat.succ_pos a)) (a + 1) (ret_mod f ha hb) ha

/-- If `τ¹⁵ θ = θ` and `τ⁵ θ ≠ θ`, every return before 15 steps is a multiple of 3,
and `θ` is not fixed: its orbit has 3 or 15 points. -/
theorem period_three_or_fifteen {x : α} (h15 : Ret f 15 x) (h5 : ¬ Ret f 5 x) :
    ¬ Ret f 1 x ∧ ∀ k, 0 < k → k < 15 → Ret f k x → k % 3 = 0 := by
  have one : ¬ Ret f 1 x := fun h => h5 (ret_mul f h 5)
  refine ⟨one, ?_⟩
  intro k hk0 hk15 hk
  have hg := ret_gcd f k 15 hk h15
  match k, hk0, hk15 with
  | 1, _, _ => exact absurd hk one
  | 2, _, _ => exact absurd hg one
  | 3, _, _ => rfl
  | 4, _, _ => exact absurd hg one
  | 5, _, _ => exact absurd hk h5
  | 6, _, _ => rfl
  | 7, _, _ => exact absurd hg one
  | 8, _, _ => exact absurd hg one
  | 9, _, _ => rfl
  | 10, _, _ => exact absurd hg h5
  | 11, _, _ => exact absurd hg one
  | 12, _, _ => rfl
  | 13, _, _ => exact absurd hg one
  | 14, _, _ => exact absurd hg one

/-- The merged census's `OrbitNine.no_short_return`, re-derived from `ret_gcd`. -/
theorem no_short_return_nine {x : α} (h9 : Ret f 9 x) (h3 : ¬ Ret f 3 x) :
    ∀ k, 1 ≤ k → k ≤ 8 → ¬ Ret f k x := by
  have one : ¬ Ret f 1 x := fun h => h3 (ret_mul f h 3)
  intro k hk1 hk8 hk
  have hg := ret_gcd f k 9 hk h9
  match k, hk1, hk8 with
  | 1, _, _ => exact one hk
  | 2, _, _ => exact one hg
  | 3, _, _ => exact h3 hk
  | 4, _, _ => exact one hg
  | 5, _, _ => exact one hg
  | 6, _, _ => exact h3 hg
  | 7, _, _ => exact one hg
  | 8, _, _ => exact one hg

end Orbit

namespace Rot5

/-- cyclic shift on 5-bit vectors -/
def rot (v : Nat) : Nat := ((v <<< 1) ||| (v >>> 4)) &&& 31

def rots (v : Nat) : List Nat := [v, rot v, rot (rot v), rot (rot (rot v)), rot (rot (rot (rot v)))]

/-- xor of the members of `l` selected by the bits of `s` -/
def comb : List Nat → Nat → Nat
  | [], _ => 0
  | a :: l, s => (if s % 2 = 1 then a else 0) ^^^ comb l (s / 2)

/-- number of distinct vectors `w < 32` in the span of the rotations of `v` -/
def spanSize (v : Nat) : Nat :=
  ((List.range 32).filter (fun w => (List.range 32).any (fun s => comb (rots v) s == w))).length

def sumRots (v : Nat) : Nat := comb (rots v) 31

theorem span_sizes : ∀ v, v < 32 →
    (spanSize v = 1 ∨ spanSize v = 2 ∨ spanSize v = 16 ∨ spanSize v = 32) ∧
    (spanSize v = 16 ↔ (sumRots v = 0 ∧ rot v ≠ v)) ∧
    (spanSize v = 2 ↔ (v ≠ 0 ∧ rot v = v)) := by
  decide

end Rot5

namespace Parity

theorem neg_pow_odd (q : Int) (hq : q % 2 = 1) : ∀ k : Nat, ((-q) ^ k) % 2 = 1
  | 0 => by rw [Int.pow_zero]; decide
  | k + 1 => by
    have ih := neg_pow_odd q hq k
    have hn : (-q) % 2 = 1 := by omega
    rw [Int.pow_succ, Int.mul_emod, ih, hn]; decide

/-- `1 − (−q)^k` is even when `q` is odd. -/
theorem even_one_sub_neg_pow (q : Int) (hq : q % 2 = 1) (k : Nat) : (1 - (-q) ^ k) % 2 = 0 := by
  have h := neg_pow_odd q hq k
  omega

end Parity

namespace Spectrum

def genus (m r : Nat) : Nat := 2 ^ (m - 2) * (r - 4) + 1

/-- The `n = 5` genus spectrum allowed by Lemmas M and E: at least 5, nothing in 6..16. -/
theorem n5_gap : ∀ r, r ≤ 20 →
    ((r = 5 ∨ r = 10) → 5 ≤ genus 4 r ∧ ¬ (6 ≤ genus 4 r ∧ genus 4 r ≤ 16)) ∧
    (6 ≤ r → 17 ≤ genus 5 r) := by
  decide

end Spectrum

#print axioms Orbit.period_three_or_fifteen
#print axioms Orbit.no_short_return_nine
#print axioms Rot5.span_sizes
#print axioms Parity.even_one_sub_neg_pow
#print axioms Spectrum.n5_gap
