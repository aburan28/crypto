import Mathlib.Algebra.Group.Basic
import Mathlib.Algebra.Group.Equiv.Defs
import Mathlib.Tactic.Group

set_option autoImplicit false

namespace WeakCurveDatabase

/-- The multiplicative form of the forward cubic Hilbert-90 construction.
The commutative group may be the nonzero elements of a field, and σ its
cyclic cubic Galois generator. The hypothesis is exactly σ³(a)=a. -/
theorem cubic_cocycle_norm_one {G : Type*} [CommGroup G]
    (σ : G ≃* G) (a : G) (hσ : σ (σ (σ a)) = a) :
    (σ a * a⁻¹) *
      σ (σ a * a⁻¹) *
      σ (σ (σ a * a⁻¹)) = 1 := by
  simp only [map_mul, map_inv]
  rw [hσ]
  simp [mul_comm, mul_left_comm, mul_assoc]

/-- In a finite-field multiplicative group, Frobenius acts by the q-th
power. Under that exact hypothesis the cocycle lies in the norm-one torus
of order dividing q²+q+1. -/
theorem cubic_cocycle_power {G : Type*} [CommGroup G]
    (σ : G ≃* G) (q : ℕ) (hq : ∀ x : G, σ x = x^q)
    (a : G) (hσ : σ (σ (σ a)) = a) :
    (σ a * a⁻¹)^(q*q+q+1) = 1 := by
  let lam := σ a * a⁻¹
  have hnorm : lam * σ lam * σ (σ lam) = 1 :=
    cubic_cocycle_norm_one σ a hσ
  have hsecond : σ lam = lam^q := hq lam
  have hthird : σ (σ lam) = lam^(q*q) := by
    rw [hq (σ lam), hsecond, pow_mul]
  rw [hthird, hsecond] at hnorm
  calc
    lam^(q*q+q+1) = lam^(1+q+q*q) := by congr 1; omega
    _ = lam * lam^q * lam^(q*q) := by simp [pow_add, mul_assoc]
    _ = 1 := hnorm

/-- If an element's annihilating exponent is three modulo four, a
fourth root is given explicitly by raising it to `(n+1)/4`. -/
theorem fourth_power_of_order_three_mod_four {G : Type*} [Monoid G]
    (x : G) (n : ℕ) (hn : n % 4 = 3) (hx : x^n = 1) :
    ∃ y : G, y^4 = x := by
  let k := (n+1)/4
  have hk : k*4 = n+1 := by omega
  refine ⟨x^k, ?_⟩
  rw [← pow_mul, hk, pow_succ, hx, one_mul]

/-- An odd base prime has square congruent to one modulo four. -/
theorem odd_square_mod_four (p : ℕ) (hp : p % 2 = 1) :
    (p*p) % 4 = 1 := by
  have hp4 : p % 4 = 1 ∨ p % 4 = 3 := by omega
  rcases hp4 with h | h <;> simp [Nat.mul_mod, h]

/-- For `q ≡ 1 (mod 4)`, every element annihilated by the cubic torus
order `q²+q+1` has an explicit fourth root. In the target family
`q=p²` with odd `p`, so this congruence holds. -/
theorem cubic_torus_fourth_power {G : Type*} [Monoid G]
    (x : G) (q : ℕ) (hq : q % 4 = 1)
    (hx : x^(q*q+q+1) = 1) : ∃ y : G, y^4 = x := by
  have hn : (q*q+q+1) % 4 = 3 := by
    simp [Nat.add_mod, Nat.mul_mod, hq]
  exact fourth_power_of_order_three_mod_four x (q*q+q+1) hn hx

/-- The fourth-power implication with the exact base-field square `q=p²`.
No finite-field object is assumed: only the torus exponent equation. -/
theorem cubic_torus_fourth_power_prime_square {G : Type*} [Monoid G]
    (x : G) (p : ℕ) (hp : p % 2 = 1)
    (hx : x^((p*p)*(p*p)+(p*p)+1) = 1) : ∃ y : G, y^4 = x := by
  exact cubic_torus_fourth_power x (p*p) (odd_square_mod_four p hp) hx

end WeakCurveDatabase
