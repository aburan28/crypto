/-!
# Multiplication by a unit on `ℤ/N` has fibres of size at most one

This is the arithmetic core of two no-go statements in the design register
(README, K7 and K8).

* **K8.** Doubling or halving fibres on the target curve are trivial. A group
  of odd order `r` that is cyclic is, as an abstract group, `(ℤ/r, +)`, and
  `[2^t]` is multiplication by `2^t`. Since `gcd(2^t, r) = 1`, every fibre
  `{P : 2^t·P = C}` has at most one element.
* **K7.** Power-map towers in `F_{2^n}^*` are trivial below the smallest
  prime factor of `2^n − 1`. The multiplicative group of a finite field is
  cyclic, of order `N = 2^n − 1`. Through a discrete logarithm,
  `x ↦ x^d` becomes multiplication by `d` on `ℤ/N`, so every fibre of a
  degree-`d` power step with `gcd(d, N) = 1` has at most one element.

**Bridge (not formalized here).** The identification of a cyclic group of
order `N` with `(ℤ/N, +)` is the textbook fact that both K7 and K8 rely on.
The same holds for the cyclicity of a group of prime order and of a finite
field's unit group (Mathlib: `isCyclic_of_prime_card`,
`instIsCyclicUnitsOfFinite`). Mathlib is not a dependency of this package (see
README, "Toolchain"), so the theorems below are stated on the model
`{a : Nat // a < N}` with multiplication mod `N`. The bridge is the one
unformalized step, and the proof map names it.
-/

namespace FactorBaseDesigns

/-- Multiplication by `k` modulo `N` is injective on `[0, N)` when
`gcd(k, N) = 1`. -/
theorem mul_mod_injective {N k a b : Nat} (h : Nat.Coprime k N)
    (ha : a < N) (hb : b < N) (e : k * a % N = k * b % N) : a = b := by
  have key : ∀ {x y : Nat}, x ≤ y → y < N → k * x % N = k * y % N → x = y := by
    intro x y hxy hy hxy'
    have h0 : (k * y - k * x) % N = 0 := Nat.sub_mod_eq_zero_of_mod_eq hxy'.symm
    have hd : N ∣ k * (y - x) := by
      rw [Nat.mul_sub]; exact Nat.dvd_of_mod_eq_zero h0
    have hd' : N ∣ y - x := Nat.Coprime.dvd_of_dvd_mul_left h.symm hd
    have hlt : y - x < N := by omega
    have : y - x = 0 := Nat.eq_zero_of_dvd_of_lt hd' hlt
    omega
  rcases Nat.le_total a b with hab | hab
  · exact key hab hb e
  · exact (key hab ha e.symm).symm

/-- Fibres of `a ↦ k·a mod N` on `[0, N)` have at most one element: any two
points of the fibre over `c` are equal. -/
theorem fibre_subsingleton {N k c : Nat} (h : Nat.Coprime k N) :
    ∀ a b, a < N → b < N → k * a % N = c → k * b % N = c → a = b :=
  fun _ _ ha hb e1 e2 => mul_mod_injective h ha hb (e1.trans e2.symm)

/-- `gcd(2^t, r) = 1` for odd `r`. -/
theorem coprime_two_pow_of_odd {r : Nat} (hr : r % 2 = 1) (t : Nat) :
    Nat.Coprime (2 ^ t) r := by
  have h2 : Nat.Coprime 2 r := by
    rw [Nat.coprime_iff_gcd_eq_one, Nat.gcd_rec, hr]
    rfl
  exact Nat.Coprime.pow_left t h2

/-- **K8 (model form).** In `(ℤ/r, +)` with `r` odd, the fibre of `[2^t]`
over any `c` has at most one point. A depth-`t` doubling tower on a cyclic
group of odd order (in particular of odd prime order) therefore cuts out at
most one point, so it is not a factor base. -/
theorem doubling_fibre_subsingleton {r : Nat} (hr : r % 2 = 1) (t c : Nat) :
    ∀ a b, a < r → b < r → 2 ^ t * a % r = c → 2 ^ t * b % r = c → a = b :=
  fibre_subsingleton (coprime_two_pow_of_odd hr t)

end FactorBaseDesigns
