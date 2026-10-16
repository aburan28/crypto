import FactorBaseDesigns.CyclicMul

/-!
# Kernel-checked number facts used by the design register

Every `decide +kernel` here is evaluated by the Lean kernel (with its
GMP-accelerated `Nat` operations), not by native code. There is no
`native_decide` in this package.

* **K7.** For each standard binary degree `n`, every `d` with
  `2 ≤ d < B_n` is coprime to `2^n − 1` (`k7_coprime_*`), via one gcd
  against the product `∏_{2 ≤ d < B_n} d`. With `CyclicMul.fibre_subsingleton`,
  every power-map tower step of degree below `B_n` has fibres of size ≤ 1 in
  the cyclic model of `F_{2^n}^*` (`k7_power_step_trivial`). Where the bound
  is the smallest prime factor, its divisibility is checked too.
* **Prime structure** used by designs #5 and #6: 2-adic valuations and small
  odd parts of `p ± 1` for the NIST primes and secp256k1.
* **Type-II normal bases** (register entry K1, experiment E-K1): orders of 2
  modulo `2n + 1`.
* **Bits per unit step degree** (#5): `q^3 < 3^q` for every integer
  `q ≥ 2` with `q ≠ 3`. Equivalently, `log₂ q / q < log₂ 3 / 3`, so among
  integer step degrees, degree 3 uniquely maximizes the factor-base bits
  bought per unit of equation degree.
-/

set_option exponentiation.threshold 1024

namespace FactorBaseDesigns

/-- `pr lo len = lo · (lo+1) ⋯ (lo+len−1)`. -/
def pr (lo : Nat) : Nat → Nat
  | 0 => 1
  | len + 1 => pr lo len * (lo + len)

/-- The same product by binary splitting (fast in the kernel), structural in
the fuel. -/
def prodR : Nat → Nat → Nat → Nat
  | 0, lo, len => pr lo len
  | f + 1, lo, len =>
    if len ≤ 64 then pr lo len else prodR f lo (len / 2) * prodR f (lo + len / 2) (len - len / 2)

theorem pr_add (lo a : Nat) : ∀ b, pr lo (a + b) = pr lo a * pr (lo + a) b
  | 0 => by simp [pr]
  | b + 1 => by
    rw [← Nat.add_assoc, pr, pr_add lo a b, pr, Nat.mul_assoc, Nat.add_assoc]

theorem prodR_eq : ∀ f lo len, prodR f lo len = pr lo len
  | 0, _, _ => rfl
  | f + 1, lo, len => by
    unfold prodR
    split
    · rfl
    · rw [prodR_eq f, prodR_eq f, ← pr_add]
      congr 1; omega

theorem dvd_pr {lo d : Nat} : ∀ {len}, lo ≤ d → d < lo + len → d ∣ pr lo len
  | 0, _, h => by omega
  | len + 1, h1, h2 => by
    unfold pr
    by_cases hd : d = lo + len
    · subst hd; exact Nat.dvd_mul_left _ _
    · exact Nat.dvd_trans (dvd_pr h1 (by omega)) (Nat.dvd_mul_right _ _)

/-- A gcd certificate on the product gives coprimality of every factor. -/
theorem coprime_of_prod {N lo len f : Nat} (h : Nat.gcd (prodR f lo len) N = 1) :
    ∀ d, lo ≤ d → d < lo + len → Nat.Coprime d N := by
  intro d h1 h2
  rw [prodR_eq] at h
  exact Nat.Coprime.coprime_dvd_left (dvd_pr h1 h2) h

/-! ## K7 instances -/

theorem cert113 : Nat.gcd (prodR 25 2 (3391 - 2)) (2 ^ 113 - 1) = 1 := by decide +kernel
theorem cert131 : Nat.gcd (prodR 25 2 (263 - 2)) (2 ^ 131 - 1) = 1 := by decide +kernel
theorem cert163 : Nat.gcd (prodR 25 2 (150287 - 2)) (2 ^ 163 - 1) = 1 := by decide +kernel
theorem cert233 : Nat.gcd (prodR 25 2 (1399 - 2)) (2 ^ 233 - 1) = 1 := by decide +kernel
theorem cert239 : Nat.gcd (prodR 25 2 (479 - 2)) (2 ^ 239 - 1) = 1 := by decide +kernel
theorem cert283 : Nat.gcd (prodR 25 2 (9623 - 2)) (2 ^ 283 - 1) = 1 := by decide +kernel
theorem cert409 : Nat.gcd (prodR 25 2 (1000000 - 2)) (2 ^ 409 - 1) = 1 := by decide +kernel
theorem cert571 : Nat.gcd (prodR 25 2 (5711 - 2)) (2 ^ 571 - 1) = 1 := by decide +kernel

/-- The bounds are attained (smallest prime factors), except `n = 409`, where
the bound `10^6` is only a search limit. -/
theorem bounds_attained :
    (2 ^ 113 - 1) % 3391 = 0 ∧ (2 ^ 131 - 1) % 263 = 0 ∧ (2 ^ 163 - 1) % 150287 = 0 ∧
    (2 ^ 233 - 1) % 1399 = 0 ∧ (2 ^ 239 - 1) % 479 = 0 ∧ (2 ^ 283 - 1) % 9623 = 0 ∧
    (2 ^ 571 - 1) % 5711 = 0 := by decide +kernel

/-- The table of (n, B_n). -/
def k7Table : List (Nat × Nat) :=
  [(113, 3391), (131, 263), (163, 150287), (233, 1399), (239, 479), (283, 9623),
   (409, 1000000), (571, 5711)]

theorem k7_coprime : ∀ e ∈ k7Table, ∀ d, 2 ≤ d → d < e.2 → Nat.Coprime d (2 ^ e.1 - 1) := by
  intro e he d h1 h2
  simp only [k7Table, List.mem_cons, List.not_mem_nil, or_false] at he
  rcases he with rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl
  · exact coprime_of_prod cert113 d h1 (by simpa using h2)
  · exact coprime_of_prod cert131 d h1 (by simpa using h2)
  · exact coprime_of_prod cert163 d h1 (by simpa using h2)
  · exact coprime_of_prod cert233 d h1 (by simpa using h2)
  · exact coprime_of_prod cert239 d h1 (by simpa using h2)
  · exact coprime_of_prod cert283 d h1 (by simpa using h2)
  · exact coprime_of_prod cert409 d h1 (by simpa using h2)
  · exact coprime_of_prod cert571 d h1 (by simpa using h2)

/-- **K7 (model form).** For every `(n, B_n)` in the table and every power
step `x ↦ x^d` with `2 ≤ d < B_n`, the step's fibres in the cyclic model
`ℤ/(2^n − 1)` of `F_{2^n}^*` have at most one element. A tower of such steps
therefore cuts out at most one point. Bridge: `F_{2^n}^*` is cyclic of order
`2^n − 1` (not formalized here; see `CyclicMul`). -/
theorem k7_power_step_trivial : ∀ e ∈ k7Table, ∀ d, 2 ≤ d → d < e.2 →
    ∀ c a b, a < 2 ^ e.1 - 1 → b < 2 ^ e.1 - 1 →
      d * a % (2 ^ e.1 - 1) = c → d * b % (2 ^ e.1 - 1) = c → a = b :=
  fun e he d h1 h2 c => fibre_subsingleton (c := c) (k7_coprime e he d h1 h2)

/-! ## Structure of `p ± 1` for standard primes (designs #5, #6) -/

def p192 : Nat := 2 ^ 192 - 2 ^ 64 - 1
def p224 : Nat := 2 ^ 224 - 2 ^ 96 + 1
def p256 : Nat := 2 ^ 256 - 2 ^ 224 + 2 ^ 192 + 2 ^ 96 - 1
def p384 : Nat := 2 ^ 384 - 2 ^ 128 - 2 ^ 96 + 2 ^ 32 - 1
def p521 : Nat := 2 ^ 521 - 1
def pk1 : Nat := 2 ^ 256 - 2 ^ 32 - 977

/-- Product of the primes up to 31. A cofactor coprime to it has no prime
factor ≤ 31. -/
def primorial31 : Nat := 2 * 3 * 5 * 7 * 11 * 13 * 17 * 19 * 23 * 29 * 31

/-- `p = 2^a·(small part)·k` with `gcd(k, primorial31) = 1`. -/
def SmallPart (p small : Nat) : Prop :=
  p % small = 0 ∧ Nat.gcd (p / small) primorial31 = 1

instance (p s : Nat) : Decidable (SmallPart p s) := by unfold SmallPart; infer_instance

theorem p256_minus : SmallPart (p256 - 1) (2 * 3 * 5 ^ 2 * 17) := by decide +kernel
theorem p256_plus : SmallPart (p256 + 1) (2 ^ 96 * 7) := by decide +kernel
theorem p224_minus : SmallPart (p224 - 1) (2 ^ 96 * 3 * 5 * 17) := by decide +kernel
theorem p224_plus : SmallPart (p224 + 1) 2 := by decide +kernel
theorem p192_minus : SmallPart (p192 - 1) 2 := by decide +kernel
theorem p192_plus : SmallPart (p192 + 1) (2 ^ 64 * 3 * 5 * 17) := by decide +kernel
theorem p384_minus : SmallPart (p384 - 1) (2 * 19) := by decide +kernel
theorem p384_plus : SmallPart (p384 + 1) (2 ^ 32 * 3 ^ 2 * 5 * 7 * 13 * 17) := by decide +kernel
theorem pk1_minus : SmallPart (pk1 - 1) (2 * 3 * 7) := by decide +kernel
theorem pk1_plus : SmallPart (pk1 + 1) (2 ^ 4) := by decide +kernel
theorem p521_minus : SmallPart (p521 - 1) (2 * 3 * 5 ^ 2 * 11 * 17 * 31) := by decide +kernel

/-- **#6.** `p521 + 1 = 2^521`: the norm-one torus of `F_{p²}` has order
`2^521`, so it carries a 2-power Dickson tower of every depth up to 521. -/
theorem p521_plus : p521 + 1 = 2 ^ 521 := by decide +kernel

/-! ## Type-II optimal normal bases (E-K1) -/

/-- `2n+1` is prime, by trial division up to its square root. -/
def isPrimeTD (q : Nat) : Bool := 2 ≤ q && (List.range (q + 1)).all fun d => d < 2 || q < d * d || q % d != 0

/-- Multiplicative order of `2` modulo `q` equals `o`: `2^o ≡ 1` and no
proper divisor works. -/
def ordTwoIs (q o : Nat) : Bool :=
  2 ^ o % q == 1 && (List.range o).all fun e => e == 0 || o % e != 0 || 2 ^ e % q != 1

theorem typeII_131 : isPrimeTD 263 && ordTwoIs 263 131 && 263 % 4 == 3 := by decide +kernel
theorem typeII_113 : isPrimeTD 227 && ordTwoIs 227 226 := by decide +kernel
theorem typeII_233 : isPrimeTD 467 && ordTwoIs 467 466 := by decide +kernel
theorem typeII_239 : isPrimeTD 479 && ordTwoIs 479 239 && 479 % 4 == 3 := by decide +kernel
/-- Small analogues for E-K1. -/
theorem typeII_11 : isPrimeTD 23 && ordTwoIs 23 11 && 23 % 4 == 3 := by decide +kernel
theorem typeII_23 : isPrimeTD 47 && ordTwoIs 47 23 && 47 % 4 == 3 := by decide +kernel

/-! ## Bits per unit step degree (#5) -/

theorem cube_lt_three_pow_shift : ∀ r, (r + 4) ^ 3 < 3 ^ (r + 4)
  | 0 => by decide
  | r + 1 => by
    have ih := cube_lt_three_pow_shift r
    have h1 : (r + 1 + 4) ^ 3 ≤ 3 * (r + 4) ^ 3 := by grind
    have h2 : 3 ^ (r + 1 + 4) = 3 * 3 ^ (r + 4) := by
      rw [show r + 1 + 4 = (r + 4) + 1 by omega, Nat.pow_succ, Nat.mul_comm]
    omega

theorem cube_lt_three_pow (q : Nat) (hq : 4 ≤ q) : q ^ 3 < 3 ^ q := by
  obtain ⟨r, rfl⟩ : ∃ r, q = r + 4 := ⟨q - 4, by omega⟩
  exact cube_lt_three_pow_shift r

/-- **Degree 3 maximizes bits per unit step degree.** For every integer step
degree `q ≥ 2` other than 3, `q^3 < 3^q`, i.e. `log₂ q / q < log₂ 3 / 3`. -/
theorem three_maximizes_bits_per_degree : ∀ q, 2 ≤ q → q ≠ 3 → q ^ 3 < 3 ^ q := by
  intro q h2 h3
  by_cases hq : q = 2
  · subst hq; decide
  · exact cube_lt_three_pow q (by omega)

end FactorBaseDesigns
