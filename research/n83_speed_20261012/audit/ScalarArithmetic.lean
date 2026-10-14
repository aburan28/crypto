import Init

namespace N83ScalarArithmetic

/-- Exact distributive expansion underlying the four-limb implementation. -/
theorem limb_product_expansion (a0 a1 b0 b1 w : Nat) :
    (a0 + w * a1) * (b0 + w * b1) =
      a0 * b0 + w * (a0 * b1 + a1 * b0) + w * w * a1 * b1 := by
  simp [Nat.add_mul, Nat.mul_add, Nat.mul_assoc, Nat.mul_comm, Nat.mul_left_comm,
    Nat.add_assoc, Nat.add_comm, Nat.add_left_comm]

/-- For any radix B and any nonzero modulus m, the floor reciprocal gives a
    quotient whose residual is nonnegative and below 2m whenever x < B. -/
theorem barrett_residual_bounds (B m x : Nat) (hB : 0 < B) (hm : 0 < m)
    (hx : x < B) :
    let mu := B / m
    let q := x * mu / B
    q * m ≤ x ∧ x - q * m < 2 * m := by
  dsimp
  let mu := B / m
  let q := x * mu / B
  let s := B % m
  let t := (x * mu) % B
  have hmu : mu * m ≤ B := Nat.div_mul_le_self B m
  have hq : q * B ≤ x * mu := Nat.div_mul_le_self (x * mu) B
  have hprod : (q * m) * B ≤ x * B := by
    calc
      (q * m) * B = (q * B) * m := by simp [Nat.mul_assoc, Nat.mul_comm, Nat.mul_left_comm]
      _ ≤ (x * mu) * m := Nat.mul_le_mul_right m hq
      _ = x * (mu * m) := by simp [Nat.mul_assoc]
      _ ≤ x * B := Nat.mul_le_mul_left x hmu
  have hqx : q * m ≤ x := Nat.le_of_mul_le_mul_right hprod hB
  constructor
  · exact hqx
  · have hs : s < m := Nat.mod_lt B hm
    have ht : t < B := Nat.mod_lt (x * mu) hB
    have hBd : m * mu + s = B := Nat.div_add_mod B m
    have hXd : B * q + t = x * mu := Nat.div_add_mod (x * mu) B
    have hid : B * x = B * (q * m) + (m * t + s * x) := by
      calc
        B * x = (m * mu + s) * x := by rw [hBd]
        _ = m * (x * mu) + s * x := by simp [Nat.add_mul, Nat.mul_add, Nat.mul_assoc, Nat.mul_comm, Nat.mul_left_comm]
        _ = m * (B * q + t) + s * x := by rw [hXd]
        _ = B * (q * m) + (m * t + s * x) := by simp [Nat.mul_add, Nat.mul_assoc, Nat.mul_comm, Nat.mul_left_comm, Nat.add_assoc]
    have hre : (x - q * m) + q * m = x := Nat.sub_add_cancel hqx
    have hsum : B * (x - q * m) + B * (q * m) = B * (q * m) + (m * t + s * x) := by
      calc
        B * (x - q * m) + B * (q * m) = B * ((x - q * m) + q * m) := by rw [Nat.mul_add]
        _ = B * x := by rw [hre]
        _ = B * (q * m) + (m * t + s * x) := hid
    have hres : B * (x - q * m) = m * t + s * x := by omega
    have hfirst : m * t < m * B := by
      simpa [Nat.mul_comm] using Nat.mul_lt_mul_of_pos_right ht hm
    have hsecond : s * x < m * B := by
      calc
        s * x ≤ m * x := Nat.mul_le_mul_right x (Nat.le_of_lt hs)
        _ < m * B := by simpa [Nat.mul_comm] using Nat.mul_lt_mul_of_pos_right hx hm
    have hless : B * (x - q * m) < m * B + m * B := by omega
    have hfact : (x - q * m) * B < (2 * m) * B := by
      rw [show 2 * m = m + m by omega, Nat.add_mul]
      simpa only [Nat.mul_comm] using hless
    exact Nat.lt_of_mul_lt_mul_right hfact

/-- A single conditional subtraction returns the exact remainder, for every
    Nat modulus/radix satisfying the explicit floor-reciprocal hypotheses. -/
theorem barrett_one_correction (B m x : Nat) (hB : 0 < B) (hm : 0 < m)
    (hx : x < B) :
    let q := x * (B / m) / B
    let r := x - q * m
    (if m ≤ r then r - m else r) = x % m := by
  dsimp
  let q := x * (B / m) / B
  let r := x - q * m
  have hb := barrett_residual_bounds B m x hB hm hx
  change q * m ≤ x ∧ r < 2 * m at hb
  have hre : r + q * m = x := Nat.sub_add_cancel hb.1
  have hmod : x % m = r % m := by
    rw [← hre]
    exact Nat.add_mul_mod_self_right r q m
  change (if m ≤ r then r - m else r) = x % m
  rw [hmod]
  split
  next hr =>
    have hsmall : r - m < m := by omega
    rw [Nat.mod_eq_sub_mod hr, Nat.mod_eq_of_lt hsmall]
  next hr =>
    have hsmall : r < m := by omega
    exact (Nat.mod_eq_of_lt hsmall).symm

/-- Canonical modular addition requires at most one subtraction. -/
theorem canonical_add_one_subtraction (m a b : Nat) (ha : a < m) (hb : b < m) :
    (if m ≤ a + b then a + b - m else a + b) = (a + b) % m := by
  split
  next hs =>
    have hsmall : a + b - m < m := by omega
    rw [Nat.mod_eq_sub_mod hs, Nat.mod_eq_of_lt hsmall]
  next hs =>
    have hsmall : a + b < m := by omega
    exact (Nat.mod_eq_of_lt hsmall).symm

end N83ScalarArithmetic
