import FactorBaseDesigns.Dickson

/-!
# Correctness of the new membership presentations

These theorems show that each presentation in the design register encodes
the set it claims to. The experiments then measure what the presentation
costs.

* **#1, exponent intervals** (`exp_interval_identity`): for bits `b_j`,
  `∏_j (1 + (g^{2^j} − 1)·b_j) = g^{Σ b_j 2^j}` in any commutative ring.
  So the bit-product system with partial products cuts out exactly
  `{g^k : 0 ≤ k < 2^t}`, with no condition on the order of `g` (no smoothness
  of `p ± 1`).
* **#2, Lucas ladder** (`lucas_double`, `lucas_double_succ`):
  `V_{2k} = V_k² − 2` and `V_{2k+1} = V_{k+1}V_k − V_1` for `V_k = γ^k + γ^{−k}`.
  These are the ladder steps of the torus exponent interval.
* **#11, norm chain** (`chain_sound`): if `G(a) = 0` and the evaluated chain
  `r_0 = a`, `r_j² = q_j·G(a) + r_{j+1}`, `r_t = 1` holds, then `a^{2^t} = 1`.
  Applied at each root `a = x_i/g` of `G(gX)/g²`, every `x_i` lies in
  `g·μ_{2^t}`. Bridge: the chain equations are polynomial identities in `X`,
  and evaluation at `a` is a ring homomorphism. This is soundness only.
  Completeness (that witnesses exist when the roots lie in the fibre) is
  Euclidean division and is checked on every planted instance by the harness.
* **#11, `m = 2` pencil** (`pencil_norm_factor`): on `y² = x³ + Ax + B`,
  the line through `−R` with slope `λ` has
  `X³ + AX + B − (λ(X − x_R) + y_R)² = (X − x_R)(X² + e₁X + e₀)`, where
  `e₁ = x_R − λ²` and `e₀ = x_R² + λ²x_R + A − 2λy_R`.
* **#11 counting** (`norm_chain_excess`): with `m − 1` free coefficients of
  `f`, `m` coefficients of `G`, and `2m − 1` chain unknowns and equations per
  step, there is exactly one more equation than unknown.
-/

namespace FactorBaseDesigns

variable {α : Type} [Lean.Grind.CommRing α]

theorem sq_pow (g : α) (v : Nat) : (g * g) ^ v = g ^ (2 * v) := by
  induction v with
  | zero => simp [Lean.Grind.Semiring.pow_zero]
  | succ v ih =>
    rw [Lean.Grind.Semiring.pow_succ, ih, show 2 * (v + 1) = 2 * v + 2 by omega,
      Lean.Grind.Semiring.pow_add, Lean.Grind.Semiring.pow_two]

theorem self_mul_pow (g : α) (e : Nat) : g ^ e * g ^ e = g ^ (2 * e) := by
  rw [← Lean.Grind.Semiring.pow_add, show e + e = 2 * e by omega]

/-- Integer value of a little-endian bit list. -/
def bitsVal : List Bool → Nat
  | [] => 0
  | b :: bs => (if b then 1 else 0) + 2 * bitsVal bs

/-- `∏_j (1 + (h_j − 1)·b_j)` where `h_0 = g` and `h_{j+1} = h_j²`. -/
def expProd (g : α) : List Bool → α
  | [] => 1
  | b :: bs => (1 + (g - 1) * (if b then 1 else 0)) * expProd (g * g) bs

/-- The base used at position `j` is `g^{2^j}`. -/
theorem iterate_sq (g : α) (j : Nat) : Nat.repeat (fun x => x * x) j g = g ^ (2 ^ j) := by
  induction j with
  | zero => simp [Nat.repeat, Lean.Grind.Semiring.pow_one]
  | succ j ih =>
    simp only [Nat.repeat, ih]
    rw [self_mul_pow, Nat.pow_succ, Nat.mul_comm]

/-- **#1.** The bit-product presentation equals `g^{bitsVal bs}`. -/
theorem exp_interval_identity (g : α) (bs : List Bool) : expProd g bs = g ^ bitsVal bs := by
  induction bs generalizing g with
  | nil => simp [expProd, bitsVal, Lean.Grind.Semiring.pow_zero]
  | cons b bs ih =>
    simp only [expProd, bitsVal, ih]
    rw [sq_pow, Lean.Grind.Semiring.pow_add]
    cases b <;> simp [Lean.Grind.Semiring.pow_zero, Lean.Grind.Semiring.pow_one] <;> grind

/-- `t` bits give exactly the exponents `0 ≤ k < 2^t`. -/
theorem bitsVal_lt (bs : List Bool) : bitsVal bs < 2 ^ bs.length := by
  induction bs with
  | nil => simp [bitsVal]
  | cons b bs ih => simp only [bitsVal, List.length_cons, Nat.pow_succ]; split <;> omega

/-- **#2.** Lucas ladder, doubling step. -/
theorem lucas_double (γ γ' : α) (h : γ * γ' = 1) (k : Nat) :
    gp γ γ' (2 * k) = gp γ γ' k * gp γ γ' k - 2 := by
  have := dickson_mul γ γ' h (i := k) (j := k) (Nat.le_refl _)
  rw [Nat.sub_self] at this
  have e : gp γ γ' 0 = 2 := by
    unfold gp; simp only [Lean.Grind.Semiring.pow_zero]; grind
  rw [e, show k + k = 2 * k by omega] at this
  grind

/-- **#2.** Lucas ladder, doubling-plus-one step. -/
theorem lucas_double_succ (γ γ' : α) (h : γ * γ' = 1) (k : Nat) :
    gp γ γ' (2 * k + 1) = gp γ γ' (k + 1) * gp γ γ' k - gp γ γ' 1 := by
  have := dickson_mul γ γ' h (i := k + 1) (j := k) (Nat.le_succ _)
  rw [show k + 1 - k = 1 by omega, show k + 1 + k = 2 * k + 1 by omega] at this
  grind

/-- **#11, soundness at a root.** -/
theorem chain_sound (a Ga : α) (r q : Nat → α) (t : Nat) (hG : Ga = 0) (h0 : r 0 = a)
    (hstep : ∀ j, j < t → r j * r j = q j * Ga + r (j + 1)) (ht : r t = 1) :
    a ^ (2 ^ t) = 1 := by
  have key : ∀ j, j ≤ t → r j = a ^ (2 ^ j) := by
    intro j
    induction j with
    | zero => intro _; rw [h0]; simp [Lean.Grind.Semiring.pow_one]
    | succ j ih =>
      intro hj
      have hs := hstep j (by omega)
      rw [hG, ih (by omega), self_mul_pow] at hs
      rw [Nat.pow_succ, Nat.mul_comm]
      grind
  rw [← key t (Nat.le_refl _), ht]

/-- **#11, `m = 2` pencil factorization.** -/
theorem pencil_norm_factor (A B xR yR lam X : α) (hR : yR * yR = xR ^ 3 + A * xR + B) :
    X ^ 3 + A * X + B - (lam * (X - xR) + yR) ^ 2 =
      (X - xR) * (X ^ 2 + (xR - lam ^ 2) * X + (xR ^ 2 + lam ^ 2 * xR + A - 2 * lam * yR)) := by
  grind

/-- Unknowns of the degree-`m` norm chain of depth `t`. -/
def chainUnknowns (m t : Nat) : Nat := (m - 1) + m + t * (2 * m - 1)
/-- Equations: `m` defining `G`, `2m − 1` per step, `m` for `r_t = 1`. -/
def chainEquations (m t : Nat) : Nat := m + t * (2 * m - 1) + m

/-- **#11 counting.** One more equation than unknown, for every `m ≥ 1`
and `t`: exactly the single `F_p`-condition that makes relations rare. -/
theorem norm_chain_excess (m t : Nat) (hm : 1 ≤ m) :
    chainEquations m t = chainUnknowns m t + 1 := by
  unfold chainEquations chainUnknowns; omega

end FactorBaseDesigns
