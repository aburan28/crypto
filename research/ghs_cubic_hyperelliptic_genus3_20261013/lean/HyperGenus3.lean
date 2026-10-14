/-!
# The `m = 3, r = 5` GHS configuration: combinatorics and the algebraic identities of the model

Points of `P¹` involved: `a, b` (τ-fixed) and `r0, r1, r2` (one τ-orbit), encoded as bits
`0..4` of a `Nat` mask. `τ` fixes `a, b` and cycles `r0 → r1 → r2 → r0`.

Part A (decide): branch-set combinatorics of `cover_census.rs` for `B = {a, b, r0, r1}`:
`r = 5`, `m = 3`, genus `3`, the rank-2 subspace with union of size 3 (hyperelliptic test),
the branch locus `{a, b}` of `C → D`, and `τB ∆ τ²B ⊆ B` (the étale double cover of `E0`
that puts `T1` in the kernel of `π0^*`).

Part B (grind, any commutative ring): the polynomial identities behind the conic model of `D`,
`F0F1F2 = Q³N²`, and the hyperelliptic equation `w² = Q(x)·xd²`.

No Mathlib. Not formalized: the function fields, the descent itself, Frobenius, the
L-polynomial, the transfer, and `π0_* ∘ φ = [4]`.
-/

namespace HyperGenus3

/-! ## Part A: combinatorics -/

def bit (i : Nat) : Nat := 2 ^ i
def a : Nat := bit 0
def b : Nat := bit 1
def r0 : Nat := bit 2
def r1 : Nat := bit 3
def r2 : Nat := bit 4

/-- τ on a point index. -/
def tauPt : Nat → Nat
  | 2 => 3 | 3 => 4 | 4 => 2 | i => i

def popcount (n : Nat) : Nat := (List.range 5).foldl (fun s i => s + (n / 2 ^ i) % 2) 0
def tauSet (s : Nat) : Nat := (List.range 5).foldl (fun t i => if (s / 2 ^ i) % 2 = 1 then t ||| bit (tauPt i) else t) 0

def B : Nat := a ||| b ||| r0 ||| r1
def tB : Nat := tauSet B
def ttB : Nat := tauSet tB

theorem tau_order_three : ∀ i, i < 5 → tauPt (tauPt (tauPt i)) = i := by decide
theorem tB_eq : tB = a ||| b ||| r1 ||| r2 := by decide
theorem ttB_eq : ttB = a ||| b ||| r2 ||| r0 := by decide
theorem tttB_eq : tauSet ttB = B := by decide

/-- `r = |B ∪ τB ∪ τ²B| = 5`. -/
theorem r_eq_five : popcount (B ||| tB ||| ttB) = 5 := by decide

/-- The seven non-trivial `F₂`-combinations of `B, τB, τ²B`. -/
def combos : List Nat := [B, tB, ttB, B ^^^ tB, B ^^^ ttB, tB ^^^ ttB, B ^^^ tB ^^^ ttB]

/-- `m = 3`: no non-trivial combination is empty. -/
theorem m_eq_three : combos.all (· ≠ 0) = true := by decide

/-- Subcover genera `(|S| − 2)/2` sum to `3 = 2^(m−2)(r−4)+1`. -/
theorem genus_sum : (combos.map (fun s => (popcount s - 2) / 2)).foldl (· + ·) 0 = 3 ∧
    2 ^ (3 - 2) * (5 - 4) + 1 = 3 := by decide

/-- The rank-2 subspace `W = {B∆τB, B∆τ²B, τB∆τ²B}` is closed and its union is `{r0, r1, r2}`. -/
theorem hyperelliptic_subspace :
    (B ^^^ tB) ^^^ (B ^^^ ttB) = tB ^^^ ttB ∧
    (B ^^^ tB) ||| (B ^^^ ttB) ||| (tB ^^^ ttB) = r0 ||| r1 ||| r2 ∧
    popcount ((B ^^^ tB) ||| (B ^^^ ttB) ||| (tB ^^^ ttB)) = 3 := by decide

/-- The triple combination is `{a, b}`: the branch locus of `C → D` lies over `a, b`. -/
theorem triple_is_ab : B ^^^ tB ^^^ ttB = a ||| b := by decide

/-- `τB ∆ τ²B = {r0, r1} ⊆ B`: the `u12` double cover of `E0` is unramified. -/
theorem etale_subcover : tB ^^^ ttB = r0 ||| r1 ∧ (tB ^^^ ttB) &&& B = tB ^^^ ttB := by decide

/-- Descent datum on indices: `y_i ↦ y_{i+1}` and `(i, i+1) ↦ (i+1, i+2)` have order 3. -/
theorem descent_cycle : ∀ i, i < 3 → (i + 3) % 3 = i % 3 ∧ ((i + 1) % 3 + 1) % 3 = (i + 2) % 3 := by decide

/-! ## Part B: polynomial identities over a commutative ring -/

section Alg
variable {R : Type} [Lean.Grind.CommRing R]

/-- The conic relation `(r2−r0)(x−r1) + (r0−r1)(x−r2) + (r1−r2)(x−r0) = 0`. -/
theorem conic_coeffs (x s0 s1 s2 : R) :
    (s2 - s0) * (x - s1) + (s0 - s1) * (x - s2) + (s1 - s2) * (x - s0) = 0 := by grind

/-- From `u01² = (x−r0)(x−r2)`, `u12² = (x−r1)(x−r0)`, `u20² = (x−r2)(x−r1)` the anticanonical
coordinates satisfy `X0² = N(x−r1)`, `X1² = N(x−r2)`, `X2² = N(x−r0)`, the conic, and
`X2²(x−r1) = X0²(x−r0)` (which recovers `x`). -/
theorem conic_model (x s0 s1 s2 u01 u12 u20 : R)
    (h01 : u01 ^ 2 = (x - s0) * (x - s2)) (h12 : u12 ^ 2 = (x - s1) * (x - s0))
    (h20 : u20 ^ 2 = (x - s2) * (x - s1)) :
    let N := (x - s0) * (x - s1) * (x - s2)
    (u12 * u20) ^ 2 = N * (x - s1) ∧ (u20 * u01) ^ 2 = N * (x - s2) ∧ (u01 * u12) ^ 2 = N * (x - s0) ∧
    (s2 - s0) * (u12 * u20) ^ 2 + (s0 - s1) * (u20 * u01) ^ 2 + (s1 - s2) * (u01 * u12) ^ 2 = 0 ∧
    (u01 * u12) ^ 2 * (x - s1) = (u12 * u20) ^ 2 * (x - s0) := by
  intro N
  refine ⟨?_, ?_, ?_, ?_, ?_⟩ <;> grind

/-- `F0 F1 F2 = Q³ N²`. -/
theorem product_of_conjugates (x s0 s1 s2 Q : R) :
    (Q * (x - s0) * (x - s1)) * (Q * (x - s1) * (x - s2)) * (Q * (x - s2) * (x - s0)) =
      Q ^ 3 * ((x - s0) * (x - s1) * (x - s2)) ^ 2 := by grind

/-- The hyperelliptic equation, in multiplied-out form: if `y_i² = F_i(x)`,
`w·Q·N = y0y1y2·xd` and `x·xd = xn`, then `(w Q N)² = (Q N)²·h` with
`h = xn² + c1 xn xd + c0 xd²` and `Q = x² + c1 x + c0`. Over a field with `Q N ≠ 0`
this is `w² = h(z)`. -/
theorem hyperelliptic_equation (x s0 s1 s2 c0 c1 y0 y1 y2 w xn xd : R)
    (hy0 : y0 ^ 2 = (x ^ 2 + c1 * x + c0) * (x - s0) * (x - s1))
    (hy1 : y1 ^ 2 = (x ^ 2 + c1 * x + c0) * (x - s1) * (x - s2))
    (hy2 : y2 ^ 2 = (x ^ 2 + c1 * x + c0) * (x - s2) * (x - s0))
    (hw : w * ((x ^ 2 + c1 * x + c0) * ((x - s0) * (x - s1) * (x - s2))) = y0 * y1 * y2 * xd)
    (hx : x * xd = xn) :
    (w * ((x ^ 2 + c1 * x + c0) * ((x - s0) * (x - s1) * (x - s2)))) ^ 2 =
      ((x ^ 2 + c1 * x + c0) * ((x - s0) * (x - s1) * (x - s2))) ^ 2 *
        (xn ^ 2 + c1 * xn * xd + c0 * xd ^ 2) := by
  rw [hw, ← hx]
  have : (y0 * y1 * y2 * xd) ^ 2 = y0 ^ 2 * y1 ^ 2 * y2 ^ 2 * xd ^ 2 := by grind
  rw [this, hy0, hy1, hy2]
  generalize hQ : x ^ 2 + c1 * x + c0 = Q
  generalize hN : (x - s0) * (x - s1) * (x - s2) = N
  have hF : Q * (x - s0) * (x - s1) * (Q * (x - s1) * (x - s2)) * (Q * (x - s2) * (x - s0)) = Q ^ 3 * N ^ 2 := by
    rw [← hN]; exact product_of_conjugates x s0 s1 s2 Q
  have hH : (x * xd) ^ 2 + c1 * (x * xd) * xd + c0 * xd ^ 2 = Q * xd ^ 2 := by rw [← hQ]; grind
  rw [hF, hH]
  grind

/-- Addendum: the norm identity behind a reduction step of the real model,
`f − v_a² = (f − V₊²) − (v_a − V₊)(v_a + V₊)`. -/
theorem step_norm (f va vp : R) : f - va ^ 2 = (f - vp ^ 2) - (va - vp) * (va + vp) := by grind

end Alg

/-! ## Part C (addendum 2026-10-14): integer bookkeeping of the balanced real-model arithmetic, genus 3

A class is `A + n∞₊ + m∞₋ − D∞` with `deg A + n + m = 3` and `D∞ = 2∞₊ + ∞₋`. -/

/-- Reduction step: if `div(y − v_a) = A + A' − e₊∞₊ − e₋∞₋` (so `degA + degA' = e₊ + e₋`), then
`A ~ ιA' + (e₊ − degA')∞₊ + (e₋ − degA')∞₋`, hence `n' = n + degA − e₋`, and the new
coefficients again sum to `3 − degA'`. -/
theorem step_bookkeeping (degA degA' ep em n m : Int) (hdiv : degA + degA' = ep + em)
    (hsum : degA + n + m = 3) :
    n + (ep - degA') = n + degA - em ∧ degA' + (n + (ep - degA')) + (m + (em - degA')) = 3 := by omega

/-- Composition: `A₁ + A₂ = A + (pairs over d)` and each pair is `~ ∞₊ + ∞₋`; with `−2D∞ = −D∞ − 2∞₊ − ∞₋`
the new coefficients are `n₁ + n₂ + deg d − 2` and `m₁ + m₂ + deg d − 1`, and they sum correctly. -/
theorem compose_bookkeeping (d1 d2 dd n1 m1 n2 m2 : Int) (h1 : d1 + n1 + m1 = 3) (h2 : d2 + n2 + m2 = 3) :
    (d1 + d2 - 2 * dd) + (n1 + n2 + dd - 2) + (m1 + m2 + dd - 1) = 3 := by omega

/-- Negation: `−(A + n∞₊ + m∞₋ − D∞) ~ ιA + (4 − degA − n)∞₊ + (2 − degA − m)∞₋ − D∞`. -/
theorem neg_bookkeeping (dA n m : Int) (h : dA + n + m = 3) : dA + (4 - dA - n) + (2 - dA - m) = 3 := by omega

/-- The transfer's correction term: twelve classes `[Q_k − ∞₊]` plus `6·[∞₊ − ∞₋]` equal
`Σ Q_k − 6∞₊ − 6∞₋` (coefficients at infinity). -/
theorem transfer_infinity : (-12 + 6 : Int) = -6 ∧ (0 - 6 : Int) = -6 := by decide

end HyperGenus3

#print axioms HyperGenus3.r_eq_five
#print axioms HyperGenus3.m_eq_three
#print axioms HyperGenus3.genus_sum
#print axioms HyperGenus3.hyperelliptic_subspace
#print axioms HyperGenus3.triple_is_ab
#print axioms HyperGenus3.etale_subcover
#print axioms HyperGenus3.tau_order_three
#print axioms HyperGenus3.descent_cycle
#print axioms HyperGenus3.conic_coeffs
#print axioms HyperGenus3.conic_model
#print axioms HyperGenus3.product_of_conjugates
#print axioms HyperGenus3.hyperelliptic_equation
#print axioms HyperGenus3.step_norm
#print axioms HyperGenus3.step_bookkeeping
#print axioms HyperGenus3.compose_bookkeeping
#print axioms HyperGenus3.neg_bookkeeping
#print axioms HyperGenus3.transfer_infinity
