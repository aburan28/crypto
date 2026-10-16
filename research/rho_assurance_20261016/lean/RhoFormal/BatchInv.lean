/-!
# Montgomery batch inversion

Functional models of `WideGf2::batch_inv` (`src/cryptanalysis/wide_gf2m.rs`,
and the identical `WideGf2::batch_inv` in `semaev_decomp.rs`) and of the
four-accumulator `Gf2::batch_inv_inline` / `batch_inv_clmul`
(`src/cryptanalysis/semaev_decomp.rs`), over an arbitrary field
(`Lean.Grind.Field` from core Lean, with decidable equality standing in for the
Rust `x != 0` test on machine words).
-/

namespace RhoFormal.BatchInv

open Lean.Grind

variable {F : Type} [Field F] [DecidableEq F]

omit [DecidableEq F] in
theorem mul_ne_zero' {a b : F} (ha : a ≠ 0) (hb : b ≠ 0) : a * b ≠ 0 := by
  intro h
  have h1 : a * b * b⁻¹ * a⁻¹ = 1 := by
    have := Field.mul_inv_cancel ha
    have := Field.mul_inv_cancel hb
    grind
  rw [h] at h1; grind

/-- The specification: invert nonzero entries, leave zeros unchanged. -/
def spec (x : F) : F := if x = 0 then x else x⁻¹

/-! ## Single accumulator (`WideGf2::batch_inv`) -/

/-- Forward pass from running product `acc`: push `acc` to `scratch`, then
multiply `x` in unless `x = 0`.  Returns `(scratch, final acc)`. -/
def fwd : F → List F → List F × F
  | acc, [] => ([], acc)
  | acc, x :: xs =>
    let r := fwd (if x = 0 then acc else acc * x) xs
    (acc :: r.1, r.2)

/-- Reverse pass over the pairs `(xs[i], scratch[i])`.  The tail (higher
indices) is processed first, as in `for i in (0..len).rev()`.  Returns
`(output, final inv)`. -/
def bwd : F → List (F × F) → List F × F
  | inv, [] => ([], inv)
  | inv, (x, p) :: rest =>
    let r := bwd inv rest
    if x = 0 then (x :: r.1, r.2) else (r.2 * p :: r.1, r.2 * x)

/-- The two-pass algorithm: forward pass from `1`, one inversion of the
total, reverse pass. -/
def batchInv (xs : List F) : List F :=
  let f := fwd 1 xs
  (bwd f.2⁻¹ (xs.zip f.1)).1

theorem fwd_length : ∀ (acc : F) (xs : List F), (fwd acc xs).1.length = xs.length
  | _, [] => rfl
  | acc, x :: xs => by simp [fwd, fwd_length _ xs]

theorem fwd_ne_zero : ∀ (acc : F) (xs : List F), acc ≠ 0 → (fwd acc xs).2 ≠ 0
  | _, [], h => h
  | acc, x :: xs, h => by
    simp only [fwd]
    apply fwd_ne_zero
    split
    · exact h
    · next hx => exact mul_ne_zero' h hx

/-- Invariant: started with the inverse of the forward pass's final product,
the reverse pass outputs `spec` of every entry and ends with `acc⁻¹`. -/
theorem bwd_fwd : ∀ (acc : F) (xs : List F), acc ≠ 0 →
    bwd (fwd acc xs).2⁻¹ (xs.zip (fwd acc xs).1) = (xs.map spec, acc⁻¹)
  | _, [], _ => by simp [fwd, bwd]
  | acc, x :: xs, h => by
    by_cases hx : x = 0
    · have ih := bwd_fwd acc xs h
      simp [fwd, bwd, hx, ih, spec]
    · have hax : acc * x ≠ 0 := mul_ne_zero' h hx
      have ih := bwd_fwd (acc * x) xs hax
      simp only [fwd, hx, ite_false, List.zip_cons_cons, bwd, ih, List.map_cons, spec]
      have h1 := Field.mul_inv_cancel h
      have h2 := Field.mul_inv_cancel hx
      have h3 := Field.mul_inv_cancel hax
      refine Prod.ext ?_ ?_
      · simp only [List.cons.injEq, and_true]; grind
      · grind

/-- **(2) Batch inversion is correct.** For every list over any field
(including the empty list and lists containing zeros), the output has the
same length and `output[i] = xs[i]⁻¹` if `xs[i] ≠ 0`, `output[i] = xs[i]`
(which is `0`) if `xs[i] = 0`. -/
theorem batchInv_correct (xs : List F) :
    (batchInv xs).length = xs.length ∧
    ∀ (i : Nat) (h : i < xs.length) (h' : i < (batchInv xs).length),
      (batchInv xs)[i] = if xs[i] = 0 then xs[i] else xs[i]⁻¹ := by
  have e : batchInv xs = xs.map spec := by
    simp only [batchInv]; rw [bwd_fwd 1 xs (by grind)]
  refine ⟨by simp [e], fun i h h' => ?_⟩
  simp only [e, List.getElem_map, spec]

theorem batchInv_eq_map (xs : List F) : batchInv xs = xs.map spec := by
  simp only [batchInv]; rw [bwd_fwd 1 xs (by grind)]

/-- In a `Lean.Grind.Field`, `0⁻¹ = 0`, so the output is just `xs.map (·⁻¹)`. -/
theorem batchInv_eq_map_inv (xs : List F) : batchInv xs = xs.map (·⁻¹) := by
  rw [batchInv_eq_map]
  apply List.map_congr_left
  intro x _
  by_cases hx : x = 0 <;> simp [spec, hx, Field.inv_zero]

/-! ## Several accumulators (`Gf2::batch_inv_inline` / `batch_inv_clmul`)

Element `i` is dealt to accumulator `lane i`; the Rust uses `lane i = i % 4`
(`chunks(4)` zipped with the four accumulators).  The proof below holds for an
arbitrary lane assignment.  Accumulators are a function `Nat → F`. -/

/-- Set accumulator `l` to `v`, leaving the others. -/
def upd (acc : Nat → F) (l : Nat) (v : F) : Nat → F := fun m => if m = l then v else acc m

def fwdL (lane : Nat → Nat) : Nat → (Nat → F) → List F → List F × (Nat → F)
  | _, acc, [] => ([], acc)
  | i, acc, x :: xs =>
    let r := fwdL lane (i + 1) (if x = 0 then acc else upd acc (lane i) (acc (lane i) * x)) xs
    (acc (lane i) :: r.1, r.2)

/-- Reverse pass with one running inverse per lane, highest index first. -/
def bwdL (lane : Nat → Nat) : Nat → (Nat → F) → List (F × F) → List F × (Nat → F)
  | _, inv, [] => ([], inv)
  | i, inv, (x, p) :: rest =>
    let r := bwdL lane (i + 1) inv rest
    if x = 0 then (x :: r.1, r.2)
    else (r.2 (lane i) * p :: r.1, upd r.2 (lane i) (r.2 (lane i) * x))

theorem fwdL_ne_zero (lane : Nat → Nat) : ∀ (i : Nat) (acc : Nat → F) (xs : List F),
    (∀ l, acc l ≠ 0) → ∀ l, (fwdL lane i acc xs).2 l ≠ 0
  | _, _, [], h => h
  | i, acc, x :: xs, h => by
    simp only [fwdL]
    apply fwdL_ne_zero
    intro l
    split
    · exact h l
    · next hx =>
      simp only [upd]; split
      · exact mul_ne_zero' (h _) hx
      · exact h l

theorem bwdL_fwdL (lane : Nat → Nat) : ∀ (i : Nat) (acc : Nat → F) (xs : List F)
    (inv : Nat → F), (∀ l, acc l ≠ 0) → (∀ l, inv l = ((fwdL lane i acc xs).2 l)⁻¹) →
    (bwdL lane i inv (xs.zip (fwdL lane i acc xs).1)).1 = xs.map spec ∧
    ∀ l, (bwdL lane i inv (xs.zip (fwdL lane i acc xs).1)).2 l = (acc l)⁻¹
  | _, acc, [], inv, _, hinv => by simp [fwdL, bwdL] at *; exact hinv
  | i, acc, x :: xs, inv, h, hinv => by
    by_cases hx : x = 0
    · simp only [fwdL, hx, ite_true] at hinv
      have ih := bwdL_fwdL lane (i + 1) acc xs inv h hinv
      simp only [fwdL, hx, ite_true, List.zip_cons_cons, bwdL, List.map_cons, spec]
      exact ⟨by rw [ih.1], ih.2⟩
    · simp only [fwdL, hx, ite_false] at hinv
      have h' : ∀ l, upd acc (lane i) (acc (lane i) * x) l ≠ 0 := by
        intro l; simp only [upd]; split
        · exact mul_ne_zero' (h _) hx
        · exact h l
      have ih := bwdL_fwdL lane (i + 1) _ xs inv h' hinv
      simp only [fwdL, hx, ite_false, List.zip_cons_cons, bwdL, List.map_cons, spec]
      have ha := h (lane i)
      have hax : acc (lane i) * x ≠ 0 := mul_ne_zero' ha hx
      have h1 := Field.mul_inv_cancel ha
      have h2 := Field.mul_inv_cancel hx
      have h3 := Field.mul_inv_cancel hax
      have hl := ih.2 (lane i)
      simp only [upd, ite_true] at hl
      refine ⟨?_, ?_⟩
      · rw [ih.1, hl]; simp only [List.cons.injEq, and_true]; grind
      · intro l
        simp only [upd]
        split
        · next he => subst he; rw [hl]; grind
        · next hne =>
          have := ih.2 l; simp only [upd, hne, ite_false] at this; exact this

/-- The four lane products are inverted with a single inversion, exactly as
the Rust does: `p01 = a0·a1`, `p23 = a2·a3`, `inv_all = (p01·p23)⁻¹`,
`i01 = inv_all·p23`, `i23 = inv_all·p01`, and the lane inverses are
`[i01·a1, i01·a0, i23·a3, i23·a2]`.  Lanes `≥ 4` are never used; for
totality they get `(acc l)⁻¹`. -/
def combine4 (acc : Nat → F) : Nat → F :=
  let p01 := acc 0 * acc 1
  let p23 := acc 2 * acc 3
  let invAll := (p01 * p23)⁻¹
  let i01 := invAll * p23
  let i23 := invAll * p01
  fun l => if l = 0 then i01 * acc 1 else if l = 1 then i01 * acc 0
    else if l = 2 then i23 * acc 3 else if l = 3 then i23 * acc 2 else (acc l)⁻¹

omit [DecidableEq F] in
theorem combine4_eq_inv (acc : Nat → F) (h : ∀ l, acc l ≠ 0) (l : Nat) :
    combine4 acc l = (acc l)⁻¹ := by
  have h0 := Field.mul_inv_cancel (h 0)
  have h1 := Field.mul_inv_cancel (h 1)
  have h2 := Field.mul_inv_cancel (h 2)
  have h3 := Field.mul_inv_cancel (h 3)
  have h01 := mul_ne_zero' (h 0) (h 1)
  have h23 := mul_ne_zero' (h 2) (h 3)
  have hall := Field.mul_inv_cancel (mul_ne_zero' h01 h23)
  simp only [combine4]
  split
  · next e => subst e; grind
  · split
    · next e => subst e; grind
    · split
      · next e => subst e; grind
      · split
        · next e => subst e; grind
        · rfl

/-- Four-accumulator batch inversion, element `i` in lane `i % 4`. -/
def batchInv4 (xs : List F) : List F :=
  let f := fwdL (· % 4) 0 (fun _ => 1) xs
  (bwdL (· % 4) 0 (combine4 f.2) (xs.zip f.1)).1

/-- **(2, four-lane variant)** Same specification for the `Gf2` algorithm. -/
theorem batchInv4_eq_map (xs : List F) : batchInv4 xs = xs.map spec := by
  have hone : ∀ l : Nat, (fun _ => (1 : F)) l ≠ 0 := fun _ => by grind
  exact (bwdL_fwdL (· % 4) 0 (fun _ => 1) xs _ hone
    (combine4_eq_inv _ (fwdL_ne_zero _ 0 _ xs hone))).1

theorem batchInv4_correct (xs : List F) :
    (batchInv4 xs).length = xs.length ∧
    ∀ (i : Nat) (h : i < xs.length) (h' : i < (batchInv4 xs).length),
      (batchInv4 xs)[i] = if xs[i] = 0 then xs[i] else xs[i]⁻¹ := by
  refine ⟨by simp [batchInv4_eq_map], fun i h h' => ?_⟩
  simp only [batchInv4_eq_map, List.getElem_map, spec]

end RhoFormal.BatchInv
