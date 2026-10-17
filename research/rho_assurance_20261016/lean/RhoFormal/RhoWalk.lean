/-!
# Rho walk algebra

**Modelling assumption.**  The prime-order cyclic subgroup `⟨G⟩` is modelled
by the scalar field `F` itself: the point `[x]G` is represented by `x`, so
`G ↦ 1`, `Q = [d]G ↦ d`, point addition is field addition, and `[m]P` is
`m * P`.  This identification is a group isomorphism `Z/rZ ≅ ⟨G⟩` when
`F = Z/rZ`, `r` prime; the bridge to actual curve points (that the curve
arithmetic in `StrongRho` implements this group, and that the
Frobenius/negation automorphism acts on `⟨G⟩` as `[±λ^k]`) is *not*
formalized here.  `F` is any `Lean.Grind.Field`.

Rust (`src/cryptanalysis/koblitz_strong_rho.rs`) keeps a state
`(point, a, b)` meaning `point = [a]G + [b]Q`.
-/

namespace RhoFormal.RhoWalk

open Lean.Grind

variable {F : Type} [Field F]

/-- A walk state: point `P` with coefficients `a, b`. -/
structure St (F : Type) where
  P : F
  a : F
  b : F

/-- Representation invariant `P = a·G + b·Q`, i.e. `P = a + b·d`. -/
def Rep (d : F) (s : St F) : Prop := s.P = s.a + s.b * d

/-! ## Individual algebraic facts -/

/-- **(3a) Jump step.** `P = a + b·d` and jump point `J = c` (a jump is
`JumpG { point: [c]G, a: c }`) give `P + J = (a + c) + b·d`. -/
theorem jump_step (d a b c P : F) (h : P = a + b * d) : P + c = (a + c) + b * d := by
  rw [h]; grind

/-- **(3b) Canonicalization.** An automorphism acting on `⟨G⟩` as `[m]`
(`m = ±λ^k`) maps `P = a + b·d` to `(m·a) + (m·b)·d`. -/
theorem canon_step (d a b m P : F) (h : P = a + b * d) : m * P = (m * a) + (m * b) * d := by
  rw [h]; grind

/-- Composition of two multipliers is their product (on points and on the
coefficient pair). -/
theorem canon_compose (m₁ m₂ P a b : F) :
    m₁ * (m₂ * P) = (m₁ * m₂) * P ∧ m₁ * (m₂ * a) = (m₁ * m₂) * a ∧
    m₁ * (m₂ * b) = (m₁ * m₂) * b := by
  refine ⟨?_, ?_, ?_⟩ <;> grind

/-- Negation is the multiplier `m = -1`. -/
theorem canon_neg (d a b P : F) (h : P = a + b * d) :
    -P = (-1 : F) * P ∧ -P = (-a) + (-b) * d := by
  rw [h]; constructor <;> grind

/-- **(3c) Fruitless-cycle escape.** `2P = 2a + 2b·d`; the Rust doubles as
`add(P, P)`, `a + a`, `b + b`. -/
theorem escape_step (d a b P : F) (h : P = a + b * d) : P + P = (a + a) + (b + b) * d := by
  rw [h]; grind

/-- **(3d) Collision.** Two representations of the same point with
`b ≠ b'` determine `d`.  The lane state is `(a, b)`, the table hit is
`(a', b') = (hit_a, hit_b)`, and the Rust computes
`candidate = (hit_a - a) * (b - hit_b)⁻¹`. -/
theorem collision (d a b a' b' : F) (h : a + b * d = a' + b' * d) (hb : b ≠ b') :
    d = (a' - a) * (b - b')⁻¹ := by
  have hne : b - b' ≠ 0 := by intro h0; apply hb; grind
  have := Field.mul_inv_cancel hne
  grind

/-- Hence the check `[candidate]G = Q`, i.e. `candidate * 1 = d`, succeeds. -/
theorem candidate_check (d a b a' b' : F) (h : a + b * d = a' + b' * d) (hb : b ≠ b') :
    ((a' - a) * (b - b')⁻¹) * 1 = d := by
  rw [← collision d a b a' b' h hb]; grind

/-- If `b = b'` the merge is wasted, not wrong: then `a = a'`.  (The Rust
`inverse_mod(b - hit_b)` fails exactly then: `b - b' = 0 ↔ b = b'`.) -/
theorem collision_same_b (d a b a' b' : F) (h : a + b * d = a' + b' * d) (hb : b = b') :
    a = a' := by
  subst hb; grind

theorem denominator_zero_iff (b b' : F) : b - b' = 0 ↔ b = b' := by
  constructor <;> intro h <;> grind

/-! ## The invariant as one theorem -/

/-- The start cursor: `cursor = [c]G + Q`, `cursor_a = c`.  Invariant
`cursor = cursor_a + d`. -/
def CursorRep (d C Ca : F) : Prop := C = Ca + d

/-- Every state the walk can produce, as an inductive closure of the Rust
operations.  `start`: lane start `(cursor, cursor_a, 1)` from a cursor
satisfying `CursorRep`; `step`: add the jump `[c]G`, `a += c`;
`canon`: multiply point and both coefficients by `m` (any `m`; the Rust uses
`±λ^k`, and leaves the point at infinity unchanged, i.e. `m = 1`);
`escape`: double point and both coefficients. -/
inductive Reachable (d : F) : St F → Prop
  | start (C Ca : F) : CursorRep d C Ca → Reachable d ⟨C, Ca, 1⟩
  | step (s : St F) (c : F) : Reachable d s → Reachable d ⟨s.P + c, s.a + c, s.b⟩
  | canon (s : St F) (m : F) : Reachable d s → Reachable d ⟨m * s.P, m * s.a, m * s.b⟩
  | escape (s : St F) : Reachable d s → Reachable d ⟨s.P + s.P, s.a + s.a, s.b + s.b⟩

/-- The initial cursor `[c]G + Q` satisfies `CursorRep`, and advancing by the
stride `[σ]G` (`cursor_a += σ`) preserves it. -/
theorem cursor_init (d c : F) : CursorRep d (c * 1 + d) c := by
  simp only [CursorRep]; grind

theorem cursor_advance (d C Ca σ : F) (h : CursorRep d C Ca) : CursorRep d (C + σ) (Ca + σ) := by
  simp only [CursorRep] at *; rw [h]; grind

/-- **(3e) Invariant preservation and collision, as one theorem.**
(i) every reachable state satisfies `Rep d`; (ii) the individual operations
preserve `Rep d`; (iii) two reachable states at the same point with
`b ≠ b'` yield `d = (a' - a)·(b - b')⁻¹`, and with `b = b'` they have
`a = a'`. -/
theorem rho_invariant (d : F) :
    (∀ s, Reachable d s → Rep d s) ∧
    (∀ C Ca, CursorRep d C Ca → Rep d ⟨C, Ca, 1⟩) ∧
    (∀ s c, Rep d s → Rep d ⟨s.P + c, s.a + c, s.b⟩) ∧
    (∀ s m, Rep d s → Rep d ⟨m * s.P, m * s.a, m * s.b⟩) ∧
    (∀ s, Rep d s → Rep d ⟨s.P + s.P, s.a + s.a, s.b + s.b⟩) ∧
    (∀ s s', Reachable d s → Reachable d s' → s.P = s'.P →
      (s.b ≠ s'.b → d = (s'.a - s.a) * (s.b - s'.b)⁻¹) ∧ (s.b = s'.b → s.a = s'.a)) := by
  have hstart : ∀ C Ca, CursorRep d C Ca → Rep d ⟨C, Ca, 1⟩ := by
    intro C Ca h; simp only [Rep, CursorRep] at *; rw [h]; grind
  have hstep : ∀ (s : St F) c, Rep d s → Rep d ⟨s.P + c, s.a + c, s.b⟩ :=
    fun s c h => jump_step d s.a s.b c s.P h
  have hcanon : ∀ (s : St F) m, Rep d s → Rep d ⟨m * s.P, m * s.a, m * s.b⟩ :=
    fun s m h => canon_step d s.a s.b m s.P h
  have hesc : ∀ (s : St F), Rep d s → Rep d ⟨s.P + s.P, s.a + s.a, s.b + s.b⟩ :=
    fun s h => escape_step d s.a s.b s.P h
  have hreach : ∀ s, Reachable d s → Rep d s := by
    intro s hs
    induction hs with
    | start C Ca h => exact hstart C Ca h
    | step s c _ ih => exact hstep s c ih
    | canon s m _ ih => exact hcanon s m ih
    | escape s _ ih => exact hesc s ih
  refine ⟨hreach, hstart, hstep, hcanon, hesc, ?_⟩
  intro s s' hs hs' hP
  have r := hreach s hs
  have r' := hreach s' hs'
  simp only [Rep] at r r'
  have heq : s.a + s.b * d = s'.a + s'.b * d := by rw [← r, ← r', hP]
  exact ⟨collision d _ _ _ _ heq, collision_same_b d _ _ _ _ heq⟩

end RhoFormal.RhoWalk
