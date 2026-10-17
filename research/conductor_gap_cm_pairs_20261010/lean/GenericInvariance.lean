/-!
# Generic-group invariance

A *generic* algorithm for the discrete logarithm sees group elements only
through handles: it may ask the oracle to multiply or invert the elements
behind handles (getting a new handle), or to test two handles for equality
(getting a `Bool`). It never sees a representation. This is Nechaev's
model: the algorithm's control flow depends only on group-operation
outcomes (equality answers) and its random tape. Baby-step giant-step is an
instance of it.

**Scope (correction, 2026-10-16, sheriff case SHF-20261016-3caa6f, L5).**
An earlier version of this comment said Pollard rho and the kangaroo are
instances of this model. The implementations that ecbench measures are
not. They choose each step by hashing the point's representation:
`mix(g.key(&pos))` in the kangaroo (`src/cryptanalysis/ecbench/generic.rs`,
line 275) and `mix(self.g.key(x))` in the rho walk
(`src/cryptanalysis/ic_boundary.rs`, `TunedWalker::index`). Here `Alg.next`
never receives an element or a representation, so the theorem below says
nothing about those walks. The sheriff case ran a toy demonstration, not
on ecbench itself: on `Z_1019` and its isomorphic image in `F_2039^*`, a
representation-hashing r-adding walk had a different cost on 197 of 200
transported instances, and BSGS had a different transcript on 0 of 200.
No theorem in this file is changed.

The theorem: a generic algorithm's whole transcript, hence its operation
count and its answer, is the same on two groups joined by an injective
homomorphism, for every random tape and every step bound. In particular
its cost distribution is the same on two isomorphic groups, so no generic
algorithm in this model can see any invariant of an elliptic curve that is not an
invariant of its group of points: not the j-invariant, not the
endomorphism ring, not its conductor.

Nothing here is specific to groups; only `mul`, `inv` and equality are
used, and the map need only be injective and preserve `mul` and `inv`.
No Mathlib is needed.
-/

namespace GenericInvariance

/-- The operations the oracle offers. No axioms are needed for the theorem. -/
structure GroupOps (α : Type) where
  mul : α → α → α
  inv : α → α

/-- A query: handles are indices into the list of elements the oracle holds. -/
inductive Query where
  | mul (i j : Nat)
  | inv (i : Nat)
  | eq (i j : Nat)
  | halt (answer : Nat)

/-- A generic algorithm chooses its next query from what it has been told:
the number of group operations so far, the equality answers so far, and a
random tape. It is never given an element. -/
structure Alg where
  next : (cost : Nat) → (answers : List Bool) → (tape : Nat → Bool) → Query

/-- The oracle's state: the elements behind the handles, the equality
answers returned, and the number of group operations charged. -/
structure State (α : Type) where
  handles : List α
  answers : List Bool
  cost : Nat

/-- One oracle step: either a new state, or the algorithm halted with
`(cost, answer)`. An out-of-range handle is charged and ignored. -/
def step {α : Type} [DecidableEq α] (ops : GroupOps α) (A : Alg) (tape : Nat → Bool)
    (s : State α) : State α ⊕ (Nat × Nat) :=
  match A.next s.cost s.answers tape with
  | .mul i j =>
    match s.handles[i]?, s.handles[j]? with
    | some a, some b => .inl { s with handles := s.handles ++ [ops.mul a b], cost := s.cost + 1 }
    | _, _ => .inl { s with cost := s.cost + 1 }
  | .inv i =>
    match s.handles[i]? with
    | some a => .inl { s with handles := s.handles ++ [ops.inv a], cost := s.cost + 1 }
    | none => .inl { s with cost := s.cost + 1 }
  | .eq i j =>
    match s.handles[i]?, s.handles[j]? with
    | some a, some b => .inl { s with answers := s.answers ++ [decide (a = b)] }
    | _, _ => .inl { s with answers := s.answers ++ [false] }
  | .halt ans => .inr (s.cost, ans)

/-- Run for at most `n` steps. `none` means the algorithm had not halted. -/
def run {α : Type} [DecidableEq α] (ops : GroupOps α) (A : Alg) (tape : Nat → Bool) :
    Nat → State α → Option (Nat × Nat)
  | 0, _ => none
  | n + 1, s =>
    match step ops A tape s with
    | .inl s' => run ops A tape n s'
    | .inr r => some r

/-- The discrete-logarithm instance `(g, h)`: handles 0 and 1. -/
def start {α : Type} (g h : α) : State α := ⟨[g, h], [], 0⟩

/-- An injective map preserving the two operations. An isomorphism is one;
so is the inclusion of a subgroup. -/
structure Hom {α β : Type} (opsG : GroupOps α) (opsH : GroupOps β) where
  toFun : α → β
  inj : ∀ a b, toFun a = toFun b → a = b
  map_mul : ∀ a b, toFun (opsG.mul a b) = opsH.mul (toFun a) (toFun b)
  map_inv : ∀ a, toFun (opsG.inv a) = opsH.inv (toFun a)

variable {α β : Type} [DecidableEq α] [DecidableEq β]
  {opsG : GroupOps α} {opsH : GroupOps β}

/-- Transport a state along the map. -/
def State.map (φ : Hom opsG opsH) (s : State α) : State β :=
  ⟨s.handles.map φ.toFun, s.answers, s.cost⟩

theorem decide_eq_map (φ : Hom opsG opsH) (a b : α) :
    decide (φ.toFun a = φ.toFun b) = decide (a = b) := by
  by_cases h : a = b
  · subst h; simp
  · have h' : φ.toFun a ≠ φ.toFun b := fun e => h (φ.inj a b e)
    simp [h, h']

/-- One step commutes with transport. -/
theorem step_map (φ : Hom opsG opsH) (A : Alg) (tape : Nat → Bool) (s : State α) :
    step opsH A tape (s.map φ) = Sum.map (State.map φ) id (step opsG A tape s) := by
  unfold step
  simp only [State.map]
  cases A.next s.cost s.answers tape with
  | mul i j =>
    simp only [List.getElem?_map]
    cases s.handles[i]? <;> cases s.handles[j]? <;>
      simp [Sum.map, State.map, φ.map_mul]
  | inv i =>
    simp only [List.getElem?_map]
    cases s.handles[i]? <;> simp [Sum.map, State.map, φ.map_inv]
  | eq i j =>
    simp only [List.getElem?_map]
    cases s.handles[i]? <;> cases s.handles[j]? <;>
      simp [Sum.map, State.map, decide_eq_map φ]
  | halt ans => simp [Sum.map]

/-- **Generic invariance.** For every algorithm, tape, step bound and
start state, the run on the transported state returns exactly the same
`(cost, answer)`, or fails to halt in exactly the same way. -/
theorem run_map (φ : Hom opsG opsH) (A : Alg) (tape : Nat → Bool) :
    ∀ (n : Nat) (s : State α), run opsH A tape n (s.map φ) = run opsG A tape n s := by
  intro n
  induction n with
  | zero => intro s; rfl
  | succ n ih =>
    intro s
    simp only [run, step_map φ]
    cases step opsG A tape s with
    | inl s' => simp [Sum.map, ih]
    | inr r => simp [Sum.map]

/-- The statement for a discrete-logarithm instance: solving `(φ g, φ h)` in
`H` costs exactly what solving `(g, h)` in `G` costs, and returns the same
answer, for every tape. Hence every moment of the cost, and the success
probability, agree. -/
theorem dlog_invariant (φ : Hom opsG opsH) (A : Alg) (tape : Nat → Bool) (n : Nat) (g h : α) :
    run opsH A tape n (start (φ.toFun g) (φ.toFun h)) = run opsG A tape n (start g h) := by
  have : start (φ.toFun g) (φ.toFun h) = (start g h).map φ := rfl
  rw [this, run_map]

end GenericInvariance
