/-!
# Fall-degree bounds: the scoring rules as Lean functions

The registered rules of `PREREGISTRATION.md` that are integer statements:

- **P3a** (E3): `d_F ≤ 5` on every draw with `n ≤ 32`.
- **P3b** (E3): between `n = 18` and `n = 32` the per-`n` median `d_F` steps
  from 4 to 5 at most once (i.e. the median sequence, in `n` order, is
  non-decreasing and changes value at most once).
- **P4a** (E4): `d_F ≤ D*` on every draw; **P4b**: `d_F < D*` on at least
  half the draws.
- **Amendment 4/6 readings** (E6, exploratory): "bounded at toy size" if
  every unsatisfiable draw has `d_F ≤ 4` and every satisfiable draw has
  `d_last ≤ 4`; "Semaev-specific" (Part D′) if every control shape has
  median `d_last ≥ 5`.

The LOG/LIN model comparison of E3 uses logarithms and is computed in Rust,
not here. A censored value (`none`) is "above the cap", recorded as the cap
plus one, so a censored draw can never satisfy a `≤` bound it might break.
-/

namespace FallDegree

/-- One measured draw: `kind` 0 = E3 (m = 2), 1 = E4/extension/Part A/B
(unsat, ladder), 2 = Part C (sat), 3 = Part D′ control. `value` is `d_F`
or `d_last`; `none` = censored above `cap`. `dstar` is the ladder's plain
`D*` where recorded (a lower bound when `dstarAtLeast`). -/
structure Row where
  kind : Nat
  m : Nat
  n : Nat
  ell : Nat
  shape : Nat     -- grouping key within a kind (cell or n)
  value : Option Nat
  cap : Nat
  dstar : Option Nat
  dstarAtLeast : Bool
  deriving Repr

def Row.v (r : Row) : Nat := r.value.getD (r.cap + 1)

/-- insertion sort: structural, so the kernel can evaluate it -/
def ins (x : Nat) : List Nat → List Nat
  | [] => [x]
  | y :: t => if x ≤ y then x :: y :: t else y :: ins x t

def isort : List Nat → List Nat
  | [] => []
  | x :: t => ins x (isort t)

/-- lower median of a list of naturals (the ⌈k/2⌉-th smallest) -/
def lowerMedian (xs : List Nat) : Nat :=
  let s := isort xs
  s.getD ((s.length + 1) / 2 - 1) 0

def keys (rs : List Row) : List Nat := (rs.map (·.shape)).eraseDups

def medianOf (rs : List Row) (k : Nat) : Nat :=
  lowerMedian ((rs.filter (·.shape == k)).map (·.v))

def e3 (rs : List Row) : List Row := rs.filter (·.kind == 0)

def p3a (rs : List Row) : Bool := (e3 rs).all fun r => decide (r.n ≤ 32 → r.v ≤ 5)

/-- medians in increasing `n`, for 18 ≤ n ≤ 32 -/
def e3Medians (rs : List Row) : List Nat :=
  let sel := (e3 rs).filter fun r => decide (18 ≤ r.n ∧ r.n ≤ 32)
  (isort (keys sel)).map (medianOf sel)

def changes : List Nat → Nat
  | a :: b :: t => (if a == b then 0 else 1) + changes (b :: t)
  | _ => 0

def nondecreasing : List Nat → Bool
  | a :: b :: t => decide (a ≤ b) && nondecreasing (b :: t)
  | _ => true

def p3b (rs : List Row) : Bool :=
  let ms := e3Medians rs
  nondecreasing ms && decide (changes ms ≤ 1) && ms.all fun x => decide (x ≤ 5)

def e4 (rs : List Row) : List Row := rs.filter fun r => r.kind == 1 && r.m == 3 && r.n ≤ 7

def p4a (rs : List Row) : Bool :=
  (e4 rs).all fun r => match r.dstar with
    | some d => decide (r.v ≤ d)
    | none => false

def p4b (rs : List Row) : Bool :=
  let strict := ((e4 rs).filter fun r => match r.dstar with
    | some d => decide (r.v < d)
    | none => false).length
  decide (2 * strict ≥ (e4 rs).length)

def boundedAtToySize (rs : List Row) : Bool :=
  (rs.filter (·.kind == 1)).all (fun r => decide (r.v ≤ 4)) &&
  (rs.filter (·.kind == 2)).all (fun r => decide (r.v ≤ 4))

def semaevSpecific (rs : List Row) : Bool :=
  let c := rs.filter (·.kind == 3)
  !c.isEmpty && (keys c).all fun k => decide (5 ≤ medianOf c k)

structure Report where
  p3a : Bool
  p3b : Bool
  p4a : Bool
  p4b : Bool
  bounded : Bool
  semaevSpecific : Bool
  deriving DecidableEq, Repr

def report (rs : List Row) : Report :=
  { p3a := p3a rs, p3b := p3b rs, p4a := p4a rs, p4b := p4b rs,
    bounded := boundedAtToySize rs, semaevSpecific := semaevSpecific rs }

-- sanity on hand-made rows
private def r0 (n v : Nat) : Row := { kind := 0, m := 2, n, ell := n / 2, shape := n, value := some v, cap := 6, dstar := none, dstarAtLeast := false }
example : p3b [r0 18 4, r0 20 4, r0 22 5, r0 24 5] = true := by decide
example : p3b [r0 18 4, r0 20 5, r0 22 4] = false := by decide
example : p3a [r0 30 6] = false := by decide
example : lowerMedian [5, 7, 6, 5] = 5 := by decide

end FallDegree
