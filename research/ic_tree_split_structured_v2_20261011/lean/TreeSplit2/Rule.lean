/-!
# Round 2: the registered decision rule, as a Lean function

Round 1's rule (`../ic_tree_split_structured_v_20261011/lean/TreeSplit/Rule.lean`)
extended to the object `Z`, whose gap is taken against the *smaller* of its
two controls (`R` and `ZR`), and to the round-2 arms. Registered (by hash)
with `PREREGISTRATION.md` before any registered cell runs. `Verdict.lean`
proves, by kernel `decide`, what it returns on the recorded cells.
-/

namespace TreeSplit2

inductive Obj | t2 | z
  deriving DecidableEq, BEq, Repr

inductive VType | poly | frob | frob2 | nbw | rand | subfield
  deriving DecidableEq, BEq, Repr

/-- One registered cell. `dsZR = none` when the cell has no `ZR` arm
(`T2`), `some none` when that arm is censored. -/
structure Row where
  obj : Obj
  vtype : VType
  n : Nat
  a : Nat
  seed : Nat
  ell : Nat
  size : Nat
  dsCurve : Option Nat
  dsR : Option Nat
  dsZR : Option (Option Nat)
  deficit1 : Nat
  t2InV : Bool
  deriving Repr

inductive Verdict | alive | closed | inconclusive
  deriving DecidableEq, BEq, Repr

/-- `D_s` of the smaller control, `none` if any control arm is censored. -/
def ctl (r : Row) : Option Nat :=
  match r.dsR, r.dsZR with
  | some d, none => some d
  | some d, some (some e) => some (min d e)
  | _, _ => none

abbrev Pt := Nat × Int × Int

def points (rs : List Row) : List Pt :=
  rs.filterMap fun r =>
    match r.dsCurve, ctl r with
    | some c, some d => some (r.ell, (c : Int), (c : Int) - (d : Int))
    | _, _ => none

def slopeOk (ps : List Pt) : Bool :=
  let k : Int := ps.length
  let sl : Int := (ps.map fun p => (p.1 : Int)).foldl (· + ·) 0
  let sd : Int := (ps.map fun p => p.2.1).foldl (· + ·) 0
  let sll : Int := (ps.map fun p => (p.1 : Int) * p.1).foldl (· + ·) 0
  let sld : Int := (ps.map fun p => (p.1 : Int) * p.2.1).foldl (· + ·) 0
  decide (0 < k * sll - sl * sl) && decide (100 * (k * sld - sl * sd) ≤ 15 * (k * sll - sl * sl))

def gapAtTop (ps : List Pt) : Option Int :=
  (ps.foldl (fun acc p =>
    match acc with
    | none => some p
    | some q => if q.1 < p.1 then some p else some q) none).map (·.2.2)

def keys (rs : List Row) : List (Nat × Nat × Nat) :=
  (rs.map fun r => (r.n, r.a, r.seed)).eraseDups

def series (rs : List Row) (k : Nat × Nat × Nat) : List Pt :=
  points (rs.filter fun r => (r.n, r.a, r.seed) == k)

/-- As round 1; `slope = false` drops the slope clause (the `Z` groups). -/
def verdict (slope : Bool) (rs : List Row) : Verdict :=
  let ss := (keys rs).map (series rs) |>.filter (· ≠ [])
  let aliveGap := ss.all fun s => match gapAtTop s with
    | some g => decide (g ≤ -2)
    | none => false
  let slopes := (ss.filter slopeOk).length
  let closed := ss.all fun s => s.all fun p => decide (-1 ≤ p.2.2)
  if !ss.isEmpty && aliveGap && (!slope || 2 * slopes ≥ ss.length) then .alive
  else if closed then .closed
  else .inconclusive

def gate (rs : List Row) : Bool :=
  let sub := rs.filter (·.vtype == .subfield)
  !sub.isEmpty && sub.all fun r => r.t2InV && decide (r.n - r.ell ≤ r.deficit1)

def sel (rs : List Row) (o : Obj) (v : VType) (p : Nat → Bool) : List Row :=
  rs.filter fun r => r.obj == o && r.vtype == v && p r.n

structure Report where
  gate : Bool
  t2Frob2 : Verdict
  t2Nbw : Verdict
  t2PolyLarge : Verdict
  zPoly : Verdict
  zNbw : Verdict
  zFrob2 : Verdict
  deriving DecidableEq, Repr

/-- The registered groups (§4). -/
def report (rs : List Row) : Report :=
  { gate := gate rs
    t2Frob2 := verdict true (sel rs .t2 .frob2 (· == 21))
    t2Nbw := verdict true (sel rs .t2 .nbw fun _ => true)
    t2PolyLarge := verdict true (sel rs .t2 .poly fun n => n == 23 || n == 25)
    zPoly := verdict false (sel rs .z .poly fun _ => true)
    zNbw := verdict false (sel rs .z .nbw fun _ => true)
    zFrob2 := verdict false (sel rs .z .frob2 fun _ => true) }

private def mk (o : Obj) (ell : Nat) (c r : Option Nat) (zr : Option (Option Nat)) : Row :=
  { obj := o, vtype := .poly, n := 13, a := 0, seed := 1, ell, size := 0,
    dsCurve := c, dsR := r, dsZR := zr, deficit1 := 0, t2InV := false }

/-- `Z` gap is against the smaller control. -/
example : verdict false [mk .z 4 (some 3) (some 5) (some (some 3))] = .closed := by decide
/-- Two below both controls at the top: ALIVE without a slope clause. -/
example : verdict false [mk .z 4 (some 3) (some 5) (some (some 5))] = .alive := by decide
/-- A censored `ZR` arm excludes the cell. -/
example : verdict false [mk .z 4 (some 3) (some 5) (some none)] = .closed := by decide
/-- `T2` rows behave as in round 1. -/
example : verdict true [mk .t2 4 (some 2) (some 2) none, mk .t2 5 (some 2) (some 4) none] = .alive := by decide
example : verdict true [mk .t2 4 (some 1) (some 2) none, mk .t2 5 (some 2) (some 4) none] = .inconclusive := by decide

end TreeSplit2
