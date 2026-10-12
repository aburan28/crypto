/-!
# The registered decision rule, as a Lean function

`PREREGISTRATION.md` §4 states the rule in prose. This file states the same
rule as a total, computable function, and is registered (by hash) with it,
before any registered cell runs. After the run, `Verdict.lean` proves, by
kernel evaluation (`decide`), what this function returns on the recorded
cells. So the verdict in `RESULTS.md` cannot differ from what the
registered rule says about the recorded data: a typo, a transcription slip
or a generous reading of the rule fails the build.

What this does **not** check: that the recorded `D_s` values are right.
That is `Replay.lean`'s job (an independent recomputation of a subset of
cells), and the Rust harness's self-test. Lean has no clock, so nothing here
is a statement about time; the metric is a degree, and the study is untimed.

No Mathlib: integers, lists and booleans only, so the kernel can evaluate it.
-/

namespace TreeSplit

/-- The factor-base subspace families of the study. -/
inductive VType
  | poly | frob | rand | subfield
  deriving DecidableEq, BEq, Repr

/-- One registered cell, as the harness records it. `none` is a censored
sampled cut-out degree. `deficit1` is the curve arm's degree-1 deficit, and
`t2InV` whether every element of `T2` lies in `V` (recorded for `subfield`
only, `false` elsewhere). -/
structure Row where
  vtype : VType
  n : Nat
  a : Nat
  seed : Nat
  ell : Nat
  t2 : Nat
  dsCurve : Option Nat
  dsR : Option Nat
  deficit1 : Nat
  t2InV : Bool
  deriving Repr

inductive Verdict
  | alive | closed | inconclusive
  deriving DecidableEq, BEq, Repr

/-- A point of a series: `(ℓ, D_s(curve), gap)`, uncensored cells only. -/
abbrev Pt := Nat × Int × Int

def points (rs : List Row) : List Pt :=
  rs.filterMap fun r =>
    match r.dsCurve, r.dsR with
    | some c, some d => some (r.ell, (c : Int), (c : Int) - (d : Int))
    | _, _ => none

/-- OLS slope of `D_s(curve)` on `ℓ` is `≤ 0.15`, cross-multiplied so it
stays in the integers: `slope = num / den` with `den > 0`. -/
def slopeOk (ps : List Pt) : Bool :=
  let k : Int := ps.length
  let sl : Int := (ps.map fun p => (p.1 : Int)).foldl (· + ·) 0
  let sd : Int := (ps.map fun p => p.2.1).foldl (· + ·) 0
  let sll : Int := (ps.map fun p => (p.1 : Int) * p.1).foldl (· + ·) 0
  let sld : Int := (ps.map fun p => (p.1 : Int) * p.2.1).foldl (· + ·) 0
  let num := k * sld - sl * sd
  let den := k * sll - sl * sl
  decide (0 < den) && decide (100 * num ≤ 15 * den)

/-- The gap at the largest uncensored `ℓ`. -/
def gapAtTop (ps : List Pt) : Option Int :=
  (ps.foldl (fun acc p =>
    match acc with
    | none => some p
    | some q => if q.1 < p.1 then some p else some q) none).map (·.2.2)

def keys (rs : List Row) : List (Nat × Nat × Nat) :=
  (rs.map fun r => (r.n, r.a, r.seed)).eraseDups

def series (rs : List Row) (k : Nat × Nat × Nat) : List Pt :=
  points (rs.filter fun r => (r.n, r.a, r.seed) == k)

/-- §4: **ALIVE** if the gap at the largest uncensored `ℓ` is `≤ −2` in
every series and at least half the series (rounded up) have slope `≤ 0.15`
(the lower median); **CLOSED** if every uncensored gap is `≥ −1`;
**INCONCLUSIVE** otherwise. Series with no uncensored cell are dropped. -/
def verdict (rs : List Row) : Verdict :=
  let ss := (keys rs).map (series rs) |>.filter (· ≠ [])
  let aliveGap := ss.all fun s => match gapAtTop s with
    | some g => decide (g ≤ -2)
    | none => false
  let slopes := (ss.filter slopeOk).length
  let closed := ss.all fun s => s.all fun p => decide (-1 ≤ p.2.2)
  if !ss.isEmpty && aliveGap && 2 * slopes ≥ ss.length then .alive
  else if closed then .closed
  else .inconclusive

/-- §4 validity gate: every `subfield` cell has `T2 ⊆ V` and a degree-1
deficit of at least `n − ℓ`, and there is at least one such cell. A run
that fails it is invalid whatever the verdicts say. -/
def gate (rs : List Row) : Bool :=
  let sub := rs.filter (·.vtype == .subfield)
  !sub.isEmpty && sub.all fun r => r.t2InV && decide (r.n - r.ell ≤ r.deficit1)

/-- The primary arms. `frob` counts only its `n = 21` series (the only `n`
with a Frobenius-stable `V` at every registered `ℓ`); the `n = 17` cells are
supplementary. -/
def arm (rs : List Row) : VType → List Row
  | .frob => rs.filter fun r => r.vtype == .frob && r.n == 21
  | v => rs.filter (·.vtype == v)

structure Report where
  gate : Bool
  poly : Verdict
  frob : Verdict
  rand : Verdict
  deriving DecidableEq, Repr

def report (rs : List Row) : Report :=
  { gate := gate rs
    poly := verdict (arm rs .poly)
    frob := verdict (arm rs .frob)
    rand := verdict (arm rs .rand) }

/-! ## Sanity: the rule on hand-made tables

These run before any data exists and pin the rule's behaviour on its
boundary cases. -/

private def mk (v : VType) (n a seed ell : Nat) (c d : Option Nat) : Row :=
  { vtype := v, n, a, seed, ell, t2 := 0, dsCurve := c, dsR := d, deficit1 := 0, t2InV := false }

/-- Equal to the controls everywhere, rising with `ℓ`: CLOSED. -/
example : verdict [mk .poly 17 0 1 4 (some 2) (some 2), mk .poly 17 0 1 5 (some 3) (some 3)] = .closed := by decide
/-- Two below the control at the top, flat: ALIVE. -/
example : verdict [mk .poly 17 0 1 4 (some 2) (some 2), mk .poly 17 0 1 5 (some 2) (some 4)] = .alive := by decide
/-- Two below, but climbing at slope 1: INCONCLUSIVE (fails CLOSED and the slope clause). -/
example : verdict [mk .poly 17 0 1 4 (some 1) (some 2), mk .poly 17 0 1 5 (some 2) (some 4)] = .inconclusive := by decide
/-- Censored cells are excluded, not counted. -/
example : verdict [mk .poly 17 0 1 4 (some 2) (some 2), mk .poly 17 0 1 5 none (some 4)] = .closed := by decide
/-- One series alive, one closed: not ALIVE (gap clause fails), not CLOSED. -/
example : verdict [mk .poly 17 0 1 5 (some 2) (some 4), mk .poly 19 0 1 5 (some 4) (some 4)] = .inconclusive := by decide
/-- Above the controls counts as closed: the lever needs a *lower* degree. -/
example : verdict [mk .poly 17 0 1 5 (some 5) (some 3)] = .closed := by decide

end TreeSplit
