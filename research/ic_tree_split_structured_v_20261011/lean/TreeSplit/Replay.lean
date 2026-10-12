import Std.Data.HashSet
import Std.Data.HashMap
import TreeSplit.Data

/-!
# Independent replay of the `poly` cells at `ℓ ≤ 5`

`Rule.lean` checks that the verdict follows from the recorded numbers. This
file checks that some of the recorded numbers are **right**: it recomputes
them from nothing but the cell's arguments, with a second implementation
written separately from the Rust harness, and compares.

What it recomputes, per cell: the field polynomial, the Koblitz curve
arithmetic, the factor base on `V = span{1, z, …, z^{ℓ−1}}`, `|T2|`, the
random control, the `splitmix64` sample streams, and the sampled cut-out
degree `D_s` of both arms. The linear algebra is deliberately different
from the harness's: the harness eliminates on the **lowest** set column;
this file keeps a basis keyed by each vector's **highest** set bit. Rank and
span membership do not depend on that choice, so agreement is a real
cross-check rather than a transcription.

Which cells: the registered `poly` cells with `ℓ ∈ {4, 5}` (twelve; the
count is checked below so the replay cannot pass vacuously). Larger cells
are within reach of the same code but take minutes of compiled Lean each,
and the protocol registered this subset.

**Trust.** `replay_ok` is closed by `native_decide`: the Lean compiler
evaluates the check and the kernel trusts the result (the proof depends on
the axiom `Lean.ofReduceBool`). That makes it an independent
re-implementation replay, not a kernel-checked proof. The verdict theorem in
`Verdict.lean` is kernel-checked (`decide`) and does not use it.
-/

namespace TreeSplit.Replay

/-! ## `splitmix64` -/

def mix (s : UInt64) : UInt64 × UInt64 :=
  let s' := s + 0x9E3779B97F4A7C15
  let z := (s' ^^^ (s' >>> 30)) * 0xBF58476D1CE4E5B9
  let z := (z ^^^ (z >>> 27)) * 0x94D049BB133111EB
  (s', z ^^^ (z >>> 31))

def streamSeed (seed n a ell tag : Nat) : UInt64 :=
  seed.toUInt64 ^^^ (a.toUInt64 <<< 40) ^^^ (n.toUInt64 <<< 32) ^^^
    (ell.toUInt64 <<< 24) ^^^ (tag.toUInt64 <<< 56)

/-! ## `GF(2)[t]` and `GF(2^n)`, on `Nat` bit patterns -/

def polyMod (a b : Nat) : Nat := Id.run do
  let db := b.log2
  let mut a := a
  for _ in [0:130] do
    if a != 0 && a.log2 ≥ db then a := a ^^^ (b <<< (a.log2 - db))
  return a

def polyGcd (a b : Nat) : Nat := Id.run do
  let mut x := a
  let mut y := b
  for _ in [0:130] do
    if y != 0 then
      let r := polyMod x y
      x := y
      y := r
  return x

structure F where
  n : Nat
  poly : Nat

def F.mul (f : F) (a b : Nat) : Nat := Id.run do
  let mut r := 0
  let mut x := a
  for i in [0:f.n] do
    if b.testBit i then r := r ^^^ x
    x := x <<< 1
    if x.testBit f.n then x := x ^^^ f.poly
  return r

def F.sqr (f : F) (a : Nat) : Nat := f.mul a a

def F.inv (f : F) (a : Nat) : Nat := Id.run do
  let e := 2 ^ f.n - 2
  let mut r := 1
  let mut base := a
  for i in [0:f.n] do
    if e.testBit i then r := f.mul r base
    base := f.sqr base
  return r

def F.trace (f : F) (a : Nat) : Nat := Id.run do
  let mut t := 0
  let mut x := a
  for _ in [0:f.n] do
    t := t ^^^ x
    x := f.sqr x
  return t

def F.halfTrace (f : F) (c : Nat) : Nat := Id.run do
  let mut z := 0
  let mut x := c
  for _ in [0:(f.n + 1) / 2] do
    z := z ^^^ x
    x := f.sqr (f.sqr x)
  return z

def irreducible (p n : Nat) : Bool := Id.run do
  let f : F := ⟨n, p⟩
  let mut t := 2
  for _ in [0:n / 2] do
    t := f.sqr t
    if polyGcd p (t ^^^ 2) != 1 then return false
  return true

/-- First irreducible trinomial, else first pentanomial in lexicographic
order (the convention of both tree-split studies). -/
def fieldPoly (n : Nat) : Nat := Id.run do
  for k in [1:n] do
    let p := 2 ^ n + 2 ^ k + 1
    if irreducible p n then return p
  for a in [1:n] do
    for b in [a + 1:n] do
      for c in [b + 1:n] do
        let p := 2 ^ n + 2 ^ a + 2 ^ b + 2 ^ c + 1
        if irreducible p n then return p
  return 0

/-! ## `K_a : y² + xy = x³ + a·x² + 1` -/

def pointsWithX (f : F) (a x : Nat) : List (Nat × Nat) :=
  if x == 0 then [] else
  let c := x ^^^ a ^^^ f.inv (f.sqr x)
  if f.trace c != 0 then [] else
  let y := f.mul x (f.halfTrace c)
  [(x, y), (x, y ^^^ x)]

def add (f : F) (a : Nat) : (Nat × Nat) → (Nat × Nat) → Option (Nat × Nat)
  | (x1, y1), (x2, y2) =>
    if x1 == x2 then
      if y1 ^^^ y2 == x1 || x1 == 0 then none else
      let l := x1 ^^^ f.mul y1 (f.inv x1)
      let x3 := f.sqr l ^^^ l ^^^ a
      some (x3, f.sqr x1 ^^^ f.mul (l ^^^ 1) x3)
    else
      let l := f.mul (y1 ^^^ y2) (f.inv (x1 ^^^ x2))
      let x3 := f.sqr l ^^^ l ^^^ x1 ^^^ x2 ^^^ a
      some (x3, f.mul l (x1 ^^^ x3) ^^^ x3 ^^^ y1)

/-- `T2` on the polynomial-basis subspace `V = [0, 2^ℓ)`. -/
def t2Poly (f : F) (a ell : Nat) : Std.HashSet Nat := Id.run do
  let base := (List.range (2 ^ ell)).flatMap (pointsWithX f a)
  let arr := base.toArray
  let mut s : Std.HashSet Nat := {}
  for i in [0:arr.size] do
    for j in [i:arr.size] do
      match add f a arr[i]! arr[j]! with
      | some (x, _) => s := s.insert x
      | none => pure ()
  return s

/-! ## Sets, samples, and the cut-out screen -/

def control (n size : Nat) (st : UInt64) : Std.HashSet Nat := Id.run do
  let mask : UInt64 := (1 <<< n.toUInt64) - 1
  let mut s : Std.HashSet Nat := {}
  let mut st := st
  -- at most this many draws are ever needed at these sizes; the bound only
  -- keeps the loop structural
  for _ in [0:size * 64 + 1000] do
    if s.size < size then
      let (st', z) := mix st
      st := st'
      let x := (z &&& mask).toNat
      if x != 0 then s := s.insert x
  return s

def samples (n k : Nat) (set : Std.HashSet Nat) (st : UInt64) : Array Nat := Id.run do
  let mask : UInt64 := (1 <<< n.toUInt64) - 1
  let mut out : Array Nat := #[]
  let mut st := st
  for _ in [0:k * 64] do
    if out.size < k then
      let (st', z) := mix st
      st := st'
      let x := (z &&& mask).toNat
      if !set.contains x then out := out.push x
  return out

/-- Monomials of degree exactly `k` in the variables `vs`, as bitmasks. -/
def combos : List Nat → Nat → List Nat
  | _, 0 => [0]
  | [], _ + 1 => []
  | v :: vs, k + 1 => (combos vs k).map (· ||| (1 <<< v)) ++ combos vs (k + 1)

/-- Multilinear monomials of degree `≤ d` in `n` variables, as bitmasks. -/
def monomials (n d : Nat) : Array Nat :=
  ((List.range (d + 1)).flatMap (combos (List.range n))).toArray

def evalRow (monos : Array Nat) (p : Nat) : Nat := Id.run do
  let mut row := 0
  for j in [0:monos.size] do
    let m := monos[j]!
    if p &&& m == m then row := row ||| (1 <<< j)
  return row

/-- Reduce `v` by a basis keyed by highest set bit. -/
def reduce (basis : Std.HashMap Nat Nat) (v : Nat) (fuel : Nat) : Nat := Id.run do
  let mut v := v
  for _ in [0:fuel] do
    if v != 0 then
      match basis.get? v.log2 with
      | some b => v := v ^^^ b
      | none => return v
  return v

/-- `some D_s`, or `none` if censored (no `D ≤ 6` within 90,000 columns). -/
def cutOut (pts : Array Nat) (n : Nat) (smp : Array Nat) : Option Nat := Id.run do
  for d in [1:7] do
    let monos := monomials n d
    let m := monos.size
    if m > 90000 then return none
    let mut basis : Std.HashMap Nat Nat := {}
    for p in pts do
      let v := reduce basis (evalRow monos p) (m + 1)
      if v != 0 then basis := basis.insert v.log2 v
    let fz := if basis.size == m then smp.size else
      (smp.filter fun p => reduce basis (evalRow monos p) (m + 1) == 0).size
    if fz == 0 then return some d
  return none

/-- Recompute `(|T2|, D_s(curve), D_s(R))` for a `poly` cell. -/
def recompute (r : Row) : Nat × Option Nat × Option Nat :=
  let f : F := ⟨r.n, fieldPoly r.n⟩
  let t2 := t2Poly f r.a r.ell
  let ctl := control r.n t2.size (streamSeed r.seed r.n r.a r.ell 2)
  let sc := samples r.n 2000 t2 (streamSeed r.seed r.n r.a r.ell 3)
  let sr := samples r.n 2000 ctl (streamSeed r.seed r.n r.a r.ell 4)
  (t2.size, cutOut t2.toArray r.n sc, cutOut ctl.toArray r.n sr)

def targets : List Row := cells.filter fun r => r.vtype == .poly && r.ell ≤ 5

def replayCell (r : Row) : Bool := recompute r == (r.t2, r.dsCurve, r.dsR)

/-- The replay is not vacuous: twelve registered cells are in scope. -/
theorem targets_count : targets.length = 12 := by decide

/-- Every in-scope recorded cell is reproduced exactly by this independent
implementation. See the trust note above. -/
theorem replay_ok : targets.all replayCell = true := by native_decide

end TreeSplit.Replay
