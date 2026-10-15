import Init

set_option maxRecDepth 8192
set_option maxHeartbeats 400000

namespace FullImageSupport

-- The scalar coordinates are the frozen native census inputs. This module
-- checks their finite combinatorics, not the elliptic-curve coordinate bridge.
def sourceScalars : List Nat := [12, 48, 55, 42, 29, 33, 43, 36]
def projectedScalars : List Nat := [12, 48, 55, 42, 29, 33, 43, 36]
def directScalars : List Nat := [50, 24, 17, 18, 45, 28, 34, 38]
def fullScalars : List Nat := [52, 46, 26, 39, 7, 29, 56, 3]

def signed (base : List Nat) : List Nat :=
  base.flatMap (fun x => [x, (61 - x) % 61])

def supported (base : List Nat) (k : Nat) : Bool :=
  (signed base).any (fun x =>
    (signed base).any (fun y => (x + y) % 61 == k))

def support (base : List Nat) : List Nat :=
  (List.range 61).filter (fun k => k != 0 && supported base k)

def only (left right : List Nat) : List Nat :=
  left.filter (fun k => !(right.contains k))

def overlap (left right : List Nat) : List Nat :=
  left.filter (fun k => right.contains k)

/-- Exact support memberships for the four frozen scalar bases. -/
theorem fullSupport :
    support fullScalars =
      (List.range 61).filter (fun k => k != 0 && !([1, 5, 56, 60].contains k)) := by
  decide

theorem projectedSupport :
    support projectedScalars =
      (List.range 61).filter (fun k => k != 0 && !([2, 28, 33, 59].contains k)) := by
  decide

theorem directSupport :
    support directScalars =
      (List.range 61).filter (fun k => k != 0 && !([24, 30, 31, 37].contains k)) := by
  decide

theorem sourceEqualsProjected :
    support sourceScalars = support projectedScalars := by
  rfl

theorem fullProjectedDifference :
    only (support fullScalars) (support projectedScalars) = [2, 28, 33, 59] ∧
    only (support projectedScalars) (support fullScalars) = [1, 5, 56, 60] ∧
    (overlap (support fullScalars) (support projectedScalars)).length = 52 := by
  decide

theorem fullDirectDifference :
    only (support fullScalars) (support directScalars) = [24, 30, 31, 37] ∧
    only (support directScalars) (support fullScalars) = [1, 5, 56, 60] ∧
    (overlap (support fullScalars) (support directScalars)).length = 52 := by
  decide

#print axioms FullImageSupport.fullSupport
#print axioms FullImageSupport.projectedSupport
#print axioms FullImageSupport.directSupport
#print axioms FullImageSupport.sourceEqualsProjected
#print axioms FullImageSupport.fullProjectedDifference
#print axioms FullImageSupport.fullDirectDifference

end FullImageSupport
