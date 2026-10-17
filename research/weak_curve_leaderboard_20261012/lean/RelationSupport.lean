import Mathlib.Data.Sym.Card
import Mathlib.Data.Finset.Card

/-!
For an unordered positive `m`-summand policy on `B` usable base classes,
`Sym (Fin B) m` is the exact set of combinations with repetition. Any map
from combinations to target classes has image cardinality at most the
stars-and-bars count. The geometry-specific summation map and the size of a
conditioned query set are separate inputs.
-/

namespace WeakCurveLeaderboard

theorem unordered_relation_support_bound {A : Type*} [DecidableEq A]
    (B m : ℕ) (target : Sym (Fin B) m → A) :
    ((Finset.univ : Finset (Sym (Fin B) m)).image target).card ≤
      (B + m - 1).choose m := by
  classical
  calc
    ((Finset.univ : Finset (Sym (Fin B) m)).image target).card ≤
        (Finset.univ : Finset (Sym (Fin B) m)).card := Finset.card_image_le
    _ = (B + m - 1).choose m := by
      simpa using (Sym.card_sym_eq_choose (α := Fin B) m)

theorem conditioned_relation_support_bound {A : Type*} [DecidableEq A]
    (B m : ℕ) (target : Sym (Fin B) m → A) (accepted : Finset A) :
    (accepted ∩ (Finset.univ : Finset (Sym (Fin B) m)).image target).card ≤
      (B + m - 1).choose m := by
  classical
  calc
    (accepted ∩ (Finset.univ : Finset (Sym (Fin B) m)).image target).card ≤
        ((Finset.univ : Finset (Sym (Fin B) m)).image target).card :=
          Finset.card_le_card Finset.inter_subset_right
    _ ≤ (B + m - 1).choose m := unordered_relation_support_bound B m target

end WeakCurveLeaderboard
