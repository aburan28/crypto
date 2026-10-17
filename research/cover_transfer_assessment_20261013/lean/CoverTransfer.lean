import Mathlib.GroupTheory.OrderOfElement

/-!
# Conditional injectivity of a cover pullback

This file isolates the group-theoretic step used after a cover's pullback and
pushforward maps have already been constructed and certified.  If their
composition on a finite subgroup is multiplication by the cover degree, and
that degree is coprime to the subgroup order, then pullback is injective on
that subgroup.

The theorem does not construct either map or establish the composition
identity.  Those remain per-cover proof obligations.
-/

set_option autoImplicit false

namespace CoverTransfer

variable {E J : Type*} [AddGroup E] [AddGroup J]

/--
Let `S` be a finite subgroup of the source group, of cardinality `r`.  If
`push (pull x) = d • x` for every `x ∈ S` and `d` is coprime to `r`, then
`pull` is injective after restriction to `S`.

For a degree-`d` curve cover, the intended maps are pullback into a Jacobian
and pushforward back to the elliptic curve.  Their existence, well-definedness
on divisor classes, and the displayed composition identity are hypotheses,
not consequences of this theorem.
-/
theorem pull_injective_on_subgroup
    (S : AddSubgroup E) [Finite S]
    (pull : E →+ J) (push : J →+ E)
    (d r : ℕ)
    (horder : Nat.card S = r)
    (hcoprime : Nat.Coprime d r)
    (hcomposition : ∀ x : S, push (pull (x : E)) = d • (x : E)) :
    Function.Injective (fun x : S ↦ pull (x : E)) := by
  intro x y hpull
  have hscaledE : d • (x : E) = d • (y : E) := by
    calc
      d • (x : E) = push (pull (x : E)) := (hcomposition x).symm
      _ = push (pull (y : E)) := congrArg push hpull
      _ = d • (y : E) := hcomposition y
  have hscaledS : d • x = d • y := by
    apply Subtype.ext
    exact hscaledE
  have hd : (Nat.card S).Coprime d := by
    rw [horder]
    exact hcoprime.symm
  exact (Nat.Coprime.nsmul_right_bijective (G := S) hd).1 hscaledS

end CoverTransfer
