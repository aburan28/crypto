/-!
# The transfer bound

The formal shape of Jao–Miller–Venkatesan ("Do all elliptic curves of the
same order have the same difficulty of discrete log?"): an isogeny
`φ : E₁ → E₂` of degree coprime to the group order restricts to an
injective homomorphism of the prime-order subgroups, and then *any* solver
for the discrete logarithm on `E₂` solves it on `E₁` at the cost of two
evaluations of `φ`: `φ(h) = φ([x]g) = [x]φ(g)`, so the logarithm of `φ(h)`
to base `φ(g)` is a logarithm of `h` to base `g`.

The solver need not be generic: it may be structural, it may read the curve
equation, the endomorphism ring, anything. It must be a `Solver`: a
deterministic `run` that is correct on every instance (`Solver.correct`),
with a self-reported cost. Beyond that, the only things used are that `φ`
is injective and preserves `one` and `mul`. The algorithm model of
`GenericInvariance.lean` is not needed here; the two files are independent,
and together they say: methods in Nechaev's model (BSGS-type) cannot see
the conductor at all (that file), and for always-correct solvers the best
cost on the two levels differs by at most the cost of the isogeny between
them (this file).

**Scope (correction, 2026-10-17, sheriff case SHF-20261017-e4d7f6, L-N2).**
An earlier version of this comment said "generic methods cannot see the
conductor at all (that file)", unqualified; in this study's vocabulary
generic methods include `rho.*` and `kangaroo.*`, which hash the point
representation and are not instances of the first file's model. This file
bounds the best cost across levels over always-correct solvers. It says
nothing about the measured `S` of a given method, and a step-capped
randomized walk is a `Solver` only once its tape is fixed and only if it
never fails. No theorem in this file is changed.

Costs are abstract naturals; `Isogeny.cost` is the price of one evaluation
in whatever unit the solver's cost is in. For the pairs of
`research/conductor_gap_cm_pairs_20261010/` the isogeny has degree `ℓ`, and
its evaluation cost is the quantity the construction makes large.
-/

namespace TransferBound

/-- A monoid's operations. No axioms are needed: `pow` is defined by
iteration and the proofs only move it through the homomorphism. -/
structure MonoidOps (α : Type) where
  one : α
  mul : α → α → α

/-- `g^k` by repeated multiplication. -/
def pow {α : Type} (ops : MonoidOps α) (g : α) : Nat → α
  | 0 => ops.one
  | k + 1 => ops.mul (pow ops g k) g

/-- An injective map preserving `one` and `mul`. The restriction of an
isogeny to a prime-order subgroup it does not kill is one. -/
structure Hom {α β : Type} (opsG : MonoidOps α) (opsH : MonoidOps β) where
  toFun : α → β
  inj : ∀ a b, toFun a = toFun b → a = b
  map_one : toFun opsG.one = opsH.one
  map_mul : ∀ a b, toFun (opsG.mul a b) = opsH.mul (toFun a) (toFun b)

variable {α β : Type} {opsG : MonoidOps α} {opsH : MonoidOps β}

theorem Hom.map_pow (φ : Hom opsG opsH) (g : α) (k : Nat) :
    φ.toFun (pow opsG g k) = pow opsH (φ.toFun g) k := by
  induction k with
  | zero => simp [pow, φ.map_one]
  | succ k ih => simp [pow, φ.map_mul, ih]

/-- A discrete-logarithm solver on `ops`: any procedure at all, returning a
cost and a candidate logarithm, correct whenever an answer exists. -/
structure Solver {α : Type} (ops : MonoidOps α) where
  run : α → α → Nat × Nat
  correct : ∀ g h x, h = pow ops g x → pow ops g (run g h).2 = h

/-- An evaluable homomorphism: a `Hom` with the cost of one evaluation. -/
structure Isogeny {α β : Type} (opsG : MonoidOps α) (opsH : MonoidOps β) extends Hom opsG opsH where
  cost : Nat

/-- **Transfer.** Solve `(g, h)` in `G` by solving `(φ g, φ h)` in `H`. -/
def transfer (φ : Isogeny opsG opsH) (S : Solver opsH) : Solver opsG where
  run g h :=
    ((S.run (φ.toFun g) (φ.toFun h)).1 + 2 * φ.cost, (S.run (φ.toFun g) (φ.toFun h)).2)
  correct := by
    intro g h x hx
    have hx' : φ.toFun h = pow opsH (φ.toFun g) x := by
      rw [hx, φ.toHom.map_pow]
    have hk := S.correct (φ.toFun g) (φ.toFun h) x hx'
    apply φ.inj
    rw [φ.toHom.map_pow]
    exact hk

/-- The transferred solver's cost is the original's plus two evaluations. -/
theorem transfer_cost (φ : Isogeny opsG opsH) (S : Solver opsH) (g h : α) :
    ((transfer φ S).run g h).1 = (S.run (φ.toFun g) (φ.toFun h)).1 + 2 * φ.cost := rfl

/-- The transferred solver returns the very answer the original found. -/
theorem transfer_answer (φ : Isogeny opsG opsH) (S : Solver opsH) (g h : α) :
    ((transfer φ S).run g h).2 = (S.run (φ.toFun g) (φ.toFun h)).2 := rfl

/-- **The bound, both ways.** With evaluable maps in both directions, every
solver on either level has a solver on the other whose cost on any instance
is the first solver's cost on the transported instance plus two
evaluations. Hence the best cost on the two levels can differ by at most
`2 · max(cost φ, cost ψ)`: the endomorphism conductor, or any other
invariant the isogeny does not preserve, is worth no more than the isogeny. -/
theorem within_isogeny_cost (φ : Isogeny opsG opsH) (ψ : Isogeny opsH opsG) :
    (∀ S : Solver opsH, ∃ S' : Solver opsG, ∀ g h,
        (S'.run g h).1 = (S.run (φ.toFun g) (φ.toFun h)).1 + 2 * φ.cost) ∧
    (∀ T : Solver opsG, ∃ T' : Solver opsH, ∀ g h,
        (T'.run g h).1 = (T.run (ψ.toFun g) (ψ.toFun h)).1 + 2 * ψ.cost) :=
  ⟨fun S => ⟨transfer φ S, fun g h => transfer_cost φ S g h⟩,
   fun T => ⟨transfer ψ T, fun g h => transfer_cost ψ T g h⟩⟩

end TransferBound
