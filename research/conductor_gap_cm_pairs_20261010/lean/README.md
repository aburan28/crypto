# Generic invariance, in Lean 4

`GenericInvariance.lean` states and proves, with no Mathlib, a result about
generic group algorithms in Nechaev's model: the algorithm holds handles,
asks for `mul`, `inv` and equality, and its next query depends only on the
operation count, the equality answers so far and a random tape. Such an
algorithm has the same transcript, cost and answer on `(g, h)` in `G` and on
`(φ g, φ h)` in `H`, for every injective `φ` that preserves `mul` and `inv`,
every random tape and every step bound (`run_map`, `dlog_invariant`).
Isomorphic groups are the case where `φ` is bijective.

**Scope** (corrected 2026-10-16, sheriff case SHF-20261016-3caa6f, finding
L5). The theorem covers BSGS-type algorithms, whose control flow depends
only on group-operation outcomes. These cannot see the endomorphism
conductor, the j-invariant, or anything else that is not an invariant of
the group of points. It does **not** cover the Pollard rho and kangaroo
that ecbench runs. Those choose each step by hashing the point
representation (`mix(g.key(&pos))`, `src/cryptanalysis/ecbench/generic.rs:275`;
`TunedWalker::index` in `src/cryptanalysis/ic_boundary.rs`). Earlier text
listed rho and the kangaroo as covered. The sheriff case gives a toy
demonstration, not on ecbench: on `Z_1019` and its isomorphic image in
`F_2039^*`, a representation-hashing r-adding walk had a different cost on
197 of 200 transported instances, and BSGS had a different transcript on
0 of 200. For the measured rho and kangaroo, conductor-blindness is an
empirical question (`../README.md` §5), and the ℓ = 8209 rho interval
there is still unexplained.

Build: `elan` (any Lean 4 toolchain; `lean-toolchain` pins v4.34.1), then
`lake build` in this directory. Both modules compile in a few seconds.

`TransferBound.lean` is the formal shape of Jao–Miller–Venkatesan: an
evaluable injective homomorphism `φ` (an isogeny restricted to the
prime-order subgroup) with a per-evaluation cost turns *any* solver on the
target into a solver on the source at two evaluations' extra cost, returning
the same logarithm (`transfer`, `transfer_cost`, `transfer_answer`). With
maps both ways, every solver on either level has a partner on the other
within that additive cost (`within_isogeny_cost`), so the conductor, or any
invariant the isogeny does not preserve, is worth no more than the isogeny.
The solver is not assumed generic; the file is independent of the first.

**Axiom audit.** `AxiomAudit.lean` runs `#print axioms` on every main
declaration: `decide_eq_map`, `step_map`, `run_map` and `dlog_invariant`
from the first file, and `Hom.map_pow`, `transfer`, `transfer_cost`,
`transfer_answer` and `within_isogeny_cost` from the second. It fails if any
of them depends on an axiom other than `propext`, `Classical.choice` or
`Quot.sound`, including `sorryAx`. To run it locally:
`lake build && lake env lean AxiomAudit.lean`. The audit log comes from CI
job `lean` in `.github/workflows/conductor-gap-lean.yml`, kept as the run
artifact `conductor-gap-axiom-audit`. That job also refuses `sorry`,
`admit`, `native_decide` and custom axioms, and re-checks both modules
with `leanchecker --fresh`. An earlier version of this README said the
audit printed only `propext`. No audit output had been kept, so read the
CI log for the actual axiom sets. The output of a green run is now
committed as [`axiom-audit.log`](axiom-audit.log) (run 38984047520, job
117251254820, PR #2026 head e5a884227, 2026-10-17): all nine declarations
depend on `propext` only, and `leanchecker --fresh` passed on both modules
(added 2026-10-17, sheriff case SHF-20261017-e4d7f6, L-N7).

Not here, and not formalizable with today's Mathlib: endomorphism rings,
conductors, isogenies as maps of curves. The two theorems bracket the
experiment from both sides. Methods in Nechaev's model (BSGS-type) see
nothing. For always-correct solvers, the best cost on the two levels differs
by at most two ℓ-isogeny evaluations.

**Scope of the second statement** (corrected 2026-10-17, sheriff case
SHF-20261017-e4d7f6, finding L-N2). An earlier version of this paragraph
said representation-hashing walks such as ecbench's rho and kangaroo "fall
under the second statement". `TransferBound`'s `Solver` is a deterministic
`run` that is correct on every instance. A rho or kangaroo walk is such a
solver only with a fixed tape and only if it never fails, and then the
theorem bounds the best cost across levels. It says nothing about the
measured `S` of rho, kangaroo or any other given method. Their
conductor-blindness stays an empirical question (`../README.md` §5).
