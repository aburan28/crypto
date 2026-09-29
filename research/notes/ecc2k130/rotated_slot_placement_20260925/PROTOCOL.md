# Frozen n19 normal-coordinate slot-placement experiment

Status: pre-outcome protocol. This asks whether allocation of the one leftover
normal-basis coordinate changes *exact projected toy support* at an equal
factor and column budget. It is a new positional question after the fixed
rank-two beta sweep [#769](https://github.com/aburan28/crypto/pull/769)
and five-base union [#775](https://github.com/aburan28/crypto/pull/775), not
a repeated beta search. Neither prior experiment uses n19 six rank-three
slots with a fourth coordinate assigned to one slot.

## Frozen model and arms

Use `E: y²+xy=x³+1` over `F2[x]/(x^19+x^5+x²+x+1)`, prime subgroup
`q=130873`, normal element beta `3`, and subgroup generator `H` selected
by the merged [#767](https://github.com/aburan28/crypto/pull/767)
rule. Let `b_j=beta^(2^j)` for `j=0..18`. In every arm `i=0..5`,
slot `i` has `V_i=span(b_i,b_(i+6),b_(i+12))`. For exactly one
position `p` in `0..5`, add `b_18` to `V_p`. These six arms partition the
19 independent normal coordinates, so rank of the joint span must be 19.
Each physical factor contains **all** rational affine points whose x lies
in its slot subspace, including both signs when present. No factor point
at infinity is added.

The positions, polynomial, beta, target rule, reference, caps, code and
SHA-256 inputs are frozen before calculating any new six-sum support.
No arm can be substituted or omitted after seeing an outcome. The prior
rank-two #767 archive is an immutable source/generator cross-check, not an
equal-budget support comparator.

For each arm record every slot's rank, physical point count, number of
distinct signed `[4]P` values, signed negation-pair column count, the
joint rank, `N=product_i |F_i|`, per-slot column budget
`B=sum_i (#distinct nonzero +/-[4]P pairs)`, and the union column count
`C`. An *equal-budget pair* has identical `N` and `B`; only such pairs
can isolate placement from a changed tuple/column budget. If none exist,
the placement-isolation gate fails and no positional support conclusion
is drawn.

## Exact same-target census and relation diagnostic

Build all `q` subgroup points `kH`, `0<=k<q`, by exact group law and map
each factor's `[4]P` to its scalar. Perform exact cyclic convolution of
the six *labelled* factor choices modulo q. Store the complete `q`-entry
little-endian uint32 multiplicity array for each arm, including zero
entries, with `sum c_k=N`. Record exact support `K=#(c_k>0)`, misses
`q-K`, second moment `M2=sum c_k²`, collision pairs `(M2-N)/2`, and
the Cauchy lower bound `N²/M2<=K`. This moment already exists for the
old rank-two n19 reference in #769; here it diagnoses *new positional*
arms and is not itself claimed as a novel formula.

The same 64 distinct nonzero target scalars `k` apply to every arm:
scan counter `j=0,1,...`, take SHA-256 of
`ECC2K130-ROTATED-SLOT-PLACEMENT-20260925-v1/target/j` as a big-endian
integer, set `k=1+(digest mod(q-1))`, and retain the first 64 unseen k.
For each supported sample, save one deterministic dynamic-programming
witness, verify its exact point sum projects to `kH`, convert its signed
`[4]P` terms to canonical nonzero columns `{min(t,q-t)}`, and verify
the scalar row identity modulo q. Measure modular rank over `F_q` of
these first-witness rows for the shared 64 targets. A miss has no row.
This rank is a toy *sampled archive-witness diagnostic*, not independent
relation yield or an implicit PDP solver.

For every equal-budget pair compare its exact all-q support intersection,
symmetric difference, `K` difference, `M2`, and the shared 64-target
membership/rank. The preregistered follow-up gate is a difference of at
least 1% of q in `K` (at least 1309 points) **and** a support symmetric
difference of at least 2% of q (at least 2618 points), with both complete
arrays independently replayed. If no equal-budget pair passes, defer
slot placement as a support lever on this rung. Report all outcomes and
budget mismatches without choosing the best pair after the fact.

The necessary counting ceiling is `K<=min(q,N)` for every arm. It is
not a support prediction. No n131 transfer, solver speed, challenge log,
end-to-end `S`, or matched-rho ratio follows from this exact n19 census.

## Independent replay, caps and cost

The producer uses the #762 extended-Euclid field/curve implementation and
a sparse-support convolution with predecessor witnesses. The verifier
uses the separate #767 bit-serial/Fermat implementation. It rebuilds all
six factor lists and the full `kH` table, independently recomputes the
complete cyclic convolution using a dense-index traversal, checks every
uint32 multiplicity, all 64 memberships and witnesses, column rows and
modular ranks, and independently reapplies the exact pair gate. Source and
parent hashes fail closed. Each producer arm has a 300-second wall and
512-MiB RSS acceptance cap; each verifier arm has a 600-second wall and
512-MiB cap. The runner preserves stdout, stderr, exit status, wall,
CPU, peak RSS, source/input/output hashes, partial output, and failures.
No replacement arm or timeout exclusion is permitted.

Charge field/basis/group setup, the q-point lookup, factor lifts,
convolution, witness/rank work, serialization, and independent replay.
Report separate operation types and host-specific times; no conversion
to full IC cost is asserted. Open a draft PR with this protocol, input,
source and exact hashes before the first census. Hash-only CI must pass
before computing any new support. Afterward commit raw arrays, receipts,
verified analysis, decision, and canonical scoreboard update in the PR.
