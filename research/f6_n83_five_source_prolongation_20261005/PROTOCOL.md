# n83 five-summand selective degree-four source prolongation

Registered before changing code or measuring the candidate. Parent is the
exact dimension-18 five-summand system in #1417 at
`e0006066706b3132b24955dc5ba2cb0f24d5015b`. Its own-degree
reduction has 332 rank, 299,558 columns and no source-only row for any
of the four ordinary T001 lifts. This gate asks whether multiplying by
a *small, frozen set of source-coordinate variables* yields a useful
linear consequence without materializing the full degree-four Macaulay
matrix. It is a structural F6-IC experiment, not a complete solver or
speed comparison.

Use registered K0 curve `icv1-f2m83-tm6151469093347-debefd74`, its
exact subgroup, the standard dimension-18 source subspace and the same
public T001 point
`(355fb5df7a905f16921eb,5900a390f42d290f1bbe)`.
Actual factor-base counts from #1417 are 261,444 subgroup-usable points
and 130,722 sign-folded columns; these do not change. For the ordinary
gate use torsion offset zero and the exact source preimage from #1417.
For the planted control use `[0,2,4,6,8]` and its full-group sum.

The frozen source variable sequence is
`[0,18,36,54,72,1,19,37,55,73,2,20,38,56,74,3]`.
For each `k` in `[1,4,8,16]`, append to the 332 original Boolean
equations one product of each original equation with each of the first
`k` source variables, using Boolean multiplication (`x²=x`) and exact
XOR cancellation. The original equations remain. Never multiply newly
added equations again. Record exact rows, term occurrences, distinct
columns, rank, linear equations and source-only status. First check that
every added equation vanishes under the planted full assignment. Then
reduce the ordinary system with the existing deterministic degree-first
monomial order and exact GF(2) row elimination.

Stop increasing `k` after a 1,500,000-column cap, a 120-second process
limit, or observed RSS above 7 GiB; retain the stopped row and cause.
Cap and timeout results are inconclusive, not contradictions. The
preregistered structural success is at least one exact contradiction
or nonconstant affine equation involving *only source-coordinate bits*
for ordinary T001 at a completed `k`, while the planted control passes.
If this occurs, repeat at all four torsion offsets before promoting the
constraint to a solver proposal. If no such row appears at `k=16`, stop
this variable set. Do not infer natural relation yield or absence of a
decomposition from any root-only miss.

Use the native Rust implementation, release build, one thread and a
frozen binary. Save source/input/binary hashes, raw stdout/stderr,
exit codes, elapsed and sampled peak RSS. This Mac lacks an auditable
exclusive CPU partition, so wall timings are feasibility diagnostics.
The primary IC metric would still require one previously unseen target
solved and verified by complete IC, then matched same-point rho under
the repository measurement contract. Those online times and any
speedup remain unknown for this experiment.
