# Full-width large-prime row adapter

The main-branch large-prime collector inspected for this study uses u64 field
coordinates and relation coefficients. It cannot import the primary F2^83
curve or its 81-bit subgroup order. The integrated source contains
src/cryptanalysis/koblitz_large_prime.rs, a separate exact BigUint row
eliminator for zero, one or two residual points. It also contains a bounded
chained-S3 partial-relation producer over two disjoint signed-Frobenius
x-domains. Neither component has been measured as a retained N83
index-calculus run.

`MixedOrbitDomains::new` checks that each supplied ambient factor base has
an exact Frobenius closure and that the small and residual x-domains are
disjoint. `search` assigns three or more small summands followed by one or
two residual summands, for five or six summands in total. It constructs the
public equation [a]G + [b]Q, searches each intermediate-identity pattern
under explicit variable, domain-clause, conflict and model caps, and checks
the original field equations, per-position domain membership, every sign
choice and subgroup membership before returning a `PartialCandidate`.
An exhausted cap or invalid SAT model is inconclusive; only exhaustive
pattern completion is marked refuted. First-pair ordering is optional and
applies to the exchangeable first two small domains. The caller must still
charge unsuccessful searches and feed every accepted candidate to the
independent graph replay.

An input partial carries a unique source ID and the point equation
[a]G + [b]Q = sum(small points) + sum(residual points). The adapter checks
the equation using curve arithmetic. It currently requires the public target,
every selected factor-base point and every residual point to lie in the
prime-order subgroup. This condition holds for the retained primary imported
base, but is a real restriction on future large-prime search policies: partials
with cofactor-torsion components are not accepted by this version.

Each small point is mapped to its signed Frobenius orbit with the exact
signed lambda-to-the-phase coefficient. Each residual point is mapped to a
canonical signed Frobenius orbit; an inverse phase coefficient is independently
checked by scalar multiplication. The graph stores at most one pivot per
residual orbit and reduces rows modulo the full subgroup order. It retains
the original source IDs and their modular combination coefficients, so a
completed row can be replayed from retained partials. It checks the completed
row as a group equation before returning it. Explicit limits on pivot count,
merge depth and row size yield UnknownCap rather than a silently truncated
relation.

`solve_completed_full_rank` is the matrix bridge for completed graph rows. It
replays every input row in the group, builds the signed-orbit columns and one
target-log column over the full subgroup modulus, and requires a unit pivot
for **every** column. The generic Gaussian routine can otherwise return
arbitrary zero values for free columns, so a merely nonempty solution is not
accepted. The bridge rejects inconsistent rows and verifies each recovered
orbit logarithm and the target logarithm by scalar multiplication. Its
`WideSolvedLogs` result carries both sets of values. This bridge uses the
adapter's subgroup-only, signed-Frobenius column basis; it does not implement
the generic driver's optional projected-orbit or non-negation column modes.
It uses dense elimination; its retained-N83 memory and runtime have not been
measured and must be charged in a cold comparison.

Small-group tests close an actual two-residual signed-orbit cycle, reject a
false point equation and duplicate source ID, and check sign cancellation.
They also combine the cycle with an independent full row to check target and
orbit-log extraction, distinguish dependent rows from full rank, and reject a
tampered completed row.
A separate 81-bit modular test checks double-large-prime cycle coefficients,
source provenance and cap classifications. A Koblitz F2^7 regression builds
disjoint two-orbit domains and obtains one SAT-discovered one-residual point
equation, then independently feeds it to the graph. It also checks an exact
two-residual sign lift and graph feed, overlap rejection and fail-closed
preflight caps. The two-residual SAT search has no positive regression or
retained N83 yield measurement. These checks establish algebra and one
small-curve producer
path. They do not measure retained N83 partial-relation yield, rank gain,
graph size, factor-base advantage or complete runtime.

The next integration gate is a source-pinned retained-N83 producer run with
the small and residual base hashes fixed. It must collect and charge
unsuccessful searches and partials, retain every source equation, and pass
independent graph and group replay. It must then measure graph filtering,
matrix rank, sparse linear algebra and one complete cold target run against
matched rho. The matrix bridge above is not a retained-N83 rank measurement.
