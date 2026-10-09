# Full-width large-prime row adapter

The main-branch large-prime collector inspected for this study uses u64 field
coordinates and relation coefficients. It cannot import the primary F2^83
curve or its 81-bit subgroup order. The isolated study branch now contains
src/cryptanalysis/koblitz_large_prime.rs, a separate exact BigUint row
eliminator for zero, one or two residual points. This is an algebra and
verification component, not a relation-discovery backend or a measured
index-calculus run.

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
source provenance and cap classifications. These tests establish the adapter's
algebra on those fixtures. They do not demonstrate a partial-relation yield,
rank gain, N83 graph size, factor-base advantage or complete DLP runtime.

The next integration gate is a source-pinned partial-relation producer for the
retained primary N83 base. It must supply exact original point equations,
collect and charge unsuccessful partials, retain each input source, and pass
independent graph and group replay. The new matrix bridge is not a retained
N83 rank measurement. Only with a producer can graph filtering, matrix rank
and a complete one-target cold run be compared with the other sweep arms.
