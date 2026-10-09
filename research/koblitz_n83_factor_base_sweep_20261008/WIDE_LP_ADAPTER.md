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

Small-group tests close an actual two-residual signed-orbit cycle, reject a
false point equation and duplicate source ID, and check sign cancellation.
A separate 81-bit modular test checks double-large-prime cycle coefficients,
source provenance and cap classifications. These tests establish the adapter's
algebra on those fixtures. They do not demonstrate a partial-relation yield,
rank gain, N83 graph size, factor-base advantage or complete DLP runtime.

The next integration gate is a source-pinned partial-relation producer for the
retained primary N83 base. It must supply exact original point equations,
collect and charge unsuccessful partials, retain each input source, and pass
independent graph and group replay. Only then can graph filtering, matrix
rank and a complete one-target cold run be compared with the other sweep arms.
