# Full small-curve F6-IC test: preregistration

The prepared n17 comparison excludes construction of the factor-base log
table. This bounded test asks whether F6-IC can recover one previously
unseen public logarithm through the **entire native IC path**: build the
base, collect ordinary relations, verify them, build and solve the final
relation matrix, descend the target and replay the recovered scalar. It
compares F6-IC, inherited F4 and matrix F5 as full IC variants on exactly
the same point. It is a correctness and cost-attribution test, not a
security-scale or IC-versus-rho speedup claim.

Freeze Koblitz `a=1`, field degree `n=9`, factor-base recipe
`standard_subspace` dimension 4, three summands, degree 3, 8192 solver
nodes, dense final linear algebra, eight ordinary queries per batch, at
most 4096 ordinary queries, algorithm seed `20261004042`, and one new
hash-to-curve target from fixture seed `20261004041`. Generate the target
without constructing a scalar. Use the same source and binary as the
eight-target panel, one Rayon worker, no prepared state, cache disabled,
and a 180-second cap per process. Run one process per arm in fixed order
F4, F6, F5. Do not modify limits or rerun on a different target after
seeing an outcome.

Before timing, run fixture and inventory modes only to freeze the public
point, exact curve, actual distinct nonidentity subgroup-usable base
points, folded columns, and complete candidate/workload identities. Keep
unknown counts null until inventory. Commit these inputs and hashes to the
PR before the timed processes. Each IC candidate has an exact manifest
with its own `PDP` solver and split rule; the final dense `LA` stage is
separate. Keep all ordinary queries, failed PDP attempts, relation yields,
rank trajectory, final matrix, target attempts and scalar replay. A
timeout, unsupported result, incomplete rank or failed replay remains an
outcome and cannot count as a win.

Report the full inside-worker wall interval and every exclusive phase,
plus the primary one-target online interval after target-independent
preparation. Require the five online phases to sum exactly when an answer
is recovered. The cold interval starts after process launch and JSON input
loading but before curve/base construction; name that boundary rather than
equating it to a whole-process stopwatch. The F6/F4/F5 online and cold
ratios are exploratory on this unisolated Mac. Cross-method IC/rho claims
require F6 in the native `ecbench` session with the same public point and
independent replay; this protocol does not make one. No scaling inference
follows from degree nine.
