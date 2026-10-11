# Next F6-IC gate at n=83: compositional high arity

The completed arithmetic component makes pair sums and residual lookups
faster on the actual K_0 confidence-gate curve. It does **not** make a
four-summand index-calculus attack viable. At the actual dimension-12 base,
`B=4,054`, four-summand coverage of a uniformly drawn K_0 target is at
most `4.662e-12`. The counting ceiling first exceeds 1% at eight
summands; nine gives more combinatorial capacity, but no guarantee of
ordinary yield. This bound fixes the next mathematical target.

The existing Boolean F4/F5/F6 path stores monomials in a `u64`. The separate
`wide_groebner` path explicitly rejects field degree above 64 and caps
monomials at 128 variables. A direct degree-83 chained Semaev system with
the actual 12-dimensional base needs 594 variables at eight summands or
689 at nine. Enumerating the 8.2 million pair sums is feasible as a
terminal index; enumerating all triple sums or recursively scanning that
index for every possible earlier summand is not a viable substitute for
elimination.

## Admission task for the next solver

Design a compositional F6-IC solver that keeps each algebraic block within
the available monomial width, passes only **proved** constraints or exact
group residuals between blocks, and uses the new pair index as a terminal
closure. It must carry explicit `found`, `refuted`, `exhausted`, and
`unsupported` outcomes; a spent budget cannot be called a proof of no
decomposition. No performance claim is attached to this design task.

Before any n=83 timing, pin source and a frozen natural-target input set.
First prove soundness and completeness on small curves against exact
enumeration, including absent relations, sign choices, infinity, repeated
points, and cofactor-projected bases. Then run planted n=83 K_0 controls
only for correctness and measured ordinary n=83 queries for yield, node
counts, timeouts, memory, and cost per useful row. A natural query that
does not finish remains a failure row. Compare complete one-target online
IC only after a full candidate manifest, reusable log preparation, target
descent, and scalar replay exist. Pair the same public point with a
one-target rho solve under the same resources; keep stage timing and any
unisolated ratios exploratory. F4 and F5 currently do not support this
n=83 instance, so a claimed F6/F4/F5 wall-time win here would be invalid.

The first performance gate for the compositional solver is whether it can
produce and independently replay **one ordinary verified relation** on the
K_0 base within a frozen resource limit while charging every failed
attempt. Only then should we optimize its block reductions or widen the
factor base. This gate is deliberately smaller than a full DLP yet directly
tests the bottleneck that the four-summand kernel cannot address.
