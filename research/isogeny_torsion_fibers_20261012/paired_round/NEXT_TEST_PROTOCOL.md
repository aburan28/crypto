# Next test: isomorphism-only transport on the frozen targets

**Status:** preregistered protocol; no solver result in this file.
**Date:** 2026-10-12 local. **Dependency:** the completed
[v4 paired round](REPORT.md) and its frozen [raw records](raw-v4.tsv)
(SHA-256 `995b7bcc33f5584d7b4fc6102b7ade2b3f19c7e5b39cb406b0b1c33dc7af11ad`).

## Hypothesis and boundary

The source/degree-three destination relation-yield difference can be
reproduced by a public base-field scaling isomorphism to the same
destination model. This is a falsifiable *mechanism* test; the v4
negation-rho arm remains the exact-source-target reference, while a
six-class rho control and fully priced `S` remain separate open gates.
The derived reference floors are `sqrt(pi/4)` for negation classes and
`sqrt(pi/12)` for all six rational automorphisms on these exact curves.

## Frozen input and new arm

Keep all 16 `Q` coordinates and target indices from `raw-v4.tsv`, eight
per source pair. Keep the source generators `G=(1,4)` over `F13` and
`G=(1,2)` over `F199`, subgroup orders 19 and 211, the destination
models `b'=11` and `b'=118`, and the v4 factor-base sizes 4 and 12.
Use the public degree-one maps
`ψ_u(x,y)=(u²x,u³y)`, with `(p,b,b',u)=(13,2,11,2)` and
`(199,3,118,14)`. The independent verifier has checked both maps on
all 44,882 ordered rational-point pairs. Apply each `ψ_u` to `G` and
the exact source `Q` before the new destination solve; no target log is
supplied to a solver.

The new `isomorphism_relations` arm uses the **same destination curve**,
first-lifted-point factor base, unordered pair table, relation and
descent SplitMix64 streams, 10,000-trial caps, modular rank/recovery,
and scalar-verification rules as the v4 `transport_relations` arm.
Run fresh source, degree-three transport and negation-rho arms on the
same targets and seeds as controls. Every arm must return the same
verified scalar. Retain failures and cap exhaustion as rows.

## Cost and decision

Charge cold `u²,u³` construction and both scaling evaluations per
workload, including field multiplications, along with all relation
setup, collection, row arithmetic, descent, verification, map and rho
counters already required by [the v4 protocol](PROTOCOL.md). Preserve
the native unpriced work separately; do not promote group-operation
`S` as a complete-cost speed ratio. Repeat the run with an independent
arithmetic checker and byte-for-byte deterministic replay. Record the
exact binary, source, inputs, host, environment, commands, raw output,
phase table, graph, PDF and failures in a new versioned run directory.

The primary comparison is the sum of counted group operations and
relation trials for `isomorphism_relations` versus
`transport_relations` on the same eight targets at each order, with
exact ratios and observed ranges. If the isomorphism arm reproduces
the degree-three arm's yield and cost, retain the coordinate effect
but do not assign it to the degree-three map. If it differs, inspect
the distinct generator/target action in the same destination model,
then test it under fully priced costs. A correctness, schema,
independent-check, or replay failure blocks that pair's comparison;
preserve the failed record and fix the cause before a new run.

The separate next reference run adds charged six-class rho on these
same source targets. Scaling claims beyond these exact groups require
growing sizes, at least one non-isomorphic pair when testing an
isogeny-specific effect, and independent-host cost replay.
