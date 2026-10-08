# Stage 29: public scalar-blind workflow targets

Stage 29 lets the staged `ic workflow` consume a public target without constructing or supplying its scalar. A `public_hash_seed` target is derived from a domain-separated BLAKE3 transcript containing the curve identity, seed, and counter. The low field-degree bits select an abscissa, another digest bit selects a canonical lift, and public cofactor multiplication places the result in the prime-order subgroup. Infinity is rejected.

The persisted target record contains the public coordinates, seed, accepted counter, construction name, and `target_scalar_constructed: false`. It contains no expected scalar. Individual descent accepts a recovered value only if `[d]G = Q`. The signed-Frobenius rho baseline applies the same group-identity check. Existing `known_log` and `random_seed` targets still require equality to their expected scalar as an additional control.

The factor-base logarithm database remains relation-derived and group-certified; no factor-base log is known by construction. A mixed-target test runs three known-answer controls and one public scalar-blind target through selection, distributed relation collection, sparse linear algebra, descent, resume, and rho. All four complete, while the public solution records `expected: "not_constructed"`.

This capability makes a larger public unknown-scalar experiment possible. It does not itself supply a validated n=31 factor base or a completed larger run, and it does not establish scaling or SOTA.
