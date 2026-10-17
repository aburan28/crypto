# Frozen isomorphism and six-class rho controls

Date: 2026-10-13 local. Conductor task T-2372. This protocol is committed
before implementing or executing either new arm. It extends the
[paired v4 protocol](../PROTOCOL.md) and the
[next-test decision](../NEXT_TEST_PROTOCOL.md); those records stay frozen.

## Question and boundaries

The mechanism question is whether the degree-three transport arm's relation
yield can be reproduced by the public degree-one isomorphism to the same
destination model. The exact-target references are the existing source
negation rho and a new six-class source rho. The walk-step floors in
`S = counted group additions and doublings / sqrt(r)` are `sqrt(pi/4)`
for two classes and `sqrt(pi/12)` for six. All cold setup, collection,
recovery, descent, and verification group operations are counted. Map field
work, factor-base square roots, table probes, row arithmetic, rho hashes,
and memory retain separate counters; no complete-cost or scaling verdict
can follow from this bounded test.

## Frozen inputs and arms

Use `y²=x³+b` with `(p,b,r,B,u,zeta)=(13,2,19,4,2,3)` and
`(199,3,211,12,14,92)`. Destination coefficients are 11 and 118.
Source generators are `(1,4)` and `(1,2)`. The 16 source targets are the
`source_relations` rows, indices 0 through 7, in [raw-v4.tsv](../raw-v4.tsv),
SHA-256 `995b7bcc33f5584d7b4fc6102b7ade2b3f19c7e5b39cb406b0b1c33dc7af11ad`.
Regenerate the prior target stream, assert its coordinates against a frozen
copy of those 16 coordinates, and independently compare every new row to v4.

Per target, rerun source relations, degree-three transported relations, and
source negation rho unchanged. Add `isomorphism_relations`: apply
`psi_u(x,y)=(u²x,u³y)` to the same source generator and target, solve in the
same destination model with the same first-lifted-point factor base, pair
table, relation and descent streams, 10,000-trial caps, rank recovery, and
verification. Charge construction of `u²,u³` and both point evaluations per
target, including six field multiplications in total. Keep the distinct
degree-three map evaluation counts. Every arm must verify its recovered
scalar on the original source `(G,Q)` as well as its working model.

Add `rho_six_class` on the original source `(G,Q)`, with the existing
`rho_walk_with` engine, the same per-target rho seed, 100,000-step cap, jump
and distinguished-point policy, and `RhoWalk::negation()` walk tuning.
Canonicalize over `{±P, ±sigma(P), ±sigma²(P)}` where
`sigma(x,y)=(zeta*x,y)`, using the smallest `PrimeCurve::key`.
For each field, select the correct root of `lambda²+lambda+1=0 mod r`
by checking the image of the public generator, counting the setup group
operations used. Multiply rho coefficients by `±lambda^i mod r` for the
selected representative. Count canonicalizations and the two coordinate
field multiplications per canonicalization separately. Preserve short-cycle,
look-ahead, walk, setup, and verification counters from the rho engine.

## Acceptance, costs, and stop conditions

Run the native release executable twice and require byte-identical raw output.
An independent Rust arithmetic checker, using its own point law and modular
inverse, must validate all source/transport/isomorphism images, all recovered
scalars, target equality, phase sums, map counts, and the six-class orbit and
multiplier property on every rational source point. A failed or capped arm
remains a row and blocks its paired comparison; fix it in a new versioned run.
Record source and binary hashes, host, isolation, commands, raw TSV, replay,
phase table, graph, report, and PDF in this directory. Copy performance
visuals to `benchmarks/` from the same source.

The primary mechanism comparison is per-field summed counted group operations
and relation trials for isomorphism versus degree-three transport across the
same eight targets, with means and observed ranges. An equal collection
trajectory in the same destination model supports a coordinate/generator
explanation; a difference requires inspection of the distinct generator
action before any stronger attribution. Six-class rho is a separate exact
target reference; its result is recorded regardless of direction. A later
full-cost and scaling test still requires larger registered orders, a
base-field non-isomorphic pair, at least four sizes, and independent-host
replay.
