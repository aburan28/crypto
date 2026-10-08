# n83 F6 wide search: derive valid coordinates from the factor base

Registered on 2026-10-05 before implementing or measuring this change.
`wide_groebner::valid_coordinates` currently loops over every coordinate
in the subspace, constructs its field element, solves the curve equation
for its point lifts, then asks whether a lift is in `index_of`. The
factor base already contains those exact curve points. The hypothesis
is that forming the set of allowed abscissae directly from indexed
factor-base points, and scanning coordinates with field-bit XORs,
returns the identical sorted coordinate list while avoiding thousands
of redundant curve lifts per query.

The baseline is commit `128cf786a60da4f4289d8680092916991cf90f8a`
(the n83 lazy-child PR), with frozen root probe binary SHA-256
`16bfcd5408a8c736269b9a2d1359375f8bd3eb2880c194bf11a0333dfd7db22f`.
Freeze the exact K0 n83 curve, standard dimension-12 source base, planted
source indices `[0,2,4]`, two-word table, default wide-search options,
node budget one, and release profile from the preceding protocol.
The source has 4,057 points and 2,029 distinct valid abscissae.

Correctness gates before performance: compare the new coordinate list
to the old curve-lift reference at n=9, 17 and 83, including a strict
subset of indexed points and a nonstandard basis if available. Run all
wide-backend release tests. The baseline and candidate probes must
both return `exhausted` with no witness, identical nodes, reductions
and matrix dimensions; the existing 128-node planted witness probe
must still find and replay an exact source and projected-group sum.
If `index_of` does not describe the supplied factor base, preserve the
old behavior rather than silently accepting external or invalid points.

Performance gate: freeze the candidate native binary; run five
baseline/candidate pairs in alternating order, each with a 120-second
timeout. Record every JSON line, status, stderr, setup/query interval
and peak RSS. Require at least 1.5× lower median query time, no more
than 1.25× candidate peak RSS, and identical outcomes. Otherwise
revert the runtime change and retain the failure. The measured root
is a planted three-summand block, not a complete eight-summand F6
decomposition or ordinary relation. The physical arm64 Mac lacks the
CPU-isolation receipt, so local ratios are exploratory stage
diagnostics; full F6/IC speedups and natural relation yield stay
unknown regardless of the result.
