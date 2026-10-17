# Exact scalar-anchor obstruction in the P-192 Type I signed base

**Result.** In the frozen 160-point signed quotient base, the 335,712
three-sums with a nonzero unknown-log row have 333,184 distinct elliptic
targets. The 2,528 repeated triples represent the *same full formal signed
base-point vector* at each repeated target. The 357,760 known-only triples
and 2,048 zero-unknown-row triples have 349,696 distinct known-scalar target
anchors. These anchors are disjoint from all 333,184 rank-bearing targets.
Consequently, selecting queries from this exact anchor set produces no
unknown-log relation row, even though it can guarantee a raw three-sum hit.

This refines the [mixed-base census](../run-20261015T1800Z/REPORT.md):
geometric incidence, a known scalar, and a nonzero unknown-log equation are
three separate conditions. The earlier 64/64 steered hit count sampled only
anchors, and its rank 0/16 follows from the exact support intersection.

## Instance and method

The group is standard NIST P-192 over
\(p=6277101735386680763835789423207666416083908700390324961279\),
with selected prime order
\(r=6277101735386680763835789423176059013767194773182842284081\).
The 160 distinct signed quotient points consist of \(\pm B_i\) for 16
unknown-log classes and \(\pm a_j(-2G)\) for 64 guided scalars whose values
are frozen in `type_i_geom_character_collector.json`. Each class has a
verified fixed-origin degree-three representative on the genus-three
quartic in the preceding run. The present audit uses the quotient points
and the same unordered three-sum policy with repetition.

For each triple, the [native producer](../../native_support/src/bin/formal_collision.rs)
computes its exact elliptic target and its vector in 80 unsigned class
columns, with signs in the coefficients. It separates targets whose first
16 coefficients are nonzero from targets with a zero unknown block. For
every repeated target key, it compares full 80-column vectors rather than
only the point key or the 16-column unknown block. The [independent verifier](../../verify_formal_collision_counts.py)
enumerates the signed multisets without elliptic arithmetic. Equality of its
formal-vector cardinalities and the native elliptic-target cardinalities
shows that the repeated keys are exhausted by formal cancellation; the
native producer also compares the vectors directly for every collision.

| Exact count | Value |
| --- | ---: |
| Unordered signed three-sums | 695,520 |
| Nonzero-unknown-row triples / targets | 335,712 / 333,184 |
| Repeated nonzero-row triples / targets with a second formal vector | 2,528 / 0 |
| Known-only triples / targets | 357,760 / 349,696 |
| Zero-unknown-row mixed triples / distinct targets | 2,048 / 128 |
| Zero-row targets already in known-only support | 128 |
| All scalar-anchor targets / repeated anchor triples | 349,696 / 10,112 |
| Rank-bearing targets also in anchor support | **0** |

The 128 zero-row targets arise from a signed unknown pair that cancels,
leaving one known point; every one already occurs as a known-only three-sum.
The separate signed-multiset verifier reproduces all target counts and the
2,528 and 10,112 repeat counts from formal vectors alone. The exact
[machine receipt](../../type_i_formal_collision_support.json) and
[native run log](../../type_i_formal_collision_support.log) retain source,
input, group-law implementation hashes and phase timings. Those timings are
shared-host diagnostics, not an isolated speed comparison.
The [run manifest](RUN.json) binds the final source, inputs, binary, Lean
module and check outputs by SHA-256; it also confirms the four overlapping
counts and two input hashes against the preceding mixed-base receipt.

## Relation criterion and proof map

Write a rank-bearing triple as \(U+K\), where \(U\) is its unknown-base
contribution and \(K\) has a known logarithm. If its target equals an anchor
\(A\) with independently known logarithm, then
\(\log U=\log A-\log K\). If two rank-bearing triples share a target,
their unknown contributions instead satisfy
\(\log U_1-\log U_2=\log K_2-\log K_1\). The new
[Lean module](../../lean/ScalarAnchor.lean) proves these two implications as
`scalar_anchor_collision_equation` and
`repeated_target_relation_equation` for an additive-group homomorphism.
Both declarations have only the standard `propext` foundation in the
`#print axioms` audit. The finite P-192 support sets and their elliptic
point identities are checked by the native and independent enumerators, not
encoded as concrete Lean finite-field objects. The quartic divisor bridge
remains at the independently checked GP level described in the preceding
run.

This gives a practical screening criterion for a fixed base and query law:
record the rank-bearing support \(S_R\), independently known-scalar anchors
\(S_A\), the distinct formal equations per target, and the rank of any
resulting equations modulo \(r\). An anchor-conditioned collector can
obtain useful rows only from \(S_R\cap S_A\), or from a new nontrivial
collision mechanism. Here both routes have zero equations. Enlarging the
known-log base with random-looking points does not by itself solve the
scalar-anchor problem: the random-set heuristic for the present support
sizes gives \(|S_R||S_A|/r\approx1.86\cdot10^{-47}\) overlaps, while
the measured intersection is exactly zero. The heuristic is not used as a
certificate. Separately, the prior accepted-scalar certificate bounds an
ideal-uniform accepted query's fixed-base rank-bearing hit probability by
\(83874/L<2^{-173}\), where
\(L=1569275433846670190958947353987675545133315652035359680560\).

The [WeakCurves criterion registry](../../../weak_curve_leaderboard_20261012/criteria.json)
now exposes this as `ic_scalar_anchor_overlap`. It is a test of a specified
base and relation policy, not a curve-wide finding from cover geometry alone.

## Reproduction and next experiment

From the repository root, build the native producer with
`CARGO_TARGET_DIR=/private/tmp/p192-formal-cargo cargo build --release --offline --manifest-path research/p192_extension_cover_20261012/native_support/Cargo.toml --bin formal_collision`.
Run the binary with `type_i_base_points_legacy_fixture.json`,
`type_i_geom_character_collector.json`, and an output JSON path, in that
order. Then run
`python3 research/p192_extension_cover_20261012/verify_formal_collision_counts.py OUTPUT.json`.
The pinned Lean project builds with `lake build ScalarAnchor` in `lean/`;
the fresh checker and axiom audit are recorded in the run manifest.

The next collector should generate **known-scalar targets outside the
three-sum anchor set** by an explicit algebraic construction and show
intersection with nonzero-row accepted support, or discover a genuinely
distinct formal collision. Once it yields independently checked novel rows,
measure rank per ordinary query and the full one-target pipeline. The
present fixed-base candidate still has `candidate_id: null` because its
relation collector and target descent are unresolved.
