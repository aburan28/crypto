# Adaptive Type I base: exact rank balance for 64 known-scalar queries

**Result.** The frozen 64 accepted known-scalar P-192 queries can each be
turned into a three-sum equation by adding a new quotient point. All 64 new
points are distinct and outside the existing 160-point signed base. The
resulting 64 equations have rank 64, but also introduce 64 private unknown
logarithm columns; the original 16 logarithms remain free. No new point has
a second three-sum decomposition in the existing base. Expanding the search
to **all 136 unordered pairs** of the 16 old unknown points for each of the
64 targets produces 8,704 distinct nonidentity candidate points, with no
duplicate, old-base point, or old-base three-sum representation. This exact
audit rules out counting these forced hits as useful rank growth under the
tested adaptive policies.

## Frozen construction

The subgroup is the standard P-192 prime-order group with
\(p=6277101735386680763835789423207666416083908700390324961279\) and
\(r=6277101735386680763835789423176059013767194773182842284081\).
The old base has 16 unsigned unknown-log points \(B_i\) and 64 unsigned
known-log guided points, with both signs available. The 64 query scalars
\(s_j\) are the independently generated `++` accepted scalars in the
[preceding mixed-base receipt](../../type_i_mixed_known_log_support.json).
For \(0\le j<64\), set

\[
 Q_j=-2s_jG,\qquad i_j=j\bmod16,\qquad
 k_j=(7j+3)\bmod16,\qquad C_j=Q_j-B_{i_j}-B_{k_j}.
\]

The indices are distinct for every \(j\). Each group equation
\(B_{i_j}+B_{k_j}+C_j=Q_j\) is independently replayed with a separate
Python affine-curve implementation. The [native producer](../../native_support/src/bin/adaptive_private_columns.rs)
checks the exact point and subgroup equations, the two-character query
law, distinctness, old-base membership, all 8,704 pair-policy candidates,
and all 695,520 unordered old-base
three-sums. The [independent verifier](../../verify_adaptive_private_columns.py)
recomputes all 64 scalar multiples, 64 point equations, 64 new-point
character cells, all 8,704 pair-policy points and their sorted-key digest,
and the complete old three-sum intersection without using
the native group-law implementation. Both read the same frozen input files,
whose hashes are in the [machine receipt](../../type_i_adaptive_private_columns.json).
The [run manifest](RUN.json) binds the source, inputs, native binary, two
independent outputs, and Lean checks by SHA-256.

| Measured property | Exact value |
| --- | ---: |
| Accepted known-scalar query targets | 64 |
| New quotient points / distinct points | 64 / 64 |
| Identity among new points | 0 |
| New points already in the old signed base | 0 |
| New points with an old-base three-sum decomposition | 0 |
| Character cells of new points `++`, `+-`, `-+`, `--` | 13, 19, 20, 12 |
| Forced equations / augmented unknown columns | 64 / 80 |
| Augmented row rank / dimension of remaining assignments | 64 / 16 |
| Full 64-target × 136-pair policy: candidate rows / distinct points | 8,704 / 8,704 |
| Full policy: duplicate, identity, or old-base points | 0 / 0 / 0 |
| Full policy: points in old three-sum support | 0 |
| Full policy if all points are adjoined: rows / columns / remaining dimension | 8,704 / 8,720 / 16 |

The `++` test on 13 new points is a necessary incidence screen for the
previous fixed-origin divisor construction. This run does not assert that
those 13 have been reduced to degree-three quartic representatives. The
full 8,704-point sweep is a quotient-group algebraic screen; those points
have not undergone quartic divisor reduction. The rank-balance conclusion
holds already in the quotient group and therefore
does not depend on that extra reduction.

## Theorem and criterion

Let \(x\) be *any* assignment of the 16 old logarithms. For each new row,
the fresh logarithm is determined by
\(z_j=\log Q_j-x_{i_j}-x_{k_j}\). The private-column block is the
\(64\times64\) identity matrix, so the augmented rows are independent but
place no restriction on \(x\). The [Lean module](../../lean/AdaptiveBase.lean)
formalizes the stronger statement for any additive commutative group and
any old-column row expressions: every old assignment has a unique extension
to all fresh columns. It also proves that adding these equations preserves
the existence of old solutions. The concrete P-192 point and support
instances are checked by the two arithmetic implementations, while the Lean
theorems do not encode the finite-field curve.

This yields the `ic_adaptive_base_rank_balance` screening criterion in the
[WeakCurves registry](../../../weak_curve_leaderboard_20261012/criteria.json):
count new logarithm columns before calling a generated equation useful.
Eliminate them from the matrix and measure the rank left on the original
unknowns. A second independent representation of a new point, a collision
cycle, or an independently known logarithm can remove a private column.
This frozen 64-query search and its entire 8,704-point pair expansion found
no old-base three-sum representation or collision cycle. Each potential
equation therefore has a private column, so its projected rank on the old
columns is zero. The
criterion is conditional on the stated base and query policy; other
geometric relation mechanisms require their own incidence and rank audit.

## Reproduction and next relation gate

Build with
`CARGO_TARGET_DIR=/private/tmp/p192-adaptive-cargo cargo build --release --offline --manifest-path research/p192_extension_cover_20261012/native_support/Cargo.toml --bin adaptive_private_columns`.
The native binary takes the base, character-density, guided-collector, and
mixed-known-log JSON inputs, then an output path. Run
`python3 research/p192_extension_cover_20261012/verify_adaptive_private_columns.py`
with those four inputs and the native output as its fifth argument.
The pinned Lean project builds `AdaptiveBase` and records fresh-kernel and
axiom checks in this run's manifest. The [native log](../../type_i_adaptive_private_columns.log)
and [independent log](../../type_i_adaptive_private_columns_independent.log)
retain the exact outcomes; timings on the shared host are diagnostics.

The next constructive experiment should force a **second** geometrically
valid representation of a new quotient point without defining that second
representation from the first. Such a collision would eliminate a private
column and yield an equation in the old logarithms. It then needs fixed-origin
quartic representatives, ordinary-query yield, rank after elimination,
and complete one-target cost accounting. This proposal still has
`candidate_id: null`.
