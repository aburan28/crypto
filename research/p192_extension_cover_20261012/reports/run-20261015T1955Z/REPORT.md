# P-192 Type I adaptive collector: short-sum collision screen

**Result.** The 160 signed Type I quotient classes in the frozen P-192 base
have independently replayed fixed-origin degree-three Jacobian
representatives. For each of 64 accepted known-scalar queries, the adaptive
three-sum equation introduces one new quotient point. An exact pair-index
search finds that none of these 64 new points is a sum of two or four old
signed classes. Combined with the previous one- and three-sum audits, each
selected new point has **old-base sum length at least five** (allowing
repetition). Thus the tested short-sum collision route yields no equation
that eliminates its private logarithm column.

## Geometry and query law

The base consists of 16 unknown-log classes and 64 known-log guided classes,
with both signs, on the explicit smooth P-192 Type I plane quartic. Each
class has a common-origin representative \([D-T]\) with \(D,T\) effective
degree-three divisors over \(\mathbb F_p\). The [mixed-base report](../run-20261015T1800Z/REPORT.md)
identifies the [80 positive classes](../../type_i_fixed_origin_reduction.json),
[16 negative base classes](../../type_i_negative_fixed_origin_reduction.json),
and [64 negative guided classes](../../type_i_guided_negative_fixed_origin_reduction.json).
Their 320 cubic-section cuts were independently checked for residual
coordinates, curve incidence, squarefreeness, and the common-origin class
identity. All 160 numerator and all 160 denominator section-coefficient
vectors are individually distinct; this is a coefficient-sharing screen,
not an assertion about every possible section in each linear system.

The new search retains the frozen accepted-scalar law. For each
\(0\le j<64\), the [adaptive receipt](../../type_i_adaptive_private_columns.json)
gives a nonzero scalar \(s_j\), its known point \(Q_j=-2s_jG\), two old
unknown-log classes \(B_{j\bmod16}\) and \(B_{(7j+3)\bmod16}\), and the
new quotient point
\[
C_j=Q_j-B_{j\bmod16}-B_{(7j+3)\bmod16}.
\]
Both sides of \(B_{j\bmod16}+B_{(7j+3)\bmod16}+C_j=Q_j\) are replayed
on the P-192 curve. The scalar law accepts only the frozen two-character
`++` cell. The earlier mixed-base census showed 83,874 accepted nonzero-row
three-sum targets, but 64 fresh accepted scalar queries missed them. The
adaptive law makes 64 equations by adjoining the \(C_j\); the present test
asks whether those points already have a short representation in the old
geometric base.

## Exact two- and four-sum test

The [native producer](../../native_support/src/bin/adaptive_short_sums.rs)
constructs all unordered pairs of the 160 signed old points, with repeated
indices allowed. It indexes their exact P-192 sums by affine point key.
For each \(C_j\), it tests direct pair membership and every complement
\(C_j-(B_u+B_v)\) against that pair index. Any four-term representation,
with any old signs and repetitions, appears in this meet-in-the-middle
enumeration. The [independent Python verifier](../../verify_adaptive_short_sums.py)
reconstructs the points with a separate affine group law, checks the frozen
adaptive equations, repeats the pair/complement enumeration, and compares
the exact hits and counts to the [machine receipt](../../type_i_adaptive_short_sums.json).

| Frozen quantity | Exact result |
| --- | ---: |
| Signed old classes / independently checked degree-three representatives | 160 / 160 |
| Unordered old pairs with repetition / distinct pair-sum points | 12,880 / 12,801 |
| Selected accepted known-scalar queries and new points | 64 |
| New points in old base; with old three-sum representation | 0; 0 |
| New points with old two-sum representation | 0 |
| Four-sum pair-complement tests / hits | 824,320 / 0 |
| Forced rows / private new columns / old-log dimensions remaining | 64 / 64 / 16 |

The native run recorded 4,413 ms of setup and 29,136 ms of search on a shared
host. The independent replay recorded 9,678 ms. These intervals are
correctness-run diagnostics, not calibrated collector costs. This exact
screen covers the 64 selected new points and the frozen 160-point base;
the earlier 8,704-point all-pairs adaptive expansion has its own one- and
three-sum audit, but no four-sum census in this run.

## Proof map and next relation step

If an adaptive point has a second representation \(C=D+E\), then its
known-target equation becomes \(B_i+B_k+D+E=Q\). If it has
\(C=D+E+F+H\), it becomes a six-summand equation in old classes.
The pinned Lean module [ShortSums.lean](../../lean/ShortSums.lean) checks both
implications in any additive commutative group. The concrete absence of
short representations uses the two independent P-192 implementations;
the Lean statements do not encode the finite-field curve. The existing
[AdaptiveBase.lean](../../lean/AdaptiveBase.lean) proves that when each row
retains a private fresh column, every assignment of the 16 old logarithms
extends uniquely. The 64 selected rows therefore add 64 rows and 64
columns without reducing the old-log dimension.

The next collector gate is an independently occurring collision or a
geometric divisor relation that removes a fresh column, followed by the
rank of the projected old-log matrix under ordinary known-scalar queries.
Any new quartic class used there needs the same fixed-origin divisor
reduction and exact quotient-map replay. The current proposal remains
`candidate_id: null` because no full relation, descent, and one-target
pipeline has been specified.

## Reproduction

From the repository root, build through the checked offline Rust manifest:

```sh
CARGO_TARGET_DIR=/private/tmp/p192-short-sums-cargo cargo build --release --offline --manifest-path research/p192_extension_cover_20261012/native_support/Cargo.toml --bin adaptive_short_sums
/private/tmp/p192-short-sums-cargo/release/adaptive_short_sums research/p192_extension_cover_20261012/type_i_base_points_legacy_fixture.json research/p192_extension_cover_20261012/type_i_geom_character_collector.json research/p192_extension_cover_20261012/type_i_adaptive_private_columns.json research/p192_extension_cover_20261012/type_i_adaptive_short_sums.json
python3 research/p192_extension_cover_20261012/verify_adaptive_short_sums.py research/p192_extension_cover_20261012/type_i_base_points_legacy_fixture.json research/p192_extension_cover_20261012/type_i_geom_character_collector.json research/p192_extension_cover_20261012/type_i_adaptive_private_columns.json research/p192_extension_cover_20261012/type_i_adaptive_short_sums.json
```

The [native pass log](../../type_i_adaptive_short_sums.log),
[independent pass log](../../type_i_adaptive_short_sums_independent.log),
and [run manifest](RUN.json) bind the exact inputs, sources, binary, receipt,
proof checks, and statuses. Run `lake build ShortSums`,
`lake env leanchecker --fresh --verbose ShortSums`, and
`lake env lean ShortSums.lean` in the pinned Lean project to replay the
structural proof and its axiom audit.
