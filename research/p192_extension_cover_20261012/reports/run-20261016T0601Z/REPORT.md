# Rational residual reuse on old-anchored P-192 quartic lines

**Result.** The complete family of lines through at least one of the 141
original points and a later point of the frozen 9,925-point P-192 Type I
quartic base contains **680,474 rational non-source residual pairs**.
Their **1,360,948 outside residual points are all distinct** as quartic
points. A second implementation replays the entire old-anchored family in
reverse pivot order and reaches the same counts and zero-reuse result.
Independent PARI arithmetic verifies 128 complete quartic cuts, while a
checked-Sage control projects their 256 residual points into the elliptic
subgroup: those 256 classes are distinct and outside the existing base.

| Old-anchored line census | Exact count |
| --- | ---: |
| Frozen selected quartic points | 9,925 |
| Original pivot points | 141 |
| Unordered pairs with an original pivot | 1,389,414 |
| Presentations of the 4,892 known source quartets | 24,460 |
| Non-source lines with split residual quadratic | 680,474 |
| Non-source lines with nonsplit residual quadratic | 684,480 |
| Rational outside residual-point presentations | 1,360,948 |
| Distinct rational outside residual points | 1,360,948 |
| Repeated geometric residual-point presentations | **0** |

The source presentation count is \(5\cdot4892\): a quartet with two old
and two new selected points contributes its old--old pair and four
old--new pairs. Every non-source pair-line has two distinct residual
points outside the frozen set by the earlier complete 49,247,850-pair
collinearity and transverse-cut certificates. This run determines which
of those residual pairs split over \(\mathbf F_p\) and whether independently
chosen old-anchored cuts reuse a rational outside point.

## Exact cut and collision test

At a selected pivot \(P=(x,y)\), write the line to a second selected
point as \((x+t,y+mt)\). The producer precomputes the directional Taylor
coefficients \(c_k(m)\) of the homogeneous quartic. The two known roots
are \(t=0\) and \(t=d\). Dividing the quartic restriction by
\(t(t-d)\) leaves

\[
 A t^2+B t+C,\qquad
 A=c_4(m),\quad B=c_3(m)+dA,\quad C=c_2(m)+dB.
\]

The full transverse certificate supplies nonzero \(A\), discriminant,
and endpoint values for every selected pair. For each non-source cut,
the [producer](../../line_nonsource_residual_reuse.py) tests the
discriminant's quadratic character. On split cuts it computes both
roots, obtains their exact \((x,y)\) coordinates, checks that neither is
selected, and hashes the full 192-bit coordinate pair. It batch-inverts
both selected-point differences and the split quadratic denominators
per pivot. It also substitutes 679 spaced residual points directly into
the quartic equation.

The [reverse verifier](../../verify_line_nonsource_reuse_full.py) visits
old pivots in descending order and their later neighbors in ascending
1,024-line blocks. It independently expands the directional coefficients,
batch-inverts each block, recomputes quadratic characters and both roots
in the opposite sign order, and constructs a separate exact residual-point
set. Its source-quartet set comes from the frozen complete collinearity
records, and it checks the same 1,389,414 pair count, split/nonsplit counts,
1,360,948 distinct outside points, and zero repeats. The
[receipt audit](../../type_i_line_nonsource_reuse_verification.json)
binds the two full results to their source and input hashes.

The [PARI replay](../../verify_line_nonsource_reuse.gp) reads the first
128 frozen split cases from the full receipt. For each, it reconstructs
the entire quartic restriction and checks that the two selected and two
reported residual points are its four distinct roots. The
[projection control](../../type_i_line_nonsource_projection_old8_128.json)
inverts those 256 residuals into the cover, applies the first elliptic
quotient and frozen cofactor projection, and checks subgroup membership
and collision with the 9,925 existing quotient classes. Its first 32
inversions and quotient images also agree with the established inverse
implementation. The first eight-pivot producer receipt supplies the same
128 cases as the full run, which the audit checks entry-for-entry.

## Dangerous-curve criterion and next experiment

The useful criterion is **projected residual reuse after leaf peeling**.
A line relation containing a subgroup-point column that occurs in no
other row can be satisfied for any assignment to the older logarithms;
the conditional statement is checked by
`P192ExtensionCover.line_private_pairs_extend` in
[LineSupport.lean](../../lean/LineSupport.lean). The complete old-anchored
screen rules out geometric point reuse in this family. The 256-point
control found no quotient-column collision, while complete projection
of the 1,360,948 outside points remains the next measurement before
assigning a projected core rank. The [catalog criterion](../../../weak_curve_leaderboard_20261012/criteria.json)
now explicitly separates geometric equality from equality after cofactor
projection, since distinct quartic points can still represent the same
subgroup logarithm column.

Under a **pairwise-uniform outside-point model**, the genus-three
Hasse--Weil lower bound on \(\#C(\mathbf F_p)\) makes the birthday union
bound for 1,360,948 presentations smaller than \(2^{-152}\). Even if
every one of the 49,247,850 selected-pair lines supplied two rational
outside points, the same model gives a bound below \(2^{-139}\). These
conditional estimates explain why unforced point reuse is a poor search
policy at this base size; they do not replace the exact finite-field
census or the projected-column test. The next implementation target is a
faster exact quotient map for outside residuals, followed by subgroup
collision hashing, a leaf-core rank computation, and the established
known-scalar query law.

| Result | Proof or replay |
| --- | --- |
| Every non-source selected-pair cut has two distinct outside geometric residuals | `nonsource_quartic_line_has_two_unselected` in [FullCollinearity.lean](../../lean/FullCollinearity.lean) and the complete collinearity and transversality receipts. |
| Two-known-root factor and simple residual quadratic | `quartic_two_known_roots_factor_identity` and `residual_quadratic_root_is_simple` in [ResidualDiscriminant.lean](../../lean/ResidualDiscriminant.lean), with the concrete field screen. |
| Full old-anchored split and reuse counts | Producer, independent reverse-order complete replay, 679 direct-substitution controls, and 128 native quartic cuts. |
| Private projected columns leave old assignments free | `line_private_pairs_extend` in [LineSupport.lean](../../lean/LineSupport.lean); its applicability to the complete new family awaits the full quotient-column identity census. |
| Actual subgroup-point projection | 256 checked-Sage points, including 32 independent inverse/quotient comparisons; the complete 1,360,948-point projection is open. |

The cited Lean modules retain their pinned toolchain, axiom audits, and
fresh kernel checks in the earlier versioned runs. The geometric and
projection observations here are computational receipts with the exact
scope shown above.

Replay from the repository root:

```sh
python3 research/p192_extension_cover_20261012/line_nonsource_residual_reuse.py --old-pivots 141 --out research/p192_extension_cover_20261012/type_i_line_nonsource_reuse_old141.json
python3 research/p192_extension_cover_20261012/verify_line_nonsource_reuse_full.py
python3 research/p192_extension_cover_20261012/prepare_line_nonsource_reuse_gp.py --screen research/p192_extension_cover_20261012/type_i_line_nonsource_reuse_old141.json --out-prefix research/p192_extension_cover_20261012/type_i_line_nonsource_reuse_old141_gp_input
gp -s 128000000 -q research/p192_extension_cover_20261012/verify_line_nonsource_reuse.gp
/Volumes/SSD990/cryptanalysis/sage -python research/p192_extension_cover_20261012/line_nonsource_projection_controls.py --screen research/p192_extension_cover_20261012/type_i_line_nonsource_reuse_old8.json --out research/p192_extension_cover_20261012/type_i_line_nonsource_projection_old8_128.json
python3 research/p192_extension_cover_20261012/verify_line_nonsource_reuse_receipt.py
```

The [run manifest](RUN.json) hashes the frozen inputs, complete receipts,
native controls, source, and logs. Elapsed durations are local correctness
diagnostics under the recorded host load.
