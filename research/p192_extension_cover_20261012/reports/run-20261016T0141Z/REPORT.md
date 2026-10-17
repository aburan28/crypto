# Complete collinearity census of the frozen P-192 Type I quartic support

**Result.** Among the 9,925 frozen affine rational points on the P-192
Type I plane quartic, exactly 4,892 lines contain at least three selected
points. Each contains exactly four: two of the original 141 points and the
two residual points of one previously frozen split cut. No other line
contains three selected points. A complete 49,247,850-pair scan and an
independent reverse-order scan agree on every line and its complete
selected-point incidence set.

The field prime and selected quotient subgroup order are

\[
\begin{aligned}
p&=6277101735386680763835789423207666416083908700390324961279,\\
r&=6277101735386680763835789423176059013767194773182842284081.
\end{aligned}
\]

The input set comprises all 141 original affine quartic points and all
9,784 distinct affine residuals from the [split-line support](../run-20261015T2201Z/REPORT.md).
The previously [replayed first-quotient projection](../run-20261016T0116Z/REPORT.md)
maps those selected points to 9,925 signed-distinct subgroup points.
Collinearity here is tested on the *quartic source coordinates*, before
projection.

| Complete frozen-set incidence | Exact count |
| --- | ---: |
| Selected affine quartic points | 9,925 |
| Unordered selected-point pairs visited | 49,247,850 |
| Lines with exactly three selected points | 0 |
| Lines with exactly four selected points | 4,892 |
| Four-point lines equal to frozen original split cuts | 4,892 |
| Novel four-point lines | 0 |
| Other old/new line compositions among the detected lines | 0 |
| Pairs whose connecting line is tangent at a selected endpoint | 0 |

The [producer](../../line_full_collinearity.py) fixes each selected point
as a pivot, batch-inverts every nonzero abscissa difference, and groups
equal slopes. It handles vertical lines separately and retains a line
only when the pivot has at least two other selected points on it. Each
unordered pair is visited once, and normalized slope/intercept keys
merge presentations of the same line across pivots. The output stores
all 4,892 keys and their complete selected-point index sets. The
least-indexed point on any line with three selected points has two later
points of the same slope, so this scan cannot miss such a line. The
[independent verifier](../../verify_line_full_collinearity.py) checks
all 9,925 points against the quartic equation, reverses their order,
uses separately implemented 1,024-point inversion blocks, and compares
the entire line-key-to-index map. Both implementations visit
\(\binom{9925}{2}\) pairs.

The [selected-endpoint tangency screen](../../line_full_selected_tangencies.py)
evaluates the quartic gradient at every selected point and checks all
49,247,850 pairs for tangent incidence at either endpoint. It finds
zero. A separate [homogeneous-gradient replay](../../verify_line_full_selected_tangencies.py)
uses all three partial derivatives and searches each tangent line
against the complete point set. Ramification at a residual point is a
different case and remains outside this selected-endpoint census.

For a non-source line whose quartic cut consists of four *distinct*
points, at most two cut points belong to this frozen set, so at least
two lie outside it. The [Lean theorem](../../lean/FullCollinearity.lean)
proves that finite-set implication from an explicit complete-census
hypothesis and four-distinct-point cut hypothesis. The concrete census,
plane-quartic equation, and source-line identities are checked by the
two integer-field programs; a Lean model of the quartic and its exact
finite-field point set remains open. Residual ramification and other
multiplicity patterns need their own incidence treatment.

This screen sharpens the geometric collector criterion: a novel
transverse line through selected points cannot close within the
9,925-point base or create a one-outside-point row. Points outside the
base could be reused by *different* new cuts, so the next collector
experiment must still build the enlarged incidence graph, peel leaves,
and measure useful rank and target overlap under a frozen known-scalar
query law. The common affine quotient offset from the prior full-line
replay must remain an explicit scalar unknown unless its logarithm is
certified.

## Proof and replay map

| Statement | Formal or computational certificate | Scope |
| --- | --- | --- |
| All 49,247,850 pairs are visited and their lines grouped | Complete producer and independent reverse-order replay | Frozen 9,925-point affine set. |
| Exactly the 4,892 original quartets have at least three selected points | Exact matching line-key and index-set receipts | Finite-field coordinate computation. |
| No line through a pair is tangent at either selected endpoint | Complete gradient census and homogeneous-gradient replay | Tangency at selected endpoints only. |
| A non-source four-distinct-point cut has at least two points outside the base | `P192ExtensionCover.nonsource_quartic_line_has_two_unselected` | Abstract finite-set theorem with complete-census and four-distinct-point hypotheses. |
| Full finite-field quartic census theorem | Open concrete Lean bridge | Coordinates and polynomial checks are independently replayed. |

The checked Lean build, axiom audit, fresh kernel replay, exact source
and input digests, independent verifier, and output records are sealed
in the versioned run manifest. The local elapsed times are correctness
diagnostics under an unisolated host.
