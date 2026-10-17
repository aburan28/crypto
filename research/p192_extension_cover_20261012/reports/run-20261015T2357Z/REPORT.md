# Second-generation line closure on the P-192 Type I quotient base

**Result.** Lines through at least one of the first 64 new quartic points
and another point in the 141-old-plus-64-new selection give 11,040
distinct pair queries. Of their residual quadratics, 5,636 split over
\(\mathbf F_p\) and 5,404 do not. Exactly 160 split cuts have both
residual points in the frozen 9,925-point set. All 160 are the five
alternative pair presentations of one of the 32 original four-point line
quartets. There is no new closed quartet and no case with exactly one
residual in the frozen set in this declared selection.

The field prime is
\(p=6277101735386680763835789423207666416083908700390324961279\).
The projected subgroup has prime order
\(r=6277101735386680763835789423176059013767194773182842284081\).
The run uses the same exact plane-quartic model, quotient map, anchor,
and point encoding as the [line-rank certificate](../run-20261015T2346Z/REPORT.md).

The [input generator](../../prepare_line_second_gen_input.py) binds the
full frozen coordinates and 32 source quartets to the GP input. The
[PARI/GP screen](../../line_second_gen_closure.gp) uses the exact
plane-quartic equation and frozen coordinates. It verifies polynomial
division by the two selected point abscissae, classifies every quadratic,
and searches both residual points in the full 9,925-point coordinate set.
The [independent Python replay](../../verify_line_second_gen_closure.py)
rebuilds the quartic from its JSON coefficients, substitutes each line,
divides over the prime field, and computes square roots using the
\(p\equiv3\pmod4\) formula. It reproduces every count and the empty
novel-witness set. Input and source digests are frozen in the run receipt.

| Screen outcome | Count |
| --- | ---: |
| Pairs containing a selected new point | 11,040 |
| Split residual quadratics | 5,636 |
| Nonsplit residual quadratics | 5,404 |
| Degree or repeated-root exceptions | 0 |
| Exactly one residual in frozen set | 0 |
| Both residuals in frozen set | 160 |
| Repeat of one of 32 original quartets | 160 |
| New closed quartet | **0** |

Every repeated quartet has the same four point columns and the common line
constant established by the [private-column analysis](../run-20261015T2346Z/REPORT.md).
The repeated pair presentation therefore adds no new row-space direction.
The [Lean theorem](../../lean/LineClosure.lean) formalizes row-span
containment when each new row is equal to an old row, over an arbitrary
semiring and module. The concrete row equality uses the complete 160-cut
coordinate census and the canonical line-section identity; a Lean model
of the plane quartic and quotient map is still open.

This test covers pair lines whose first endpoint is among the first 64
new points, with the second endpoint among 141 old plus those 64 new
points. It does not screen pairs involving later new points as endpoints
or other cover families. The next targeted search should select endpoints
from different independently generated families and measure repeated
projected points and core rank before investing in random target queries.

## Proof map

| Statement | Lean declaration | Finite-field bridge |
| --- | --- | --- |
| Duplicate rows stay in the old row span | `P192ExtensionCover.duplicate_relation_rows_span_le` | All 160 closed cuts have one of 32 frozen quartets. |
| The second-generation selection adds zero novel closed quartets | Open concrete finite-field bridge | Complete GP census and independent Python prime-field replay. |

The checked Lean build, axiom audit, and fresh kernel replay are stored with
the source and run manifest. A redundant first fresh-check process was
cancelled after a second audit command was launched; the retained second
process completed with exit zero. Local elapsed time is a correctness
diagnostic.
