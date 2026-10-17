# Transverse line cuts across the complete frozen P-192 quartic support

**Result.** Every line through two of the 9,925 frozen affine rational
points on the P-192 Type I plane quartic has four distinct affine
intersection points over an algebraic closure. The complete
49,247,850-pair screen finds no vertical selected pair, no loss of
quartic degree, no repeated residual root, and no tangency at either
selected endpoint. Together with the [complete collinearity census](../run-20261016T0141Z/REPORT.md),
this means that a non-source selected-pair line has two distinct
residual points outside the frozen set. Their rationality and reuse
across different lines remain separate questions.

The field prime and selected first-quotient subgroup order are

\[
\begin{aligned}
p&=6277101735386680763835789423207666416083908700390324961279,\\
r&=6277101735386680763835789423176059013767194773182842284081.
\end{aligned}
\]

| Exact full-pair screen | Count |
| --- | ---: |
| Frozen affine quartic points | 9,925 |
| Unordered selected-point pairs | 49,247,850 |
| Vertical pairs | 0 |
| Restricted quartic of degree below four | 0 |
| Repeated roots of the residual quadratic | 0 |
| Tangencies at either selected endpoint | 0 |

For a pivot \(P=(x,y)\), a nonvertical pair line has slope \(m\),
second-point abscissa difference \(d\ne0\), and restricted quartic

\[
 F(x+t,y+mt,1)=c_1(m)t+c_2(m)t^2+c_3(m)t^3+c_4(m)t^4.
\]

Both selected points are roots, so division by \(t(t-d)\) leaves
\(A t^2+B t+C\), where

\[
 A=c_4(m),\qquad B=c_3(m)+dA,\qquad C=c_2(m)+dB.
\]

The [full producer](../../line_full_residual_discriminant.py)
precomputes each \(c_k(m)\) as a degree-\(k\) polynomial in the slope
at each pivot. It batch-inverts the abscissa differences, visits every
unordered pair exactly once, checks \(c_1+dC=0\), and tests \(A\),
\(B^2-4AC\), \(C\), and \(Ad^2+Bd+C\). These respectively detect a
degree drop, repeated residual root, tangent at the pivot, and tangent
at the second selected point. The [independent verifier](../../verify_line_full_residual_discriminant.py)
uses reversed point order and separate 1,024-point inversion blocks.
It additionally compares Taylor coefficients against direct quartic
substitution on every one of the 11,040 frozen second-generation
control pairs.

The [Lean module](../../lean/ResidualDiscriminant.lean) proves the
two-known-roots factor identity over a commutative ring. It also proves
that a quadratic root cannot be repeated when its discriminant is
nonzero, and that nonzero constant and endpoint values keep residual
roots away from the two selected roots. The concrete finite-field
quartic, exact points, and all-pair census are checked by coordinate
programs; they are not yet represented as one Lean theorem.

This transversality screen makes the earlier closure criterion exact
for the declared point set: aside from the 4,892 original split
quartets, no selected-pair line closes inside the 9,925-point base or
has only one distinct residual point. A residual quadratic with
nonzero discriminant can still be nonsplit over \(\mathbf F_p\); the
next collector experiment must classify usable rational residuals,
hash their projected subgroup classes across independent cuts, peel
fresh columns, and measure known-scalar target overlap and useful rank.

## Proof and replay map

| Result | Certificate | Boundary |
| --- | --- | --- |
| All pair-line residual quadratics have degree two and nonzero discriminant | Complete producer and independent reverse-order replay | 49,247,850 frozen affine point pairs. |
| Taylor quotient agrees with direct quartic restriction | Independent 11,040-pair substitution control | Previously frozen second-generation pair selection. |
| Two-known-roots quotient and simple residual roots | `quartic_two_known_roots_factor_identity`, `residual_quadratic_root_is_simple`, and `residual_root_avoids_known_points` | General commutative-ring algebra with explicit nonzero hypotheses. |
| Concrete finite-field all-pair statement | Open Lean bridge | Coordinate and subgroup inputs are frozen and independently replayed. |

The versioned run manifest retains the inputs, exact source, complete
receipts, replay logs, Lean build, axiom audit, fresh kernel check, and
artifact hashes. Local elapsed time is an unisolated correctness
diagnostic.
