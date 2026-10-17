# P-192 Type I: all-character reduction of 141 quotient classes

**Result.** All 141 line-derived P-192 quotient classes now have explicit,
independently replayed, fixed-origin degree-three representatives in the
genus-three Jacobian. Their quotient points occupy every two-fiber character
cell: ++:30, +-:26, -+:40, --:45. Every cell supplied a prime-field cubic
input pair, a successful two-cubic reduction, and a nonspecial final
divisor. Thus the ++ condition used to select known-scalar queries is not
a necessary condition for this cubic-section reduction on these classes.

The [previous line-support run](../run-20261015T2201Z/REPORT.md) froze the 141
distinct nonzero first-quotient points as exact \(r\)-torsion P-192
coordinates. The field prime and subgroup order are
\[
p=6277101735386680763835789423207666416083908700390324961279,\qquad
r=6277101735386680763835789423176059013767194773182842284081.
\]

## Fiber descent and class identity

For each quotient point \(Q\in E(\mathbf F_p)[r]\), let
\(H=[2^{-1}\bmod r]Q\), so \(2H=Q\). Above each of \(Q,H\), take both
points in the degree-two inverse cover fiber and their three \(p\)-Frobenius
images. Call the resulting six-point divisors \(F_Q,F_H\). If the relevant
quartic right-hand side is a square in \(\mathbf F_{p^3}\), the third
Frobenius iterate fixes the individual fiber points. If it is a nonsquare,
that iterate swaps the two signs. In either case, Frobenius permutes the full
six-point divisor. The [Lean theorem](../../lean/LineFiberDescent.lean)
formalizes this fixedness for both third-iterate actions on an additive
divisor group; it keeps the concrete cover action as an explicit geometric
input.

The [cubic constructor](../../line_quotient_cubic_probe.py) interpolates a
prime-field cubic through each \(F_P\) and the same three rational anchor
points \(T\). The cubic cut has degree twelve and leaves a squarefree
degree-three residual \(R_P\). All 282 input residuals for the 141 \(Q,H\)
pairs pass curve incidence, cubic incidence, prime-field descent, and
squarefreeness checks. The cubic cuts imply
\[
  [R_H-R_Q]=[F_Q-F_H].
\]

Each conjugate fiber pair is the pullback of \(Q\) or \(H\) through its
matching elliptic quotient. Hence the right-hand class is a sum of
pullbacks of \(Q-H=H\) and is killed by \(r\). The
[independent quotient-image replay](../../verify_line_quotient_images.py)
checks both points in every matching fiber and all three pushforward
components for every class; each image vector is \((Q,Q,Q)\). Since \(Q\ne0\)
and \(r\) is prime, these are nonzero classes of exact order \(r\).

## Fixed-origin reduction

The [reducer](../../line_quotient_fixed_origin.py) uses one fixed rational
degree-three origin divisor \(T\). For each class it constructs a denominator
cubic cut \(G=B+T+U\), then a numerator cubic cut \(F=A+U+D\), with
\(A=R_Q\), \(B=R_H\), and \(D\) effective of degree three. Because the cuts
belong to the same cubic linear system,
\[
  [D-T]=[B-A]=[F_Q-F_H].
\]

The [independent PARI replay](../../verify_line_fixed_origin_native.gp)
checks all 282 cubic cuts from the frozen decimal
[input](../../type_i_line_fixed_origin_native_input.gp): resultants,
factor quotients, residual coordinates, squarefreeness, overlap exclusions,
and numerator interpolation rank. The [receipt audit](../../verify_line_fixed_origin_receipt.py)
binds the Sage and PARI inputs and results and checks that every final
\(y(x)\) interpolant has a nonzero quadratic coefficient. A squarefree
degree-three \(x\)-polynomial and a nonzero quadratic \(y\)-coefficient make
the three geometric support points noncollinear. On the smooth
nonhyperelliptic plane quartic, Riemann--Roch then gives \(h^0(D)=1\), so
each \(D\) is the unique effective representative at the chosen origin.

| Character cell for \(Q,H\) | Quotient classes | Prime-field input pairs | Reduced, nonspecial classes |
| --- | ---: | ---: | ---: |
| ++ | 30 | 30 | 30 |
| +- | 26 | 26 | 26 |
| -+ | 40 | 40 | 40 |
| -- | 45 | 45 | 45 |
| Total | **141** | **141** | **141** |

## Four-cell known-scalar collector

The [collector](../../line_all_cell_collector.py) uses the same 141 reduced
quotient classes and a single deterministic stream of nonzero scalars
\(a\). Every queried target remains \(T=-2[a]G\), so its scalar equation is
known. For this target the inverse fiber pair is \(H=[a]G\) and
\(Q=[2a]G\). The control selects the first 64 queries in the ++ cell;
the broadened policy selects the first 16 in each of ++, +-, -+,
and --. Both laws use the same seed and trial stream.

The first 16 base classes have 816 distinct three-sum targets. The exact
character partition is 206, 218, 199, and 193 in the four cells. Thus
the four-cell policy admits all 816 support targets, compared with 206
under the ++ filter. In the frozen stream, it filled its balanced
64-query arm after 76 scalar trials; the 64-query ++ control required
241 trials. The 141-class exact pair-index lookup found zero matches and
rank zero in each arm. These counts measure eligibility and finite
query yield at this base size; the extra admitted support has not yet
produced relations. The [independent prime-field replay](../../verify_line_all_cell_collector.py)
recomputes the scalar targets, pair index, and both arms' hits and rank.

| Policy | Exact first-16 eligible targets | Scalar trials to fill 64 queries | Full-base matches | Rank |
| --- | ---: | ---: | ---: | ---: |
| ++ | 206 | 241 | 0 | 0 |
| Balanced four-cell | 816 | 76 | 0 | 0 |

The frozen [run manifest](RUN.json) hashes the sources, all 141 representative
records, GP and Sage replay logs, four-cell collector records, Lean checks,
and toolchain receipts. The next controlled experiment should enlarge the
usable base and measure verified relation yield per ordinary scalar trial,
then rank per verified row under both query laws.

## Lean proof map

| Result | Lean declaration | Concrete bridge |
| --- | --- | --- |
| Six-point fiber divisor fixed when third Frobenius swaps signs | P192ExtensionCover.six_fiber_fixed_of_three_swap | Sage checks all concrete fiber points and quotient maps. |
| Six-point fiber divisor fixed when third Frobenius fixes signs | P192ExtensionCover.six_fiber_fixed_of_three_fixed | Same concrete checks; the abstract result has no character assumption. |
| Four character cells total 141 | P192ExtensionCover.p192_line_quotient_character_cells | Sage computes the 141 elliptic-character cells. |
| Two cubic cuts give the fixed-origin class | P192ExtensionCover.fixed_origin_cubic_class_identity in CrossRatio.lean | Native PARI verifies all 282 concrete section cuts. |
| Nonzero quadratic residual interpolant excludes collinearity | P192ExtensionCover.quadratic_triple_noncollinear in NonspecialSupport.lean | Native PARI checks support and Python checks all 141 coefficients. |

The Lean modules are kernel checked with the pinned Lean 4.34.1 and Mathlib
revision. The finite-field curve, the concrete quotient maps, cubic-section
divisor equations, and Riemann--Roch application are written and independently
replayed here; they have not all been defined as Lean objects.
