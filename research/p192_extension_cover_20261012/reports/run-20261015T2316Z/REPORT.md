# P-192 signed short-support relation screen

**Result.** The 141 frozen, reduced Type I quotient classes have no
nontrivial equality between positive or negative supports of lengths one,
two, or three. All 12 zero-target, cross-length, and opposite-sign
intersection cells are empty. In particular, the screen produces zero
homogeneous short-relation rows and rank zero modulo the P-192 subgroup
order. This extends the [exact three-sum census](../run-20261015T2255Z/REPORT.md),
which established that the 10,011 two-sums and 477,191 three-sums are
individually distinct.

Let \(S_n=\{Q_{i_1}+\cdots+Q_{i_n}:i_1\le\cdots\le i_n\}\), with repetition
allowed, for \(n=1,2,3\). The exact computation gives
\[
(|S_1|,|S_2|,|S_3|)=(141,10011,477191),\qquad 0\notin S_1\cup S_2\cup S_3.
\]
For \(n<m\), \(S_n\cap S_m=\varnothing\). For every
\(1\le n\le m\le3\), \(S_n\cap(-S_m)=\varnothing\). Combined with
injectivity within each \(S_n\), this excludes a relation obtained by
equating two different signed unordered sums having at most three base
terms on each side. The result is specific to this 141-class base; it
does not screen four-term supports or ordinary known-scalar queries.

| Screen cell family | Cells | Nontrivial overlaps |
| --- | ---: | ---: |
| Sum equals identity, lengths 1, 2, 3 | 3 | 0 |
| Same-sign cross-length, 1:2, 1:3, 2:3 | 3 | 0 |
| Opposite-sign, \(1\le n\le m\le3\) | 6 | 0 |
| Total | **12** | **0** |

The [integer-field producer](../../line_signed_short_relation_screen.py)
enumerates the complete supports, checks all source points satisfy the
P-192 curve equation and are killed by \(r\), and binds \(S_3\) to the
previous sorted-coordinate support digest. The
[independent Sage replay](../../verify_line_signed_short_sage.py) uses
Sage's elliptic-curve group law and reconstructs all 12 intersections
from the frozen coordinates. Both preserve the 141-point source receipt,
the exact field representation, and the runtime/toolchain records.

The [Lean short-sum module](../../lean/ShortSums.lean) proves that a
same-sign list-sum equality yields a zero difference row and an
opposite-sign equality yields a zero signed-sum row for arbitrary
additive commutative groups. The concrete finite-field nonintersection
certificate is independently enumerated here, not represented as
141 Lean elliptic-curve points.

| Result | Lean declaration | Concrete bridge |
| --- | --- | --- |
| Equal list sums give a zero difference row | P192ExtensionCover.same_sign_lists_relation | Complete same-sign support joins use the frozen P-192 point list. |
| Opposite list sums give a zero signed row | P192ExtensionCover.opposite_sign_lists_relation | Complete negated-support joins use independent integer and Sage group laws. |
| No signed short overlap on the 141 classes | Open concrete finite-field bridge | Exact 12-cell enumeration and support digest, replayed independently. |

The run manifest hashes the source, receipts, both group-law replays,
the Lean build, axiom audit, and fresh kernel check. The first Sage
replay read a producer receipt before its digest field was written;
its failed log is preserved beside the successful rerun. Local elapsed
time is a correctness diagnostic.

The next collector test should grow the actual usable subgroup base
under a frozen construction and measure ordinary known-scalar
query overlap and useful rank. Any proposed geometry-conditioned
advantage must use the support count and ambient subgroup population
for each selected character cell, as required by the
[capacity criterion](../../../weak_curve_leaderboard_20261012/criteria.json).
