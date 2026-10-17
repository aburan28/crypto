# Exact three-sum coverage of the 141-class P-192 quotient base

**Result.** The 141 fixed-origin P-192 Type I quotient classes have
477,191 distinct nonidentity unordered three-sums, exactly the maximum
\(\binom{141+2}{3}\). Their 10,011 unordered pair sums are also all
distinct. Independent Sage and integer-field group-law enumerations agree
on the complete three-sum support digest
\(\texttt{95a111382b5b3088a11f438e564864a462e46fb83619ed9b170afda03e067fb0}\).
Consequently no two- or three-sum collision yields a homogeneous relation
on this frozen base.

The [preceding all-character reduction](../run-20261015T2235Z/REPORT.md)
certified nonspecial degree-three Jacobian representatives for all 141
classes and compared two known-scalar query laws. This run measures the
complete support against which those laws query. The subgroup order is
\[
r=6277101735386680763835789423176059013767194773182842284081.
\]
Every base point is nonidentity and killed by \(r\); the independent
verifier rechecks the exact P-192 curve equation and each \(rP=O\).

## Support theorem and uniform-query criterion

For \(B\) base classes in an abelian group, unordered three-sums with
repetition are images of size-three multisets on a \(B\)-element set.
Therefore
\[
|\{Q_i+Q_j+Q_k:i\le j\le k\}|\leq\binom{B+2}{3}.
\]
The [Lean theorem](../../lean/ThreeSumCapacity.lean) formalizes this
image-cardinality bound for arbitrary additive commutative monoids, and
its \(B=141\) specialization. Equality for the frozen P-192 points is
an exact computational result, replayed with two different group-law
implementations.

For ideal uniform \(a\in\{1,\ldots,r-1\}\) and generator \(G\), the map
\(a\mapsto -2[a]G\) permutes the nonidentity \(r\)-torsion points. The
complete support excludes identity, so the exact one-query hit
probability is
\[
\frac{477191}{r-1}\simeq 2^{-173.136}.
\]
For 64 independent ideal-uniform queries, the union bound is
\(64\cdot477191/(r-1)\simeq2^{-167.136}\).
As a capacity-only necessary condition for expected hit count at least
one in 64 ideal-uniform queries, a base would need at least
8,379,989,631,760,464,847 usable classes: that is the smallest integer
\(B\) with \(64\binom{B+2}{3}\ge r-1\). It is a counting threshold,
not a construction or a sufficiency result.

The character-conditioned query laws in the preceding run draw from
different subsets of scalars. Their hit probabilities require each
cell's full subgroup population as denominator and each cell's full
support count as numerator. The measured first-16 support partition
(206, 218, 199, 193) changes eligibility but alone does not establish
a conditional hit-rate gain.

| Exact frozen base result | Count or rank |
| --- | ---: |
| Base classes | 141 |
| Unordered pairs with repetition | 10,011 |
| Distinct pair targets | 10,011 |
| Pair-collision relation rank | 0 |
| Unordered triples with repetition | 477,191 |
| Distinct nonidentity three-sum targets | 477,191 |
| Triple-collision relation rank | 0 |
| Distinct first-16 three-sum targets | 816 |

The [Sage producer](../../line_three_sum_coverage.py) enumerates all
477,191 triples and saves a sorted 24-byte-coordinate support hash. The
[independent verifier](../../verify_line_sum_coverage.py) uses integer
modular arithmetic for point addition, doubling, and scalar multiplication;
it rebuilds all 10,011 pairs and 477,191 triples, checks every base point's
curve equation and \(r\)-annihilation, and matches the support hash. The
[pair-rank producer](../../line_pair_collision_rank.py) constructs any
homogeneous collision rows over \(\mathbf F_r\). Its initial assertion
confused the number of pair-lines, 9,870, with the number of unordered
pair sums, 10,011; the failed log is retained and the corrected replay
checks the latter quantity.

The run manifest hashes the input, sources, receipts, independent replay,
Lean build and fresh-kernel logs, and toolchain records. Local elapsed times
are correctness diagnostics under an unisolated host. The next measurement
is a frozen, independently sampled full-base support-cell census and
query-cell population estimate, followed by rank-per-ordinary-query
comparison on a substantially larger usable base.

## Character-cell enrichment diagnostic

A frozen [512+512 sample](../../line_cell_enrichment_sample.py) draws
uniformly from the 477,191 distinct three-sum targets without replacement
and from nonzero scalar targets with replacement. Exact cubic-field
quadratic characters give these counts:

| Fiber cell | Support sample | Uniform-scalar sample | Ratio of sampled shares |
| --- | ---: | ---: | ---: |
| ++ | 131 | 125 | 1.048 |
| +- | 103 | 130 | 0.792 |
| -+ | 140 | 133 | 1.053 |
| -- | 138 | 124 | 1.113 |

Both samples contain 512 targets and no exceptional fibers. The ratio
estimates the multiplier \(p_{\mathrm{support},c}/p_{\mathrm{ambient},c}\)
that a cell-conditioned scalar law would apply to the uniform hit
probability. The ++ point estimate is 1.048; these samples do not
resolve a useful enrichment. A balanced four-cell law has point
estimate \((1.048+0.792+1.053+1.113)/4\approx1.002\) relative to
ideal uniform sampling. These are exploratory estimates with
binomial-scale sampling error, not exact full-cell populations or a
measured relation rate. The [receipt audit](../../verify_line_cell_sample.py)
replays the seed, triple indices, scalar stream, and count tables;
the exact character labels are produced by the checked Sage geometry
job.

## Lean proof map

| Mathematical result | Lean declaration | Verification boundary |
| --- | --- | --- |
| Three-sum support has at most \(\binom{B+2}{3}\) elements | P192ExtensionCover.three_sum_support_card_le | General additive commutative monoid, no curve assumptions. |
| \(B=141\) capacity is 477,191 | P192ExtensionCover.p192_base_141_three_sum_capacity | Exact arithmetic in Lean. |
| Frozen base reaches the capacity | Open concrete bridge | Sage and independent integer-field enumeration of all points and sums. |
| Uniform scalar map is bijective and gives exact hit probability | Written proof from \(r\) odd prime and generator \(G\) | Concrete elliptic group/generator is checked computationally; no full Lean elliptic bridge. |
