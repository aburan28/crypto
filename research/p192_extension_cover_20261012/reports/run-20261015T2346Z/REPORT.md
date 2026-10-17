# Split-line relations have private columns on the P-192 Type I base

**Result.** The frozen line screen yields 4,892 split line cuts and 9,784
distinct new subgroup points, each new point attached to exactly one cut.
Adding the 141 old points gives 9,925 distinct, sign-distinct usable quotient
points. Every split cut supplies a four-point equation with one common
constant, but its two new columns occur in no other split-cut row. Hence the
4,892 line rows are linearly independent while their span contains no
nonzero equation involving only the old 141 columns and the constant. This
is a structural rank criterion: a large raw relation count can coexist with
zero elimination of the old logarithms.

Let (A) be the fixed plane-quartic anchor, (K) its canonical divisor,
and φ the first elliptic quotient homomorphism on the Jacobian, followed
by the frozen cofactor projection (h). For a split line with old points
(P_i,P_j) and residual rational points (R_k,S_k), define
(q(P)=h\phi([P-A])). The plane-quartic line section is canonical, so

\[
 q(P_i)+q(P_j)+q(R_k)+q(S_k)
   = h\phi([K-4A])=:C.
\]

The right side does not depend on the line. The independent GP line screen
checks all 4,892 split incidences and point distinctness; the independent
Sage projection replay checks the 9,784 projected points are nonzero,
pairwise distinct even modulo sign, and disjoint modulo sign from the old
141. The [Sage probe](../../line_leaf_relation_probe.py) evaluates this
identity on the first 32 split cuts, and an [independent prime-field
replay](../../verify_line_leaf_relation.py) checks those 32 equations with
a separate elliptic-curve group law. The common (C) in this sample is
nonidentity. The full family equation follows from the canonical-section
identity and quotient homomorphism; the concrete full-family equality has
not been enumerated in quotient coordinates.

For coefficients (c_k\in\mathbf F_r), the fresh column belonging to
(R_k) in the sum of line rows has coefficient (c_k). Therefore vanishing
of all fresh columns forces every (c_k=0). The
[Lean private-column theorem](../../lean/LineLeafRank.lean) formalizes
this implication over any semiring. It uses no assumption about the old
columns or the common constant. The 4,892-row rank conclusion depends on
the concrete incidence and quotient-distinctness checks just named.

| Evidence or bound | Value | Basis |
| --- | ---: | --- |
| Old quotient points | 141 | Frozen line-guided collector |
| Split line cuts | 4,892 | Complete GP line screen |
| New projected points | 9,784 | Full independent Sage replay |
| Sample constant-equation checks | 32/32 | Sage and independent prime-field group laws |
| Structural rank of split-line rows, with (C) as a column | 4,892 | Private fresh columns |
| Old-only rows in their span | 0 | Private-column theorem |
| Combined usable point count | 9,925 | Disjoint projected sets |
| Unordered three-sum capacity ceiling | 162,993,977,475 | (\binom{9927}{3}) |
| Uniform 64-query union-bound numerator | 10,431,614,558,400 | (64\binom{9927}{3}) |

For a uniformly sampled nonzero subgroup target, the 64-query union bound
is (10{,}431{,}614{,}558{,}400/(r-1)\), approximately (2^{-148.754}).
This is a capacity ceiling, not a measured hit rate. A character-conditioned
query law needs its own support-in-cell and cell-population counts; the
uniform denominator does not carry over to that law.

The criterion suggested by this result is to record the **private-column
or leaf fraction and the rank after eliminating leaf variables**, not just
the number of verified geometric relations. The next construction should
make new quotient points recur across independently sampled cuts or
different cover families, then measure the remaining core rank and ordinary
known-scalar target yield. A 9,925-point base is still far below the
three-sum capacity needed for frequent uniform P-192 target coverage.

## Proof map and status

| Statement | Lean declaration | Concrete bridge |
| --- | --- | --- |
| Private fresh columns force all line coefficients to vanish | `P192ExtensionCover.private_column_coefficients_zero` | Full GP incidence and independent Sage quotient-distinctness receipts. |
| 141 old plus two points per 4,892 cuts equals 9,925 | `P192ExtensionCover.full_split_line_point_count` | GP cut count and Sage point count. |
| Every split line has common (C=h\phi([K-4A])) | Open geometric formalization | Written canonical-section proof; 32 Sage and prime-field checks. |
| Full concrete P-192 4,892-row rank | Open finite-field bridge | Private-column Lean theorem plus full input receipts; no encoded 9,784-point Lean curve model. |

The run manifest pins sources, projection and incidence inputs, sample
receipts, toolchain, Lean build and axiom audit. Local elapsed times are
correctness diagnostics.
