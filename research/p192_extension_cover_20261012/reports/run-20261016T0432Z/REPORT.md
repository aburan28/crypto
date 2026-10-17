# P-192 Type I: new reduced classes and full character-guided support

**Result.** The first 64 of the 9,784 newly projected P-192 quotient points
have explicit fixed-origin degree-three representatives on the descended
smooth genus-three plane quartic. The two cubic sections for every class
passed an independent PARI replay, and all 64 final divisors are nonspecial.
Together with the previously replayed 141 line-derived classes, this gives
205 distinct reduced quotient classes under the frozen projection policy.
The 64 new classes occupy every two-fiber character cell: `++` 13, `+-` 19,
`-+` 11, and `--` 21.

The first 64 points also expose a **31-rank three-sum collision family**.
They are paired residuals from 32 line sections through one old quartic
point. Cancelling that shared point makes the other old point and the two
new points on each line a decomposition of one common quotient target.
All 32 decompositions and rank 31 are checked in the
[pencil certificate](../../type_i_line_pencil_three_sum.json). Each line's
new-point pair is private, so these rows do not yet determine an old-base
logarithm.

The enlarged 205-class base has **1,456,904 distinct three-sum targets** from
1,456,935 unordered triples. Its only 31 repeated decompositions are exactly
the common-point pencil above, and their incidence rank is 31. The complete
character census places 364,036 targets in `++` and 1,092,868 in the other
three cells. Thus allowing all four cells increases **admissible support**
from 364,036 to 1,456,904 targets while every query still has the known
equation \(T_a=-2[a]G\). On the same frozen 64-query arms, both the `++`
control and the balanced four-cell law found zero exact matches and rank
zero. The balanced arm required 76 scalar trials; the control required 241.
The earlier 141-class base has 477,191 targets, partitioned as shown below.
These are unisolated coverage and collector diagnostics for the declared
bases.

| Character cell, ordered as \((Q,H)\) | 64 new reduced classes | 141-base targets | 205-base targets |
| --- | ---: | ---: | ---: |
| `++` | 13 | 119,253 | 364,036 |
| `+-` | 19 | 119,179 | 364,258 |
| `-+` | 11 | 119,403 | 363,786 |
| `--` | 21 | 119,356 | 364,824 |
| Total | **64** | **477,191** | **1,456,904** |

The field is the standard P-192 prime field and the quotient generator has
prime order
\[
p=6277101735386680763835789423207666416083908700390324961279,
\qquad
r=6277101735386680763835789423176059013767194773182842284081.
\]
The exact quartic coefficients, quotient-point encoding, deterministic
scalar stream, and checked Sage runtime are retained in the linked inputs.

## Geometric reduction

The [frozen quotient point list](../../type_i_line_full_new_quotient_points.json)
contains 9,784 exact \(r\)-torsion coordinates, all distinct from the older
141 and from their negatives under the earlier projection replay. This run
selects entries 0 through 63, in that order. For each quotient point \(Q\),
it computes the unique \(r\)-torsion half \(H=[2^{-1}\bmod r]Q\). The
archived inverse-cover constructor takes both signs in each fiber and their
three Frobenius conjugates, then cuts the quartic with a cubic through those
six points and three fixed rational anchors. The resulting prime-field
degree-three residuals are \(A=R_Q\) and \(B=R_H\).

The [new wrapper](../../line_new_quotient_reduction.py) runs the archived
[fiber constructor](../../line_quotient_cubic_probe.py) and
[fixed-origin reducer](../../line_quotient_fixed_origin.py) on this new
slice, retaining their source hashes. With the same effective degree-three
origin \(T\), the reducer finds cubic sections cutting
\(G=B+T+U\) and \(F=A+U+D\). Both cuts lie in
\(|\mathcal O_C(3)|\), so \([D-T]=[B-A]\), the intended quotient class.
The [receipt](../../type_i_line_new64_fixed_origin.json) records every
input cubic, both section coefficients, the intermediate degree-six divisor,
and the final triangular \((h_D(x),y_D(x))\) representative.

The independent [GP input](../../type_i_line_new64_fixed_origin_input.gp)
is generated from that receipt by a separate serializer. Its
[native verifier](../../verify_line_new_fixed_origin_native.gp) recomputes
all 128 degree-twelve resultants, exact residual quotients, residual
\(y\)-coordinates, squarefreeness and overlap conditions, and the
rank-nine numerator interpolation. The [native log](../../type_i_line_new64_fixed_origin_native_verification.log)
ends `NEW_LINE_FIXED_ORIGIN_NATIVE_PASS cases=64 prime_bits=192`.
The [receipt auditor](../../verify_line_new64_receipt.py) checks input and
source hashes, the checked Sage runtime, each quotient coordinate and
order, the cell labels, and nonzero quadratic \(y_D\)-coefficients for all
64 representatives. For a squarefree degree-three support on a smooth
plane quartic, this noncollinearity condition gives \(h^0(D)=1\).

## Query law and exact support

The [full support program](../../line_full_character_support.py) enumerates
all \(\binom{143}{3}=477191\) unordered three-sums of the 141 old quotient
classes. It checks that all resulting targets are distinct. For each target
\(T\), it computes \(H=-[2^{-1}\bmod r]T\), \(Q=2H\), and the two
prime-field character signs of the frozen sextic at \(x(Q),x(H)\). PARI's
Kronecker symbol supplies the complete census, with 948 spaced Python Euler
checks and all 141 archived one-point eligibility cells as controls. Its
[receipt](../../type_i_line_full_character_support.json) reproduces the
earlier first-16 cell counts `206, 218, 199, 193` and the full counts in the
table above.

The same program independently reconstructs the archived 64-query arms.
For each recorded scalar \(a\), it verifies \(T_a=-2[a]G\), the corresponding
\((Q,H)\) cell, lookup in the full three-sum support, and rank of the
resulting incidence rows. On the 141-class base, the `++` arm admits only
119,253 targets; the balanced arm draws 16 known scalars from each cell,
making all 477,191 targets eligible under its law. Both specific 64-query
samples have zero hits.

The [205-class collector](../../line_new64_character_collector.py) repeats
the full enumeration after appending the 64 reduced quotient classes. It
checks every one of the 1,456,935 unordered triples, records the first
decomposition and every repeated decomposition, and computes collision-row
rank modulo \(r\). The [receipt](../../type_i_line_new64_character_collector.json)
has exactly 31 duplicates and rank 31; since the independently checked
pencil already provides 31 duplicate decompositions, it accounts for all
collisions in this base. Its four-cell census uses PARI characters with 727
spaced Python Euler controls and checks every new class's archived cell.
The same 64 scalar queries are tested against the enlarged support with
\(T_a=-2[a]G\). The `++` law admits 364,036 targets; the four-cell law
admits 1,456,904. Both frozen arms still have zero matches and rank zero.
The changed query law increases eligible support, while a larger base or a
different relation mechanism is needed to test useful rank and complete
single-target cost.
At fixed *accepted-query* count, the four nearly equal cell supports give
nearly the same ideal-uniform hit probability per query as the `++` control.
The gain is per raw scalar trial: the balanced arm admits four cells and
filled its 64-query quota in 76 trials, versus 241 for the `++` control.
Neither count is a wall-time comparison.

The existing genus-seven and genus-thirty-one character-cover estimates
also bound **each** sign cell. Doubling permutes the odd prime-order
elliptic group, so the character sum for the second fiber equals the
first. For signs \(\epsilon,\eta\in\{\pm1\}\), the full-group cell count
obeys
\[
4(A_{\epsilon,\eta}+\delta_{\epsilon,\eta})
 = r+(\epsilon+\eta)S_1+\epsilon\eta S_{12},
\qquad 0\le\delta_{\epsilon,\eta}\le1.
\]
The prior certified bounds \(|S_1|\le14U+|t|\) and
\(|S_{12}|\le62U+|t|\), with \(U=\lceil\sqrt p\rceil\), therefore give the
same nonzero-scalar lower bound
\(L=1569275433846670190958947353987675545133315652035359680560\)
for every cell. The [four-cell arithmetic certificate](../../type_i_line_four_cell_bound.json)
checks \(364824\cdot2^{171}<L\). Under ideal uniform sampling conditional
on a chosen cell, the 205-class support hit probability is below
\(2^{-171}\) per query; the union bound for 64 such queries is below
\(2^{-165}\). The [Lean theorem](../../lean/FourCellDensity.lean) proves the
four-sign implication with the count identity, character-sum bounds, and
identity exception explicit. The concrete cover genera and Hasse--Weil
application remain in the earlier written geometric proof and independently
checked finite-field inputs.

## Common-point pencil and useful-rank boundary

The [pencil checker](../../line_pencil_three_sum.py) uses the frozen old and
new quartic coordinates, their source-pair indices, and both sets of
projected quotient coordinates. For each of the first 32 split source
lines, it verifies that the two new quartic points are collinear with old
point 1 and the line's other old point. In the quotient group the 32 triples
\(Q_{j_k}+N_{2k}+N_{2k+1}\) all equal one point. Comparing any 31 against
the first gives 31 independent incidence rows over the prime subgroup order:
row \(k\) has its own new-point pair. A separate
[PARI replay](../../verify_line_pencil.gp) reads independently serialized
[prime-field points](../../type_i_line_pencil_gp_input.gp), checks all 32
elliptic additions and the matrix rank, and ends with
`PENCIL_GP_REPLAY_PASS lines=32 triples=32 rank=31 prime_bits=192` with empty
stderr. Its two preliminary parser failures are retained as failed logs.
This supplies verified relation
collisions, but the same private columns allow every assignment of the 141
old logarithms to extend to the 64 new columns. The 31-row system has
nullity \(205-31=174\): after fixing all old columns, 33 new-column
degrees of freedom remain. The next collector must connect
these line pencils through reused residual classes or other divisor
relations before their rank constrains old logs.

## Proof and replay

| Mathematical or computational step | Check |
| --- | --- |
| \([D-T]=[B-A]\) from the two cubic section cuts | `P192ExtensionCover.fixed_origin_cubic_class_identity` in [CrossRatio.lean](../../lean/CrossRatio.lean); the current build, axiom audit, and fresh kernel replay passed. |
| Concrete section divisors for entries 0–63 | Checked Sage producer and independent GP replay of 128 cuts. |
| Nonspecial degree-three support | GP checks squarefree divisors; receipt audit finds 64 nonzero quadratic \(y_D\)-coefficients. |
| Three-sum character partition | Complete 477,191-target enumeration; 948 Euler controls and 141 archived cell controls. |
| Enlarged-base support and exact duplicate set | Complete 1,456,935-triple enumeration, 727 Euler controls, 64 new-class cell controls, and 31 duplicate decompositions accounted for by the pencil. |
| All four character cells have the same certified density floor | `four_cell_character_lower_bound` and `p192_four_cell_nonzero_scalar_lower_bound` in [FourCellDensity.lean](../../lean/FourCellDensity.lean); build, axiom audit, and fresh kernel replay passed. The exact P-192 inputs and integer hit inequality are checked in the arithmetic certificate. |
| Known-scalar collector equation | Independent prime-field point arithmetic and direct membership/rank replay for both frozen 64-query laws. |
| Two lines through one selected point give equal residual three-sums | `P192ExtensionCover.pencil_three_sum_eq_of_line_sections` in [PencilThreeSum.lean](../../lean/PencilThreeSum.lean); build, axiom audit, and fresh kernel replay passed. The 32 concrete line and quotient identities are replayed in Python and PARI. |
| Private new columns leave old assignments free | `P192ExtensionCover.line_private_pairs_extend` in [LineSupport.lean](../../lean/LineSupport.lean); the 31 concrete rows have distinct private new-point pairs. |

Reproduce the geometric run through the repository's checked Sage launcher:

```sh
/Volumes/SSD990/cryptanalysis/sage --runtime-info
/Volumes/SSD990/cryptanalysis/sage -python research/p192_extension_cover_20261012/line_new_quotient_reduction.py --start 0 --count 64 --out research/p192_extension_cover_20261012/type_i_line_new64_fixed_origin.json
python3 research/p192_extension_cover_20261012/prepare_line_new_fixed_origin_gp.py
gp -q research/p192_extension_cover_20261012/verify_line_new_fixed_origin_native.gp
python3 research/p192_extension_cover_20261012/verify_line_new64_receipt.py
python3 research/p192_extension_cover_20261012/line_full_character_support.py --out research/p192_extension_cover_20261012/type_i_line_full_character_support.json
python3 research/p192_extension_cover_20261012/line_new64_character_collector.py
python3 research/p192_extension_cover_20261012/line_four_cell_bound.py
python3 research/p192_extension_cover_20261012/line_pencil_three_sum.py
python3 research/p192_extension_cover_20261012/prepare_line_pencil_gp.py
gp -q research/p192_extension_cover_20261012/verify_line_pencil.gp
```

The Sage [runtime receipt](../../type_i_line_new64_runtime_info.json),
[producer log](../../type_i_line_new64_fixed_origin.log),
[first failed wrapper log](../../type_i_line_new64_fixed_origin.failed.log),
[GP log](../../type_i_line_new64_fixed_origin_native_verification.log), and
[character-census log](../../type_i_line_full_character_support.log) preserve
the command outcomes. The failed wrapper run stopped after computing the
64 cubic pairs because it omitted the archived reducer's `cell` field; it
produced no reduction receipt. The corrected source reran all 64 inputs.
These timings are local correctness diagnostics. The exact finite-field
geometry and group calculations are independently replayed, while Lean
formalizes the conditional section-divisor identity; the concrete quartic
and quotient maps are not yet fully represented as Lean objects.
