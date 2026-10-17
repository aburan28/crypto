# Frobenius-trace projection and cross-pivot residual columns at P-192

**Result.** A Frobenius-trace projector gives the same P-192 prime-order
quotient point as direct extension-field cofactor multiplication on all 256
frozen rational residual points. The common line-section sum then determines
the partner residual column from one direct projection. Using those two
identities, the first 4,096 split cuts at old pivot 0 and the first 16 split
cuts at each of all 141 old pivots give **6,336 distinct cuts and 12,672
distinct projected subgroup columns** after accounting for their 16-line
overlap. The columns remain distinct after sign folding and have no hit in
the existing 9,925-column base. An independent prime-field replay checks
the complete selected cut prefixes, all quartic residual coordinates,
line-section group equations, and the cross-policy union.

| Frozen selection law | Split cuts | Projected columns | Direct partner checks | Exact repeats | Signed repeats |
| --- | ---: | ---: | ---: | ---: | ---: |
| First-pivot prefix | 4,096 | 8,192 | 143 | 0 | 0 |
| First 16 at each of 141 pivots | 2,256 | 4,512 | 149 | 0 | 0 |
| Union, subtracting 16 shared cuts | **6,336** | **12,672** | -- | **0** | **0** |

The earlier [complete geometric screen](../run-20261016T0601Z/REPORT.md)
has 680,474 split old-anchored non-source cuts and 1,360,948 distinct
outside quartic points. This run projects the declared 6,336-cut union into
the actual P-192 subgroup. It measures projected identity, exact equality,
sign-folded equality, and selected-base overlap rather than equating
quartic-coordinate distinctness with subgroup-column distinctness.

## Trace criterion and exact line shortcut

Let \(E/\mathbf F_p\) be the P-192 curve, let
\(r=\#E(\mathbf F_p)\) be prime, and write
\(\#E(\mathbf F_{p^6})=rh\) with \(r\nmid h\). For every
\(P\in E(\mathbf F_{p^6})\), set
\(\operatorname{Tr}_6(P)=\sum_{i=0}^{5}\pi^i P\). The trace is in
\(E(\mathbf F_p)\), hence is annihilated by \(r\). The point \([h]P\)
has order dividing \(r\) and lies in the unique \(r\)-primary subgroup,
which is \(E(\mathbf F_p)\). Thus

\[
 [h]\operatorname{Tr}_6(P)
 =\operatorname{Tr}_6([h]P)= [6h]P,
 \qquad
 [h]P=[h\,6^{-1}\bmod r]\operatorname{Tr}_6(P).
\]

For a fixed cover anchor \(A\), the producer precomputes its trace and
projects \(P-A\) as
\([h6^{-1}](\operatorname{Tr}_6(P)-\operatorname{Tr}_6(A))\) on the
prime-field curve. The [paired native run](../../type_i_line_trace_projection_256.json)
checks the trace and direct \([h](P-A)\) coordinates for all 256 frozen
residuals. In its unisolated component clocks, inverse-cover construction
took 1.072 s, direct extension-field projection 13.249 s, and the trace
projection 8.771 s on those same points. The 1.51 ratio compares only the
two projection components; it is not a full collector speedup.

Every transverse quartic line with selected endpoints \(A,B\) and outside
residuals \(X,Y\) has the same projected line-section sum \(S\), already
replayed on all 4,892 source lines. Consequently
\(Q_Y=S-Q_A-Q_B-Q_X\). The new producers project \(X\) by trace and
compute \(Q_Y\) with prime-field addition. They also project 143 and 149
partner points directly, respectively, including the first partner at
every old pivot. All comparisons agree. The first 128 cuts of the prefix
match the frozen complete-screen coordinates exactly.

The [prefix producer](../../line_trace_projected_residual_sample.py) and
[stratified producer](../../line_trace_projected_residual_stratified.py)
retain every selected pair, both quartic residual coordinates, and both
projected subgroup coordinates. The
[independent verifier](../../verify_line_trace_projected_residual_sample.py)
reconstructs the residual quadratic by evaluating the quartic at four
fixed line offsets, an arithmetic path distinct from the producers'
directional Taylor expansion. It enumerates every candidate before the
last selected cut at each pivot, checking that the claimed first-4,096 and
first-16 laws omit no eligible split cut. It checks the four quartic
intersections and the subgroup line sum for every row. The
[union verifier](../../verify_trace_residual_union.py) checks that the
16 shared rows agree byte-for-byte and hashes both exact and signed
subgroup columns across policies. Its
[receipt](../../type_i_line_trace_residual_union.json) records zero
cross-policy repeats and zero identities. The separate prime-field
replays record zero exact and signed selected-base hits.

## Rank interpretation and proof map

Each sampled non-source line contributes one equation containing two
outside subgroup columns that occur nowhere else in the sampled union or
the selected base, even after sign folding. Therefore the 6,336 new line
rows have raw rank 6,336 when those columns are admitted, while their
elimination leaves zero constraints on the older 9,925 columns. This is
the same private-pair mechanism as the source pencil, now tested across
two old-anchored selection laws. A dangerous-curve screen must measure
the **projected leaf core**, not the number of geometric relations or
their raw row rank.

| Result | Formal or computational check |
| --- | --- |
| Six-step trace commutes with scalar multiplication, and fixed points trace to six times themselves | `sixTrace_nsmul` and `sixTrace_fixed` in [TraceProjection.lean](../../lean/TraceProjection.lean). |
| Cofactor projection equals scaled trace under the stated order and fixed-point hypotheses | `sixTrace_cofactor_projection` and the concrete scalar congruence `p192_trace_scalar_certificate` in [TraceProjection.lean](../../lean/TraceProjection.lean), with the P-192 group-order and Frobenius conditions checked in the quotient receipts and 256 paired native projections. |
| Line constant determines the partner column | `line_partner_from_constant` in [TraceProjection.lean](../../lean/TraceProjection.lean), the full source-line constant replay, and 292 direct partner controls across the two laws. |
| Private columns extend any old-column assignment | `line_private_pairs_extend` in [LineSupport.lean](../../lean/LineSupport.lean), applied to the independently replayed 6,336-line union. |

The pinned Lean module built with `lake build TraceProjection`; its
`#print axioms` audit lists only standard `propext`, and
`lake env leanchecker --fresh TraceProjection` exited successfully.
The full concrete elliptic-curve group-order facts and the fixed-origin
Jacobian-to-quotient map remain coordinate and arithmetic certificates,
not one combined Lean theorem. The local component clocks are correctness
diagnostics under the recorded host conditions. The next exact measurement
is a resumable full 680,474-cut projected-column census and leaf-core
calculation; the conditional uniform birthday estimates in the prior
report are not a substitute for it.

Replay from the repository root:

```sh
/Volumes/SSD990/cryptanalysis/sage --runtime-info
/Volumes/SSD990/cryptanalysis/sage -python research/p192_extension_cover_20261012/verify_frobenius_trace_projection.py --out research/p192_extension_cover_20261012/type_i_line_trace_projection_256.json
/Volumes/SSD990/cryptanalysis/sage -python research/p192_extension_cover_20261012/line_trace_projected_residual_sample.py --split-lines 4096 --out research/p192_extension_cover_20261012/type_i_line_trace_residual_4096.json
/Volumes/SSD990/cryptanalysis/sage -python research/p192_extension_cover_20261012/line_trace_projected_residual_stratified.py --per-pivot 16 --out research/p192_extension_cover_20261012/type_i_line_trace_residual_all141x16.json
python3 research/p192_extension_cover_20261012/verify_line_trace_projected_residual_sample.py --sample research/p192_extension_cover_20261012/type_i_line_trace_residual_4096.json --out research/p192_extension_cover_20261012/type_i_line_trace_residual_4096_replay.json
python3 research/p192_extension_cover_20261012/verify_line_trace_projected_residual_sample.py --sample research/p192_extension_cover_20261012/type_i_line_trace_residual_all141x16.json --out research/p192_extension_cover_20261012/type_i_line_trace_residual_all141x16_replay.json
python3 research/p192_extension_cover_20261012/verify_trace_residual_union.py
```

The [run manifest](RUN.json) hashes the inputs, code, proofs, receipts, and
logs. The earlier 128-cut producer and replay are retained as the initial
control and are included in that manifest.
