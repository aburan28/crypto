# P-192 cross-pivot subgroup-column screen

**Result.** The first two residual-pivot fans over the descended P-192
genus-three quartic have no exact or sign-folded projected-column overlap
with each other. Pivot 001 has 4,971 rational split cuts and 9,941 distinct
projected subgroup columns from 9,942 residual presentations. Its only
overlap with the complete 1,360,948-column old-anchored family is the
partner on the source line `(1, 143)`. The corresponding result for pivot
000 is the source line `(0, 285)`. Every other row of both fans retains
private residual columns, so these two fans add no constraint on the
9,925 selected-class logarithms after their source-line restatements are
removed.

| Quantity | Pivot 000 | Pivot 001 | Cross-fan |
| --- | ---: | ---: | ---: |
| Rational split cuts | 5,037 | 4,971 | 10,008 total |
| Residual presentations | 10,074 | 9,942 | 20,016 total |
| Distinct signed projected columns | 10,073 | 9,941 | 0 overlap |
| Old-family exact / signed overlaps | 1 / 1 | 1 / 1 | Both source-line partners |
| Source-line restatements | 2 | 2 | 4 total |
| Non-source rows with private residuals | 5,035 | 4,969 | 10,004 total |

The [first fan](../../type_i_residual_pivot000_three.json) was frozen and
independently replayed in the earlier run. The new
[parameterized producer](../../line_residual_pivot_cross_fan.py) takes the
first residual from old-anchored pivot 001's first source cut. For each of
the 9,925 selected quartic points it reconstructs the line section,
projects one residual by three-step Frobenius trace, and obtains its
partner from the constant four-point line sum. Every 256th split row
projects both residuals directly through the six-step method; pivot 001
has 20 such controls. Its
[pivot-001 receipt](../../type_i_residual_pivot001_cross.json) contains all
4,971 rows and records row SHA-256
`bd3323100978ccec9bef128007c31f34bc8777f6c1d6451b4a72ea40007d7ea7`.
That digest matches a separate in-memory execution of the same fan.
The independent [prime-field replay](../../verify_line_residual_pivot001_cross.py)
passed and retained its [receipt](../../type_i_residual_pivot001_cross_replay.json).

The selected elliptic quotient inputs receive a direct modular check of
their affine P-192 curve equation before Sage point construction. All
9,925 are on the curve and distinct. Construction with
`E0.point(..., check=False)` then avoids repeating the expensive
equation check. In an unisolated 64-point setup probe the checked
constructor took 13.043 seconds and the preflighted unchecked
constructor took 0.270 milliseconds; the resulting point objects agreed.
These are setup-stage diagnostics, not an isolated full-producer
wall-time comparison. The accepted Sage runtime for the versioned run
is frozen in the [runtime receipt](../../type_i_residual_pivot001_cross_runtime_info.json).
The producer guards the two import-time receipt writes and proceeds only
when their proposed bytes match the frozen files.

An independent prime-field geometry screen found 4,971 pivot-001 split
cuts and 9,941 distinct residual quartic points. Among all 1,360,948
old-family residual presentations, only the source-line partner
intersects that set. The independent replay reconstructs every line
polynomial from four separate quartic evaluations, checks every
residual point on the quartic, checks every reported subgroup point on
P-192, replays the constant line sum, and recomputes both old-family and
first-fan exact and sign-folded key intersections. The Sage producer's
20 direct controls check the sampled projection bridge.

For the rank implication, delete the two presentations of each pivot's
source line, since they restate that line's old-family equation. Each
remaining row has two columns occurring in no other retained fan row,
old-family row, or selected-base position. Given arbitrary selected-class
logarithms, choose one residual log per row and solve the other from
the line equation. The general implication is the kernel-checked Lean
theorem `P192ExtensionCover.line_private_pairs_extend` in
[LineSupport.lean](../../lean/LineSupport.lean), whose recorded axiom
audit uses only standard Lean foundations. The concrete column
distinctness is supplied by the finite-field receipts and the
[two-pivot private-core verifier](../../verify_two_pivot_private_core.py).
Its [certificate](../../type_i_two_pivot_private_core.json) checks all
10,004 non-source rows, both signed columns per row, and all 20,014
distinct signed fan columns. The entire quartic incidence graph is not
encoded as one Lean object.

This is a necessary screen for danger from this line collector:
projected-column reuse must survive leaf peeling and close a relation
cycle before it can constrain the selected logarithms. A single
residual-column collision in a row with another private column is
insufficient. Further pivots should be admitted for expensive
projection only when their quartic geometry or quotient symmetry
offers a plausible cross-family reuse mechanism; their exact signed
column core must still be measured.

The [run manifest](RUN.json) freezes the producer, independent verifier,
runtime, inputs, receipts, report, and Lean proof dependency. Replay from
the repository root:

```sh
/Volumes/SSD990/cryptanalysis/sage --runtime-info
/Volumes/SSD990/cryptanalysis/sage -python research/p192_extension_cover_20261012/line_residual_pivot_cross_fan.py --pivot 1 --compare research/p192_extension_cover_20261012/type_i_residual_pivot000_three.json --out research/p192_extension_cover_20261012/type_i_residual_pivot001_cross.json
python3 -B research/p192_extension_cover_20261012/verify_line_residual_pivot001_cross.py
python3 -B research/p192_extension_cover_20261012/verify_two_pivot_private_core.py
```
