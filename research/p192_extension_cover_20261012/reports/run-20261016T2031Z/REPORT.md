# P-192 second-generation residual-pivot line fan

**Result.** Pivoting on the first projected residual of the complete
old-anchored Type I census produces 5,037 rational split cuts through the
9,925 selected quartic points. Their 10,074 projected residual presentations
give 10,073 exact and sign-folded subgroup columns. The sole repeated column,
the two selected-base hits, and the sole overlap with the 1,360,948-column
old-anchored family all occur on the original line that created the pivot.
The other 5,035 split rows each have two new private residual columns. Thus
this complete second-generation fan adds no constraint on the selected-class
logarithms after its private columns are eliminated.

| Quantity | Count |
| --- | ---: |
| Selected points paired with the residual pivot | 9,925 |
| Rational split / nonsplit cuts | 5,037 / 4,888 |
| Projected residual presentations / distinct signed columns | 10,074 / 10,073 |
| Selected-base hits / within-fan repeats | 2 / 1 |
| Old-family exact and signed column overlaps | 1 / 1 |
| Additional split rows with private residual pairs | 5,035 |
| Six-step direct projection controls in each producer | 20 |

The pivot is the first residual of old-anchored row `(0, 285)` in the frozen
[`pivot-000` receipt](../../type_i_line_trace_full_chunks/pivot-000.json).
The [fully checked producer](../../line_residual_pivot_three.py) takes a
three-step Frobenius trace for every first residual and verifies its
`p^3`-fixed coordinates. It derives the second subgroup image from the
constant line-section sum. At every 256th split cut it projects **both**
residuals directly by the archived six-step method. Its
[receipt](../../type_i_residual_pivot000_three.json) freezes all 5,037
quartic intersections, subgroup coordinates, source hashes, checked Sage
[runtime](../../type_i_residual_pivot000_three_runtime_info.json), and
old-family membership results.

The [independent Python replay](../../verify_line_residual_pivot_three.py)
reconstructs every line polynomial from four separate quartic evaluations,
checks all residual points against the quartic, checks every projected point
against the P-192 elliptic curve, replays every four-point line sum, and
streams all 141 archived old-family pivot receipts for exact and signed
column intersections. Its [verification receipt](../../type_i_residual_pivot000_three_replay.json)
reports `RESIDUAL_PIVOT_THREE_REPLAY_PASS 5037 [[0, 0, 1, [0, 0]]]`.
The independent replay checks reported subgroup points and line equations;
the 20 native direct controls check projection against the six-step inverse.

On the source line, the residual pivot `P` and its partner `P'` complete
`Q_0 + Q_285 + P + P' = S`. Lines from `P` to selected points 0 and 285
restate this same equation, with `P'` as their shared residual. Every other
split row has its own two residual subgroup columns, disjoint even under
sign reversal from the selected and old-family columns. Given any selected
logarithms, choose one value per private pair and solve the other from its
line equation. The old source row determines `P'` after choosing `P`, and
its two restatements add no equation. The existing Lean theorem
[`LineSupport.line_private_pairs_extend`](../../lean/LineSupport.lean)
formalizes the private-pair extension implication. The concrete disjointness
and source-line identities are the point-level receipts above; the entire
finite-field quartic and incidence graph are not encoded as one Lean object.

The [sampled-closure producer](../../line_residual_pivot_three_fast.py)
checks `p^3`-closure at the same 20 rows as its direct six-step controls.
It compares every projected row, count, and old-family overlap with the
fully checked [reference](../../type_i_residual_pivot000_three.json) before
writing its [receipt](../../type_i_residual_pivot000_three_fast.json).
Both runs produced exactly the same 5,037 rows. The fully checked run took
541.189 s total after model initialization, including 508.694 s in its
projection interval; the sampled-check run took 912.475 s and 854.442 s,
respectively. These unisolated runs had different host contention. They
establish exact output parity but do not establish a wall-time improvement.
The abstract `Tr_6(P)=2 Tr_3(P)` and cofactor-scalar implication is proved
in [`TraceThree.lean`](../../lean/TraceThree.lean) under its stated field
closure hypotheses; the full-check producer verifies closure for every
first residual in this pivot fan.

This pilot rules out this entire residual-pivot fan as a source of new
selected-log rank. The next geometry search should choose pivots from
different source lines and test **cross-fan projected-column reuse** before
expanding relation collection. Its preliminary gate is a signed-column
overlap or a closed cycle that survives elimination, followed by ordinary
known-scalar query yield and target-descent accounting.

Replay from the repository root:

```sh
/Volumes/SSD990/cryptanalysis/sage --runtime-info
/Volumes/SSD990/cryptanalysis/sage -python research/p192_extension_cover_20261012/line_residual_pivot_three.py
python3 -B research/p192_extension_cover_20261012/verify_line_residual_pivot_three.py
/Volumes/SSD990/cryptanalysis/sage -python research/p192_extension_cover_20261012/line_residual_pivot_three_fast.py
```
