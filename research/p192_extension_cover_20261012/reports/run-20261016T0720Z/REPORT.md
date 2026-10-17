# Full-family P-192 projection gate and maximal three-sum capacity

**Result.** A complete geometric census bounds the largest possible
three-sum base formed from the 9,925 selected P-192 Type I classes and all
1,360,948 rational outside residual presentations by **1,370,873** classes.
Even if every presentation projects to a new usable class, its unordered
three-sum capacity is at most **429,379,561,999,436,125**. Against the
certified four-cell accepted-scalar population floor, any such support has
ideal-uniform overlap below \(2^{-131}\) per accepted query. The
[Lean theorem](../../lean/FullLineCapacity.lean) proves the bound for every
actual base size and support count below those maxima, with the query-cell
denominator hypothesis explicit. It is a capacity screen, independent of
the unfinished complete projected-column census.

The first **two complete old-pivot families** have also been projected and
independently replayed. They contain 9,603 split cuts and 19,206 distinct
subgroup columns, with no identity, selected-base hit, exact repeat, or
sign-folded repeat. Their line rows all peel on private residual columns,
leaving a zero projected core for this **two-pivot prefix**. This is an
incremental result. The remaining 139 pivot files and their independent
replays are required before assigning a full-family core.

| Frozen boundary | Value |
| --- | ---: |
| Selected quotient classes | 9,925 |
| Complete geometric outside-point presentations | 1,360,948 |
| Maximal hypothetical usable base | 1,370,873 |
| Maximal unordered three-sum capacity | 429,379,561,999,436,125 |
| Four-cell accepted-scalar population floor per cell | 1,569,275,433,846,670,190,958,947,353,987,675,545,133,315,652,035,359,680,560 |
| Conditional ideal-uniform hit bound per accepted query | \(<2^{-131}\) |
| Complete and independently replayed pivots at this checkpoint | 2 of 141 |
| Split cuts and subgroup columns in that prefix | 9,603 and 19,206 |
| Signed projected repeats and signed selected-base hits in that prefix | 0 and 0 |
| Signed leaf-core rows in that prefix | 0 |

For a base of size \(B\), let \(M\) be its distinct three-sum support and
\(D\) the accepted scalar population of one chosen character cell. The
structural bound is \(M\le\binom{B+2}{3}\). The complete
[geometric receipt](../../type_i_line_nonsource_reuse_old141.json) gives
\(B\le9925+1360948=1370873\) for this construction, regardless of
projected collisions or identity exclusions. The certified four-cell
character estimate supplies \(D\ge L\), where \(L\) is the floor in the
table. Therefore \(M2^{131}<D\), exactly as checked in
`p192_old_anchor_three_sum_hit_bound`. The inequality is conditional on
ideal-uniform sampling within an accepted cell. It does not infer a relation
yield, useful rank, or end-to-end cost from the capacity bound.

## Resumable projected-column census

The [per-pivot producer](../../line_trace_projected_residual_chunk.py)
enumerates every selected pair with a declared old pivot, skips the frozen
source quartets, and tests the quartic residual discriminant. On a split
cut it computes both rational residual coordinates. It maps the first
through the checked Frobenius-trace projector and derives its partner from
the independently replayed common line-section sum
\(Q_Y=S-Q_A-Q_B-Q_X\). It directly projects the first partner and every
1,024th later partner at each pivot. Each complete pivot is written
atomically to `pivot-NNN.json`, with source and input hashes, counts, full
coordinates, a residual-stream digest, and its own timing boundary.

The [independent pivot verifier](../../verify_line_trace_projected_chunk.py)
reconstructs every cut from four quartic evaluations, checks the complete
split/nonsplit sequence, all residual coordinates, and every prime-field
group equation. Its [replay driver](../../replay_line_trace_projected_chunks.py)
can verify newly completed pivot files without changing prior receipts.
The [aggregator](../../aggregate_line_trace_projected_chunks.py) hashes
exact and sign-folded subgroup columns across pivots. On the complete
141-pivot set it will require the same selected-pair and split counts and
the identical ordered geometric residual-stream hash as the earlier full
screen; it marks a full result verified only when every pivot also has an
independent replay. It peels rows with private residual columns and retains
any surviving core for a separate rank computation.

The [pivot-0 receipt](../../type_i_line_trace_full_chunks/pivot-000.json)
has 4,829 split cuts, matching the earlier complete geometric run's first
pivot count. The [pivot-1 receipt](../../type_i_line_trace_full_chunks/pivot-001.json)
adds 4,774. Both passed the independent verifier. The frozen
[two-pivot aggregate](../../type_i_line_trace_full_partial_2.json) has
status `partial_p192_projected_leaf_core_census`; it reports 9,603 split
cuts, 19,206 distinct exact and signed subgroup columns, and zero core
rows. Subsequent pivots are separate immutable files and are not silently
included in this checkpoint.

| Mathematical or computational result | Proof or replay |
| --- | --- |
| Maximal base and three-sum count | `p192_old_anchor_maximal_base_size` and `p192_old_anchor_maximal_three_sum_capacity` in [FullLineCapacity.lean](../../lean/FullLineCapacity.lean). |
| Uniform per-cell capacity bound with explicit \(B,M,D\) assumptions | `p192_old_anchor_three_sum_hit_bound` in [FullLineCapacity.lean](../../lean/FullLineCapacity.lean). |
| Cofactor projection by six-step Frobenius trace | [TraceProjection.lean](../../lean/TraceProjection.lean), 256 paired native projections, and direct partner controls. |
| Complete projected pivot-0 and pivot-1 cut and group equations | Atomic producer receipts and independent four-offset interpolation replays. |
| Two-pivot private-column core | Exact and signed global column hashing in the partial aggregate; private-column extension theorem in [LineSupport.lean](../../lean/LineSupport.lean). |

`lake build FullLineCapacity` and its `#print axioms` audit passed using the
pinned Lean toolchain. The fresh kernel replay result is recorded in the
[proof receipt](../../lean/full-line-capacity-verification-20261016.json).
The [run manifest](RUN.json) hashes this checkpoint's frozen artifacts.
