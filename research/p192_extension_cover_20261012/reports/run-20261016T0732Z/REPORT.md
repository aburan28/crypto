# Thirteen-pivot projected subgroup-column audit at P-192

**Result.** The complete old-anchored pivot families
0–4, 50–53, and 95–98 contain 62,505 rational split non-source cuts.
Their 125,010 cofactor-projected subgroup columns are pairwise distinct,
also after sign folding, and have no identity or selected-base overlap.
Each of the 13 producer files passed the independent four-offset quartic
interpolation and prime-field line-equation replay. The resulting line
rows all have private residual columns: raw row rank is 62,505 and the
exact and sign-folded leaf cores each contain zero rows. This is an exact
statement for the 13 listed pivots; 128 of the 141 old pivots remain outside
this checkpoint.

| Pivot range | Split cuts | Projected columns | Direct partner controls |
| --- | ---: | ---: | ---: |
| 0–4 | 23,993 | 47,986 | 25 |
| 50–53 | 19,219 | 38,438 | 20 |
| 95–98 | 19,293 | 38,586 | 20 |
| **Union** | **62,505** | **125,010** | **65** |

The [sparse aggregate](../../type_i_line_trace_full_partial_sparse_20261016T0732Z.json)
lists all covered and missing pivot IDs, records hashes of every producer
and independent replay file, and counts exact and sign-folded column
collisions, selected-base hits, and leaf-core rows. It has zero unreplayed
pivots. The earlier [eight-pivot checkpoint](../run-20261016T0726Z/REPORT.md)
remains frozen. A prior sparse aggregate captured two newly completed
producer pivots before their replays; that [receipt](../../type_i_line_trace_full_partial_sparse_20261016T0731Z.json)
is preserved with `unreplayed_pivots` explicit.

The [producer](../../line_trace_projected_residual_chunk.py) constructs
quartic residuals and trace-projects the first one; the second follows
from the line-section sum and is directly checked on the first and every
1,024th split cut per pivot. The [independent verifier](../../verify_line_trace_projected_chunk.py)
recomputes the residual quadratic by interpolation, checks every selected
pair and all residual coordinates, and replays every prime-field quotient
sum. The [aggregator](../../aggregate_line_trace_projected_chunks.py)
counts subgroup-column reuse across the selected union and peels its
private columns. The [Lean private-pair result](../../lean/LineSupport.lean)
formalizes why a zero leaf core leaves the older-column assignments free.

The complete-family gate is unchanged: all 141 producer files and replays
must exist, their pair/source/split counts must agree with the
[full geometric screen](../../type_i_line_nonsource_reuse_old141.json), and
their ordered residual-coordinate stream must reproduce its SHA-256
digest. The three checked-Sage producer processes continue in disjoint
starting ranges. The [run manifest](RUN.json) freezes only completed files
for this checkpoint.
