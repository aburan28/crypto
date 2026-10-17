# Complete P-192 old-anchored projected line-incidence census

**Result.** The complete 141-pivot old-anchored Type I line family on
standard NIST P-192 has **680,474 rational split non-source cuts**.
Its **1,360,948 cofactor-projected residual subgroup columns are all
distinct**, also after sign folding, and none is the identity or one of
the 9,925 selected classes. Every line row therefore has private
residual columns. The line-row matrix has raw row rank 680,474, while
eliminating those private columns leaves **zero constraints on the
selected-class logarithms**. The exact and sign-folded leaf cores both
contain zero rows.

| Full-family quantity | Count |
| --- | ---: |
| Selected quartic points | 9,925 |
| Old pivot points | 141 |
| Ordered pivot-to-later-point pairs | 1,389,414 |
| Source-line presentations | 24,460 |
| Nonsplit non-source cuts | 684,480 |
| Rational split non-source cuts | 680,474 |
| Geometrically distinct outside residual points | 1,360,948 |
| Distinct exact and sign-folded projected columns | 1,360,948 |
| Direct second-residual projection controls | 705 |
| Replayed complete producer pivots | 141 |
| Signed or exact repeated columns; selected hits; identity columns | 0 |
| Exact and sign-folded leaf-core rows | 0 |

The [full aggregate](../../type_i_line_trace_full_verified_20261016T1849Z.json)
hashes all 141 producer and independent replay files. It checks the
pair/source/split/nonsplit totals against the earlier
[complete geometric screen](../../type_i_line_nonsource_reuse_old141.json)
and reproduces its ordered residual-coordinate SHA-256 stream,
`8590422e190d7339ee8f6377e1730997de92e67f3904a36de59b05901894f551`.
Its status is `verified_full_p192_projected_leaf_core`, with no missing
or unreplayed pivot. The [global auditor](../../aggregate_line_trace_projected_chunks.py)
then hashes each subgroup point and its negative, checks selected-base
membership, and peels residual columns from the incidence rows.

The [independent pivot verifier](../../verify_line_trace_projected_chunk.py)
enumerates every selected pair for its pivot, interpolates the restricted
quartic from four evaluations, rederives each residual quadratic and
point, checks the quartic equation, and replays the prime-field quotient
group sum. It checks all 141 complete producer chunks. The native
[producer](../../line_trace_projected_residual_chunk.py) computes the
first residual's subgroup image through the six-step Frobenius trace and
derives its partner from the common line-section sum; it directly
projects the partner at the first and every 1,024th split cut per pivot.
The prior [trace controls](../run-20261016T0654Z/REPORT.md) compare
direct and trace projection on 256 frozen points. The independent
verifier checks the point and group equations for every reported column;
it does not recompute the extension-field inverse for every column.
An independent [SQLite full-family replay](../../type_i_line_trace_full_sqlite_replay_20261016T1852Z.json)
inserts all 1,360,948 sign-canonical column keys into a unique database
index, checks every reported point against the P-192 curve equation,
checks selected-base overlap, and independently reproduces the
geometric stream digest. It finds zero duplicate, identity, off-curve,
or selected columns; its sorted signed-column digest is
`13bf9e49198a9fa9fd3421cbcf1bb34d8e25431333a0c7276fdd147c3ddc4a65`.

The rank conclusion follows from the full distinct-column certificate:
each row has two residual columns occurring in no other row or selected
base position. For any assignment of the selected-class logarithms, one
can solve each line equation by choosing a value for one of that row's
private residual logs. The [Lean private-pair theorem](../../lean/LineSupport.lean)
formalizes this implication; its hypotheses are discharged here by the
complete projected-column census. The 4,892 source-line relations are a
separate family and are not included in this non-source matrix.

This yields an exact **danger criterion for this collector design**:
repeating more old-anchored selected-pair lines cannot create useful
selected-log rank while their projected residual columns remain private.
A successor geometry must create subgroup-column reuse or another
identity that closes the residual columns, then demonstrate nonzero
eliminated core rank and ordinary known-scalar target overlap. The
[capacity theorem](../run-20261016T0720Z/REPORT.md) separately bounds
the maximal three-sum support under an ideal-uniform four-cell query law.
The next concrete search is to pivot on newly found residual classes,
construct cross-family line cycles, and audit their signed core before
charging relation collection and target descent.

The [run manifest](RUN.json) freezes this complete source, input,
producer, replay, aggregate, and proof evidence. The new
[three-step trace theorem](../run-20261016T0747Z/REPORT.md) offers a
subfield specialization for a future implementation; this census used
its original frozen six-step projector throughout.
