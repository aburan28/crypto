# Twenty-six-pivot projected column and leaf-core census at P-192

**Result.** The complete old-anchored pivot families 0–9, 50–57, and
95–102 contain 125,326 rational split non-source cuts and 250,652
cofactor-projected subgroup columns. Every column is distinct, including
after sign folding, and none is the identity or a selected-base column.
All 26 producer receipts have independent quartic-interpolation and
prime-field quotient-line replays. Their line equations have 125,326
private residual-column pivots, giving raw row rank 125,326 and zero rows
in the exact and sign-folded leaf cores. This conclusion covers exactly
the listed 26 of 141 old pivots.

| Pivot range | Split cuts | Projected columns | Direct partner controls |
| --- | ---: | ---: | ---: |
| 0–9 | 47,922 | 95,844 | 50 |
| 50–57 | 38,690 | 77,380 | 40 |
| 95–102 | 38,714 | 77,428 | 40 |
| **Union** | **125,326** | **250,652** | **130** |

The [sparse aggregate](../../type_i_line_trace_full_partial_sparse_20261016T0750Z.json)
has hashes for all 26 completed producer and replay files, zero
`unreplayed_pivots`, exact and sign-folded collision counts, selected-base
overlap checks, and leaf-peeling counts. The [independent pivot verifier](../../verify_line_trace_projected_chunk.py)
recomputes the quartic residuals from four polynomial evaluations and
checks every selected pair and quotient sum. The [global auditor](../../aggregate_line_trace_projected_chunks.py)
reconstructs the projected-column incidence rows and peels degree-one
columns. The [Lean private-pair theorem](../../lean/LineSupport.lean)
states the relevant row-rank and older-column freedom implication.

The 141-pivot full census remains the next gate. Its completed receipt
must match the [full geometric screen](../../type_i_line_nonsource_reuse_old141.json)
in all pair, source, split, and nonsplit counts and in the ordered
residual-coordinate SHA-256 stream. The three checked-Sage producer
processes are continuing; their source and parameters are frozen in
the [prior run](../run-20261016T0720Z/RUN.json). The [new run manifest](RUN.json)
hashes only the 26 completed pivots, their replays, and this audit.
