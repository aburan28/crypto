# Cross-range projected residual-column audit at P-192

**Result.** Eight complete old-anchored pivot families at P-192—pivots
0–3, 50–51, and 95–96—contain **38,383 rational split non-source cuts**.
Their **76,766 cofactor-projected subgroup columns are distinct**, remain
distinct after sign folding, and have no identity or signed overlap with
the 9,925 selected classes. All eight pivot files passed independent
four-offset quartic interpolation and prime-field line-equation replay.
Across their union, every new line row has a private projected residual
column, so the raw new-row rank is 38,383 and the projected leaf core has
zero rows. This is an exact result for the declared eight-pivot selection;
the complete family contains 141 old pivots and 680,474 split cuts.

| Frozen pivot subset | Split cuts | Projected columns | Direct partner controls | Signed repeats |
| --- | ---: | ---: | ---: | ---: |
| 0–3 | 19,156 | 38,312 | 20 | 0 |
| 50–51 | 9,643 | 19,286 | 10 | 0 |
| 95–96 | 9,584 | 19,168 | 10 | 0 |
| **Union** | **38,383** | **76,766** | **40** | **0** |

The [sparse aggregate](../../type_i_line_trace_full_partial_sparse_20261016T0725Z.json)
records the exact pivot IDs, a hash for every complete producer chunk and
independent replay, separate exact and signed collision counts, selected
overlaps, and leaf-peeling counts. Its status is
`partial_sparse_p192_projected_leaf_core_census`; the absent 133 pivots
are explicit. The four-pivot contiguous-prefix
[aggregate](../../type_i_line_trace_full_partial_4_retry.json) separately
records 19,156 split cuts and a zero core. The first attempt to aggregate
that prefix encountered already completed distant pivots and stopped;
its failed log is retained. The corrected partial mode stops at the first
missing pivot, while `--include-sparse` deliberately audits all available
pivots without promoting the result to a complete census.

Each [producer](../../line_trace_projected_residual_chunk.py) receipt is
atomic and immutable once complete. The first residual is projected by the
proved six-step Frobenius-trace method; its partner follows from the common
line-section sum and is directly checked on the first and every 1,024th
split cut of its pivot. The [pivot verifier](../../verify_line_trace_projected_chunk.py)
independently recovers the residual quadratic from four evaluations of the
quartic, enumerates all selected pairs for that pivot, checks all reported
quartic coordinates, and replays every prime-field quotient sum. The
[global auditor](../../aggregate_line_trace_projected_chunks.py) hashes
the actual projected columns, folds signs, tests selected-base overlap,
and peels private residual columns. The [Lean private-pair theorem](../../lean/LineSupport.lean)
explains why the zero core leaves older-column assignments free for this
selection. The maximal three-sum capacity bound remains the separate
conditional theorem in [FullLineCapacity.lean](../../lean/FullLineCapacity.lean).

The full-census gate remains precise: every pivot file and independent
replay must exist; the sum of pair, source, split, and nonsplit counts must
match the [complete geometric screen](../../type_i_line_nonsource_reuse_old141.json);
and concatenating its residual coordinates in pivot/neighbor order must
reproduce that screen's SHA-256 stream. Only then will the aggregator issue
`verified_full_p192_projected_leaf_core`. Three checked-Sage producer
processes are continuing across the remaining ranges, with completed
pivots durable in the chunk directory. The [run manifest](RUN.json)
hashes the eight-pivot checkpoint independently of those live files.
