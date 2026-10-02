# Pre-measurement route-guard amendment

The merged protocol requested `KIC_F5_DIRECT_PACK=1` and originally required
`direct_pack_used=true` on every call. A native development fixture (n=12,
seed 17) showed `direct_pack_used=false`, `direct_unpack_used=true` before
any registered discovery or holdout seed was run. The inherited
`f5_rows_packed_full_columns` path returns `None` when selected rows do not
occupy all ambient degree-four columns; the timed F5 source then uses its
existing correct sorted-row fallback. Sparse quadratic systems can reach
that guard without a correctness or resource failure.

`PROTOCOL.md` and `protocol.json` now explicitly allow and record this exact
source-defined fallback while still requiring direct unpacking and the same
fixed environment. Every original cell remains in the grid. The outer-time
definition, optimistic ceiling, four primary groups, 2x upper-bound stop
rule, caps, discovery seeds and untouched holdout seeds are unchanged. No
performance sample from the registered grid preceded this amendment.

The measurement text also now states explicitly that exact returned-row
digest validation and destruction are charged inside outer time `T`, while
immutable input construction and the independent small F4 cross-check are
common supplied-reference work outside the call clock. This clarifies the
implemented metric without changing the ceiling formula or selecting data.
