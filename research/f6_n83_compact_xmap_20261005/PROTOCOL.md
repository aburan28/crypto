# n83 F6 signed-pair compact x-key lookup gate

Registered before changing code or measuring the candidate. The retained
compact pair-sum branch #1416 stores each full sum in 64 bytes, but its
4,108,723-entry lookup still uses the 32-byte `PointXKey` enum as a hash
key. Replace that map with raw `u128` affine x keys and a separate
identity slot. This preserves exact identity handling even on fields
where all 128 x-coordinate bits are used. Keep the existing random-state
hashing, insertion order, target-query order, sign check, and independent
full-group witness replay. This is a four-summand F6-IC component only.

Baseline: #1416 head `ee93b179a` and its frozen binaries/receipts in
`research/f6_n83_compact_pairs_20261005`. Freeze the registered K0
curve `icv1-f2m83-tm6151469093347-debefd74`, standard cofactor-projected
dimension-8/10/12 bases (258, 1,048, 4,054 actual usable points), and
the same public T001 point
`(355fb5df7a905f16921eb,5900a390f42d290f1bbe)`.
The full base has 8,219,485 unordered pairs and 4,108,723 signed-sum
representatives. Build matched release binaries with the same compiler
and flags, then freeze their SHA-256 values before running the panel.

Pass the nine focused geometry tests, including exhaustive small-base
four-sum membership, identity and same-x exceptions. Run the full-base
planted `[0,2,4,6]` control through both portable and PMULL queries and
replay any witness in the group. Require identical representative counts
and the same exact ordinary T001 outcome. Run the small-base three-repeat
probe once per arm, then full-base baseline, candidate, candidate,
baseline, with no build between arms. Each process has a 120-second
timeout. Preserve stdout, stderr, exit codes, code/input/binary hashes,
index-build and query intervals, and peak RSS.

Retain the raw-key map only if exactness passes, full-base median query
time falls at least 10%, full-base peak RSS falls at least 10%, and
full-base median build time does not rise more than 10%. Otherwise
restore the baseline runtime code and archive the rejected patch and
raw results. This host lacks auditable exclusive CPU isolation, so wall
ratios are exploratory diagnostics. Index construction is reusable
target-independent preparation. The K0 four-summand natural-target
coverage ceiling is `4.662e-12`; no complete F6 or IC/rho speedup can be
inferred from this gate.
