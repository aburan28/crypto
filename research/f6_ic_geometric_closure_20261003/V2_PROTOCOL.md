# F6-IC v2: packed exact geometry, preregistered before execution

V1 cut total inherited-F4 reductions from 379 to 163 on its new public
target, but used 17,859 general-representation point additions and changed
the full online call only slightly. Version 2 keeps the **same** exact
one-fixed branch logic and 256-point cap, replacing its inner arithmetic
with the repository's independently tested single-word `FastCurve` on fields
it supports. At gate creation, build a target-dependent packed point vector,
exact packed point-to-index map and packed target. Charge this construction
to target PDP. Each branch adds packed points, negates and subtracts the
target in packed form, and performs an exact packed residual lookup. It still
checks partial assigned x-code bits and independently replays a proposed
relation with the general group arithmetic. If the packed representation is
unavailable, use v1's general path. The source, candidate manifest and code
digest change; v0 and v1 receipts stay immutable.

Before timing, compare packed and general branch decisions on exhaustive n9
positive and negative cases, including both partial-bit values, permuted
coordinates, missing support and affine placeholders. Run the prepared n17
golden scalar control. For the bounded pilot, use the same certified
target-independent n17 state, 62 actual usable points, 29 folded columns,
`m=3`, degree 3, 8192 nodes, one Rayon worker, one **new** public point from
fixture seed `20261003021`, algorithm seed `20261003036`, one query per trial,
and at most 32 queries. Freeze the exact point, source/binary/input hashes,
candidate IDs, workload ID and resource envelope in a commit before timing.
Run six fresh processes under 600 seconds each in order F4/F6, F6/F4,
F4/F6. Preserve all exits and raw phase/attempt ledgers.

The correctness gate requires all six verified calls to recover the same
scalar, all five phase sums to close, and packed and general small-base
decisions to agree. The mechanistic gate requires the new fast path to be
used, retain fewer F4 reductions than inherited F4, and charge every packed
addition and lookup. The exploratory engineering target is at least 2x lower
complete online time in each of the three paired runs. Any controlled CPU
gain needs the isolated benchmark service; a rho comparison additionally
needs a sealed `ecbench` session against strong rho on the same public point.
If full-call work regresses or the arithmetic does not produce a material
gain, preserve the negative receipt and stop this branch of F6.
