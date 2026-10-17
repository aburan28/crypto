# Frozen non-isomorphic isogeny and rational torsion-fiber control

Date: 2026-10-13 local. Conductor T-2405. Commit this protocol before
executing the native search or solver. Follow the prior
[same-target protocol](../paired_round/PROTOCOL.md), its
[isomorphism control](../paired_round/iso_control_20261013/REPORT.md), and
the repository's complete-cost and evidence gates.

## Question and exact candidate

Test whether an explicit degree-three isogeny changes complete *relation
solver* work on the same source targets when its destination is not
base-field isomorphic to its source. Also exercise a nontrivial rational
torsion fiber. The candidate is the first prime `p=6r-1` with prime
`r>=101`; direct trial division through 1009 must confirm `p=617, r=103`.
Use `E/F_617: y²=x³+1` and `E'/F_617: y²=x³-27=x³+590`.
The proposed map is

`phi(O)=O`, `phi(0,±1)=O`, and for `x != 0`,
`phi(x,y)=((x³+4)/x², y(x³-8)/x³)`.

Before any solver run, native code must verify primality of 617 and 103,
enumerate both curves and confirm their order 618, enumerate the map on
every rational source point, check its group homomorphism on every ordered
pair, and confirm rational kernel exactly `{O,(0,1),(0,616)}`. It must
enumerate sixth powers in `F_617^*` and confirm `590` is not one, which
proves these `j=0` models are not base-field isomorphic via a scaling.
Any failure blocks the solve; preserve the failed output.

## Subgroup, targets, and paired arms

The target subgroup has prime order 103 and cofactor 6. Choose the first
affine source point in `(x,y)` order for which `[6]P != O`; set `G=[6]P`,
and verify `[103]G=O`. Set `G'=phi(G)` and verify its order 103. Generate
eight distinct nonidentity source targets by SplitMix64 initialized at
`0x2026101300010000 + 617`: draw `x=next()%617`, lift the smaller square
root, flip its sign if the next word is odd, project by `[6]`, reject `O`
and duplicates. Freeze their coordinates in the raw run and confirm them
independently; no generating scalar is supplied to a solver.

Each target receives three same-target arms: source relations on `(E,G,Q)`,
degree-three destination relations on `(E',G',phi(Q))`, and source negation
rho on `(E,G,Q)`. The source and destination relation solvers use the same
rank, pair-table, random-stream, 10,000-trial collection/descent caps, and
scalar-recovery rules as the prior paired round. In each factor base, scan
abscissae upward and use the smaller square root; project the lifted point
by `[6]`, skipping infinity and duplicate subgroup points, until twelve
points are retained. Charge every lift, projection, pair-table operation,
relation trial, rank update, recovery, descent, and verification. Seeds are
`0x52454c4154494f4e ^ p ^ (index<<16)` for relations,
`0x44455343454e5400 ^ p ^ (index<<16)` for descent, and
`0x52484f4e45470000 ^ p ^ (index<<16)` for rho. Rho uses the prior
negation-class walk and a 100,000-step cap; a six-class rational
automorphism control is not available for this `p=617` pair because
`F_617` has no nontrivial cube root of unity.

Every arm must recover and independently verify the same scalar on the
original source `(G,Q)`. The destination arm also verifies it on
`(G',phi(Q))`. Retain capped or failed rows and block their paired cost
comparison. A change in cap/restart policy requires a separately frozen
version and retained earlier run.

## Boundary, accounting, and acceptance

The exact-target rho reference is the charged negation-class arm; the
generic two-class walk-step floor is `sqrt(pi*r/4)`. Report
`S=(counted group additions + doublings)/sqrt(103)` for complete counted
phases, and the ratio to matched rho and to the floor in one table. This
is a group-operation diagnostic: separately record map multiplications,
squarings and inversions, factor-base square roots, modular row work,
lookups, rho hashes/probes, memory, construction/preflight, and failed
attempts. Do not infer full native-cost speed or a scaling exponent from
this one subgroup order. Record cold setup per target and do not amortize
it across eight targets in the primary comparison.

Acceptance requires all 24 arm records verified, a byte-identical same-host
native replay, and an independent Rust checker using its own point law and
inverse for curve order, map, subgroup, targets, recovered scalars, and
phase totals. Retain source/binary/input hashes, exact commands, host and
isolation, raw records, phase table, run report, graph, and PDF. If this
bounded control succeeds, the next performance gate is fully priced native
work at four growing sizes, five interleaved rounds, and independent-host
replay. A failed search or timeout is bounded evidence only.
