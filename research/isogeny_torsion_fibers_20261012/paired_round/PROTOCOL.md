# Same-target isogeny transport and relation-solver test

Date: 2026-10-12. Conductor task: T-2223. This protocol is frozen before
implementation and execution. It follows the structural control in the parent
directory and uses public synthetic instances only.

## Question and boundary

Does changing coordinates through the explicit degree-three isogeny change the
cost of a complete relation-based discrete-log solve on the same source target?
The generic reference is the repository's `rho.negation` walk on the source
curve and exact target. The generic walk-step floor is
`sqrt(pi*r/4)` for negation classes; the full solve is reported separately.
No performance or scaling claim follows from this small correctness round.

## Frozen curves and targets

Use the existing exact pairs `p=13, b=2, r=19` and `p=199, b=3, r=211`,
both with `a=0`, source `y²=x³+b`, destination
`y²=x³+(-27b mod p)`, and geometric kernel
`{O,(0,+sqrt(b)),(0,-sqrt(b))}`. In both cases `b` is nonsquare, the
rational kernel has order one, and the curve orders are prime `r` per the
parent census. Recount both sides and verify the map on all rational points
before solving. Pick the source generator as the first affine point in
lexicographic `(x,y)` order; transport it by the map and check order `r`.

Generate eight distinct nonidentity public source targets per curve by a
fixed SplitMix64 stream with initial state `0x2026101200010000 + p`:
draw `x=next()%p`, lift using the smaller square root, flip its sign when
the next stream word is odd, reject an unliftable `x` or a repeated point,
and continue. No scalar used to generate a target is supplied to any solver.
Transport each exact `Q` to `phi(Q)`; rho receives source `(G,Q)`, the
source relation solver receives `(G,Q)`, and the destination solver receives
`(phi(G),phi(Q))`. The runner verifies all returned logs on the source and
destination independently and requires equality modulo `r`.

## Frozen relation path and caps

For each solver arm and target, build a fresh factor base consisting of the
first `B` abscissae that lift, with the smaller root, where `B=4` at `p=13`
and `B=12` at `p=199`. Build every unordered pair sum `P_i+P_j`, `i<=j`,
and store all representations per sum. For each relation trial, draw a fresh
nonzero `k mod r` from a separate SplitMix64 stream seeded by
`0x52454c4154494f4e ^ p ^ (target_index<<16)`, compute `[k]G` through
the counted group, and emit every pair representation with that sum as
`log(P_i)+log(P_j)=k mod r`. Incremental modular elimination retains only
independent rows. Stop collection at full rank or after 10,000 trials.
Recover factor-base logs by back substitution, then verify every log by
scalar multiplication. For target descent, draw `t` from a separate stream
seeded `0x44455343454e5400 ^ p ^ (target_index<<16)`, compute
`Q+[t]G`, and look for a pair representation. Return
`log(P_i)+log(P_j)-t mod r`; stop after 10,000 descent trials. Verify
the candidate on both curves. All streams and caps are unchanged by outcome.

Run the repository's existing `rho_walk_with(..., NegationClasses,
RhoWalk::negation())` with seed
`0x52484f4e45470000 ^ p ^ (target_index<<16)` and a 100,000-step cap.
Retain exhausted or failed rows. This is the matched source-target reference,
but the custom pair runner is not an admitted ecbench session.

## Accounting and acceptance

Record separately: exact model/map checks; map evaluations and their field
operations; factor-base lift and pair-table setup; relation attempts,
accepted independent rows and lookups; modular row operations; descent
attempts and lookups; verification; rho's group and field operations and
its uncharged counters. Every group add, double and scalar multiply goes
through `CountedGroup`; prime group field operations come from that ledger.
The map's explicit field work, square roots, lookups and row arithmetic
remain separate counters. Since there is no pinned conversion for all these
native operations, group-op `S` is diagnostic and no end-to-end speed ratio
will be promoted. Host wall times are descriptive only on this shared macOS
host, with no isolation admission.

Success requires every arm to return the same verified logarithm for all
16 workloads, exact source/destination map checks, and deterministic
independent replay from the frozen executable and inputs. A correctness
failure, cap exhaustion, or replay difference is retained and blocks that
workload's comparison. Reproduce once by rerunning the executable and
comparing byte-for-byte outputs; this is same-host replay, so independent
host validation remains an open later gate. Do not infer a general advantage
or an exponent change from this round.
