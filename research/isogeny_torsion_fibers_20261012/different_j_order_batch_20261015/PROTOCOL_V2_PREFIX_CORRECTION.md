# Correction before execution: cold costs at every target prefix

**Status:** frozen before target freeze and before any order-229 batch
measurement. This versioned correction supersedes the process-count,
scalar-replay-count and prefix-cost sentences in
[PROTOCOL_V1.md](PROTOCOL_V1.md). All other instance, target, base,
relation, control, cap, memory, seed, and verification choices remain
as preregistered there.

The V1 plan placed the 1, 4 and 16-target prefixes inside each
64-target child. Their internal timers would omit final output and
process shutdown, so they could not be called complete cold prefix
costs. The measurement will instead launch a **fresh child process for
each** prefix `N=1,4,16,64`, arm and round. That is
`4·7·3=84` metered child processes, each receiving the same committed
64-target public input and a frozen prefix length, then solving exactly
its first `N` targets. The output of one child has `N` target rows.
The aggregate is `7·3·(1+4+16+64)=1,785` independently replayed
scalar rows, including the V1-planned `1,344` at `N=64`.

For each child, use the native parent meter to record complete child
user+system CPU, wall, peak RSS, stdout/stderr bytes and exit status
from spawn through termination. The child also records common parsing,
construction, setup, every query and verification, and field/group
counters. Pair arms by exact `(prefix,round,target index)`. Interleave
arms by the V1 round rotation **within each prefix**; preserve the
prefix order `1,4,16,64`. Run the three-repetition, 100,000-addition
source-curve CPU/field calibration once before the 84 children, with
its exit, stdout and stderr retained. It is a unit-conversion
instrument and excluded from attack setup, while each child includes
its own factor-base or BSGS table setup.

The checker must require one result and meter receipt for each of the
84 `(prefix,round,arm)` cells, verify the exact prefix target sequence
and 1,785 recovered scalars, and independently reconstruct every
relation witness and rank-nine recovery. A capped or failed child
keeps its raw receipt and makes the corresponding paired cell
incomplete. Whole-child CPU ratios and any derived local
group-addition-equivalent `S` use the appropriate cold child totals;
the native internal checkpoints are secondary diagnostics. Portable
counted `S` remains bounded if non-group primitives lack measured
weights. The L0 shared-host wall and CPU observations do not become an
isolated performance claim from this correction.
