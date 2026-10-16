# N83 K3114 exact-root online frontier screen

The [Lean capacity result](../n83_online_finegrid_20261016/README.md)
identifies K3114 as the final folded-column count using the 2³⁰-slot,
five-byte dynamic root table. This screen measures its one-target online
interval against K3072 on the same previously exposed public point. Under
the current 166-points-per-column construction, their expected usable
factor-base sizes are 516,924 and 509,952 points. The actual sets, base
hashes, index-entry counts and effective columns must come from the producer
and independent verifier.

The exact producer source SHA-256 is
`3ee6bc11399b382e72789b152043e5a67abc5b783c4e998ad8c534d60baebc03`;
the frozen point, curve, subgroup, rho seed, source files and lockfile are
checked by `run.py --preflight`. The arm order is K3072 A/B, K3114 A/B,
then one same-point rho solve, all on one GitHub-hosted Ubuntu 24.04 runner.
Every IC arm is a fresh process with an empty base and index, 900-second
cap and 15.25 GiB process RSS guard. Preserve cap, failure and resource
outcomes as rows. The rho arm uses the same public point and its frozen
single-target walk configuration.

Before interpreting online timing, independently replay each base, target
relation, scalar and rho answer. Require both repetitions of each K to
have identical base hashes and exact index-entry counts, admitted resource
receipts, and online intervals differing by at most 10% of that pair's
mean. Select a faster K from this screen only if its slower interval is
faster than the other K's faster interval by more than the largest observed
absolute A/A difference. An incomplete, noisy or tied comparison retains
`BOUNDED_UNKNOWN` and a null selection. A selected row is an exposed-point
screen result; the separately frozen five fresh-point Tier 1 panel remains
the single-target speedup gate.

The workflow archives the full source/host receipts, binaries, raw producer
rows, usage, base dumps, verifier analysis, stdout and stderr, including
timeouts and failed arms. The downloaded artifact must be hashed and
independently replayed before exact IC1 candidate manifests or results are
versioned. Treat cloud-runner timing as exploratory until its full Tier 1
telemetry and same-point panel gates pass.
