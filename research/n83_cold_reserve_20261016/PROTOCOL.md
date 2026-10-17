# N83 four-point cold factor-base confirmation

This panel compares four exact factor-base pipelines selected before any
reserve-point timing: K1024/root22 and K1536/root22 from the repeated cold
screen, plus K1152/root22 and K1152/dynamic root11 from the completed
neighborhood screen. The four public targets and answer sidecars were frozen
before either screen. Their SHA-256 values are respectively
`6ab08db41c5e86c96162f538a0c45d53331c2f9f823a66eea9bf0b458dd3554d`
and `b21d3a1d54c701255a079a3f1d88940f46563a2f857303b4e3bb98ac40483f39`.
All four have distinct points and scalars and are outside the earlier public
target manifest. Each producer receives a one-line public-point file, never
the answer sidecar. Each arm starts a new process with an empty factor base,
exact root index and factor-log system.

The curve is `EC1N83Ckb1a1haacb69ee9f37`, polynomial-basis `GF(2^83)`
with modulus bitmask `9671406556917033397649557`, model
`y²+xy=x³+x²+1`, subgroup order `8569786107849059`, cofactor
`1128547018`, and rank seed 29. The exact IC producer source SHA-256 is
`d094345474605075566d646ddecef8f694cdd644b0f9190a6b1ef1551fd06888`.
All layouts use dynamic exact five-byte root slots, symmetric S3 and sparse
canonicalization without a Bloom filter. K1024, K1152 and K1536 contain
169,984, 191,232 and 254,976 distinct subgroup-usable base points before
orbit folding, and the corresponding effective column counts are K.

Run five fresh IC processes per point on one Ubuntu 24.04 cloud allocation:

| Point | First control | Middle arm 1 | Middle arm 2 | Middle arm 3 | Last control | Rho order |
| ---: | --- | --- | --- | --- | --- | --- |
| 1 | K1024/root22 A | K1152/root22 | K1152/root11 | K1536/root22 | K1024/root22 B | after IC |
| 2 | K1024/root22 A | K1152/root11 | K1536/root22 | K1152/root22 | K1024/root22 B | before IC |
| 3 | K1024/root22 A | K1536/root22 | K1152/root22 | K1152/root11 | K1024/root22 B | after IC |
| 4 | K1024/root22 A | K1152/root22 | K1536/root22 | K1152/root11 | K1024/root22 B | before IC |

Each IC child has a 900 s cap; each fresh strong-rho process has a 120 s
cap and a point-specific seed `2026101601 + point_index`, with one-based
point index. The allocation must expose at least four CPUs, 15.5 GiB host
memory and 2 GiB free disk. Require peak RSS at most 15.25 GiB and record
source/binary hashes, topology, affinity, cgroup limits, per-process outer
wall and CPU time, major faults, throttle/OOM deltas and all failures.
Independently replay every completed base, four-point relation, IC scalar
and rho scalar against the frozen sidecar. Preserve capped and invalid rows.

The cold metric is producer `process_total`: first setup operation through
producer-verified one-target answer. Retain exclusive setup, rank, target
and five target-online phase costs, plus outer wall and independent replay.
For each point, use the mean of K1024 A/B as its K1024 time. Admit a cold
ranking only if all 20 IC and four rho runs replay, every resource gate
passes, each control has identical base hash and exact root-entry count, and
each A/B cold-process difference is at most 10% of its mean. Do not replace
an excluded point or arm after seeing its timing; a whole-panel retry uses
a fresh recorded allocation and retains the first panel.

Rank the four pipelines by their arithmetic mean cold process over the four
points. For every pair, retain the four paired time differences, their mean,
sample standard deviation and a two-sided 95% Student-t interval using
`t_(0.975,3) = 3.182446305284263`. A winner is selected only if it has the
lowest mean and, against each other pipeline, its advantage is positive on
all four points, exceeds half the observed K1024 A/B spread on each point,
and has a positive lower 95% paired bound. Otherwise keep the selected
candidate null and carry the full paired panel into a larger fresh-point
comparison. The A/A gate and paired interval are descriptive controls for
this cloud allocation, not a global optimum proof.

After a candidate is selected, generate and freeze at least five new distinct
public points for a separate alternating same-point IC/rho Tier 1 online
panel. The present four-point cold confirmation does not stand in for that
timing gate.
