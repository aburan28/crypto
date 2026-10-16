# N83 K1104/K1120 cold precision panel

The five-point Tier 1 holdout measured K1120 0.401 s faster than K1104, with
paired 95% interval [-0.236, 1.038] s. This panel fixes 16 new public points
and a mirrored fresh-process schedule to resolve that comparison. The panel
belongs to the same Koblitz `GF(2^83)` curve, root11 dynamic index layout,
four-summand decomposition, factor-base policy, and source-bound producer as
the five-point holdout. It compares exact K1104 and K1120 candidate manifests
from that holdout. No K value or point is changed after this protocol is
committed.

## Frozen input and schedule

Generate indices 0 through 15 using
`s_i = 1 + (little_endian(SHA256(domain || u32_le(i))) mod (r-1))`,
`domain = n83-e2e-cold-precision-sixteen-20261016-v1`, and
`r = 8569786107849059`. The producer receives only `Q_i = s_i G`; the
verifier-only answer sidecar supplies `s_i` after execution. The generator
replays all earlier committed N83 public-point files, checks subgroup
membership and distinctness, and rejects any repeated point. Commit the
point file and sidecar before timing. Retain every point and failed arm.

On each point, run four fresh IC processes with empty base, index, and log
caches. Odd points use `K1104a, K1120a, K1120b, K1104b`; even points reverse
that order. Run one same-point rho solve, before IC on even points and after
IC on odd points. Use the same exact source, compiled binaries, caps,
verification, Tier 1 host receipts, and resource envelope as the sealed
five-point holdout. Each IC arm has a 480 s cap; rho has a 120 s cap. Require
at least four exposed CPUs, 15.5 GiB memory, and 2 GiB free disk.

## Frozen decision

For each point, average the two K1104 cold process times and separately
average the two K1120 times. Define paired saving as K1104 minus K1120.
Compute the mean and two-sided 95% Student-t interval over all 16 paired
savings (`t_(0.975,15) = 2.131449545559323`). Select K1120 within this
two-candidate comparison only if the lower bound is strictly positive;
select K1104 only if the upper bound is strictly negative. Otherwise record
`BOUNDED_UNKNOWN` and a null selected candidate. Both candidates must pass
exact-base, relation, scalar, source, resource, topology, order, and replay
checks on every arm. Each candidate's A/B cold spread must be at most 2% of
its A/B mean on every point. A failed, missing, noisy, or timed-out arm makes
the aggregate selector null while preserving all rows and statuses.

Report one-target cold process and target-online intervals separately, with
same-point rho online timing and complete phase costs. This panel identifies
the faster of these two exact factor bases under the measured Tier 1 envelope;
other geometries and index layouts require their own paired comparisons.
