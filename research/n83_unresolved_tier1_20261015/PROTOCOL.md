# N83 one-target Tier 1 diagnostic for both unresolved finalists

The final tournament and the fresh confirmation each placed
K3425/root4/standard ahead of K3450/root11/wide3 by mirrored mean, while
their frozen two-sided 99% paired intervals crossed zero. Their verified
confirmation selection is `BOUNDED_UNKNOWN`. This diagnostic keeps that
selection status and asks a separate primary question: how quickly does each
complete pipeline solve the **same one previously unseen public point**
online relative to a matched strong-rho solve on one Tier 1 cloud instance?
The two IC times also form a descriptive same-point comparison; one target
does not establish a distribution-wide factor-base ranking.

The single public point was generated with domain
`n83-unresolved-tier1-one-fresh-target-20261015-v1` after excluding all
prior N83 panels, the tournament, confirmation, and the unused five-point
Tier 1 fixture. Its public-point SHA-256 is
`dfea95ee4e423dd627adbd116828e9f478779c51476f0fa172cf4c160e645d12`;
its verifier-only scalar sidecar SHA-256 is
`ef4d939a30fbf26809f97acc0607fcd4fe04661e658e8b353c1ec77f5075d335`.
The curve is `EC1N83Ckb1a1haacb69ee9f37`, with subgroup order
`8569786107849059`. The IC source SHA-256 is
`f2745227d95e36ef45c73c5104a23ada6cebb0794405368581d2be17044845dc`;
the strong-rho source SHA-256 is
`fb14072f4fa5d08e18d8b79130b6d58a220adae0f76369cc9be949691c83eaa0`.
The frozen finalist admission SHA-256 is
`803df5758dad2e7d4769554a6d6eb741c61ce7b8c8b1887c843f1616fcc9e1d8`;
the confirmation analysis SHA-256 is
`81a945b5711cf599996aa127aa3fbec3128ca8f523a8bad0b2fc8a4244e7359a`.

Use one cloud VM or Pod for both candidate pairs. Freeze its allocation,
topology, four execution CPUs, source and binary hashes, affinity, cgroup
limits, memory and disk. Require at least four visible CPUs, 14.5 GiB
observed RSS capacity plus 256 MiB host headroom, and 2 GiB free disk.
Each IC child has a 1,500-second cap; each rho child has a 120-second cap.
Run fresh processes and empty IC indexes and logs. Candidate order is
K3425/root4/standard then K3450/root11/wide3. For the first candidate run
`IC,rho`; for the second run `rho,IC`. Both pairs receive the identical
public point. Use rho seeds 2026101401 and 2026101402.

Before and after each pair, run three fixed one-million-step SHA-256
calibration samples with identical checksums. Require each triple's maximum
divided by minimum at most 1.15 and after/before median ratio between 0.85
and 1.15. Require zero new cgroup CPU-throttled periods and OOM kills,
observed peak RSS at most 14.5 GiB for each child, nonnegative system
counter deltas, and steal ticks no more than 1% of total CPU ticks over
the pair. Preserve every timeout, failure, noisy pair, raw target output,
base dump and calibration receipt. No target replacement is permitted.

Independently replay the public point, subgroup scalar, factor-base digest,
four-point relation, phase sum, rho scalar, source receipt, host and noise
envelope. Start each target-online IC clock at its first target-dependent
operation after reusable base/index/log setup and stop after scalar recovery
and producer verification. Start rho at its first target-dependent walk and
stop after scalar recovery and producer verification. Record independent
external replay separately. Report each verified same-point
`rho_online_ms / IC_online_ms`, all five exclusive IC online phases, cold
setup, RSS and resource counters in one table. A speedup for a candidate
requires both arms of its own pair to pass. Preserve the factor-base
ranking as `BOUNDED_UNKNOWN` regardless of this single-point observation.
