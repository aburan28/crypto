# Prospective N83 root11 capacity-edge holdout

The 24 public N83 points in `targets.jsonl` were sealed before any
K3088/K3104/K3114/K3136 producer outcome. Their SHA-256 is
`725a36ccc58ecbe2904faeaf59a492ecb0acb58b9078d6589ab8e21efb0e404e`;
the verification-only answer sidecar SHA-256 is
`f3c74059d1fda5335ac5dbd50859d0d33f9daaaed660b18788841acc25c593a3`.
All 24 points are distinct from prior N83 study fixtures checked at sealing.
The exact field, curve, generator and prime subgroup are identified by
`EC1N83Ckb1a1haacb69ee9f37` and the shared canonical curve record.

## Frozen arm rule

K3072/root11 is the fixed control, with candidate ID
`IC1N83Ckb1a1fb509952PDP4rootRCguidedLAgaussTDdirectISO0hec7ae07d8b9a`.
Among independently replayed, completed K3088, K3104, K3114 and K3136
arms in the sealed 128-point screen, require at least two complete arms and
at least one geometric target-online mean below the frozen K3072 mean. Select
the two completed arms with the lowest geometric means, breaking exact ties
by smaller K. If the gate fails, leave selection null and do not dispatch.
New candidate IDs remain null until their exact base dumps and executable
configuration are materialized and checked. `selection_plan.json` and
`choose.py` freeze this rule before any new screen outcome.

## Tier 1 comparison after arm materialization

Assign points 0–5, 6–11, 12–17 and 18–23 to four independent cloud
instances or Pods. Run the three selected exact IC1 candidates on each host's
six points in fresh processes with empty base, root-index and factor-log
caches. Pair every IC run with a separately seeded strong rho solve of the
same public point. Rotate candidate order by point index, alternate IC/rho
order by point and candidate, and run an exposed-point rho calibration before
and after each pair. Freeze the source, binaries, input and commands, and
retain topology, affinity, cgroup limits, throttling, steal, frequency,
memory and raw output. Apply 900-second IC, 120-second rho and calibration,
300-minute producer and 330-minute host job caps.

Admit each stratum only after exact base, rank, four-point relation and
scalar replay, five target-phase reconciliation, same-point rho scalar
replay, source and resource checks, calibration noise at most ten percent,
and provider host identity. Preserve failed or noisy pairs as censored rows.
For admitted strata, compare the 24 paired log online intervals for every
candidate pair using an exact two-sided sign-flip test and Bonferroni
threshold `0.05/3`. Select an online winner only when its observed mean is
lower than both other arms and both tests pass. Apply the same rule separately
to complete cold-process intervals. Report rho/IC online ratios on each
identical public point and the observed Tier 1 envelope.
