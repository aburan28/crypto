# Fresh N83 edge confirmation

This protocol and its public-point fixture are fixed before the root11
K3450/K3500/K3550 edge screen reports. The 4,000 points use
`n83-edge-confirm-fresh-4000-20261015-v1` and are excluded from every earlier
N83 screen, the final tournament and the five Tier 1 same-point targets.
`answers.jsonl` is only a fixture sidecar for independent scalar replay.
The generator SHA-256 is
`eed2e5994628da77ab501597e7541f6d8d47fe0502edadd6ead4f1afe66087c9`;
the public points SHA-256 is
`66c4683e252f1cd321b4b6c7c9d06744c82140222ddce160faa4e7a0e0984b7e`;
the sidecar SHA-256 is
`e1207d1be9a7689c42209ed50f048325812436d5f7e6f313b5e5e8a04a5730b1`.

## Admission

Read the independently replayed edge screen and final tournament artifacts.
Require every edge arm to pass exact base, four-point relation and scalar
replay under its 14 GiB guard, with zero new throttled periods and OOM kills.
Require every final tournament arm to pass the same checks under its 14.5 GiB
guard. If either set is incomplete, record `INCOMPLETE` and do not select.

Let the edge challenger be the lower mirrored mean of K3500/root11/wide3
and K3550/root11/wide3; ties break toward lower K. Admit it only when its
mirrored mean is strictly lower than K3450/root11/wide3 in that edge screen.
Let the tournament comparator be its 99% selected pipeline, or, when that
decision is `BOUNDED_UNKNOWN`, its descriptive lowest eligible pipeline.
Deduplicate an exact candidate identity. If no edge challenger qualifies or
both identities are the same, record `NO_NEW_CHALLENGER` and stop.

## Confirmatory measurement

Freeze exact candidate IDs, manifests, source hashes, base hashes, index
sizes, target and sidecar hashes, child caps and RSS guard before running.
Execute both complete pipelines in A/B/B/A order on the same cloud allocation
and all 4,000 fresh points. Use fresh processes with empty indexes and logs,
identical source and resource limits. Archive source, binary, host, cgroup,
calibration, usage, full target output and base dumps for every arm. Account
for setup and final relation linear algebra separately from each target's
online interval. Preserve failed, timed-out and over-cap arms.

Independently replay each subgroup point, base digest, four-point relation,
scalar, target-online phase sum and source/host receipt. A comparison is
eligible only if every arm passes under its frozen RSS guard with zero new
throttled periods and OOM kills. Use the target-wise average of the two
repetitions for each pipeline, then 80 contiguous blocks of 50 public points
and 10,000 block bootstrap resamples with seed 20261015601. The challenger
advances only if the upper endpoint of the two-sided 99% interval for
`challenger minus comparator` is below zero. Otherwise record
`BOUNDED_UNKNOWN` with the measured interval. This is a secondary factor-base
selection panel; a separate Tier 1 same-point IC/rho result is required for
the primary one-target online comparison.
