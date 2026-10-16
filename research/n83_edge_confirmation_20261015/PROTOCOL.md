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
Require every final tournament arm to pass exact base, four-point relation and
scalar replay under its 14.5 GiB guard, with zero new throttled periods and OOM
kills. If that set is incomplete, record `INCOMPLETE` and do not select. An
edge candidate is eligible only when both of its mirrored arms pass the same
checks under the 14 GiB edge guard. Preserve every failed edge arm; it does
not prevent another complete edge candidate or a confirmation among
unresolved tournament finalists.

When K3450/root11/wide3 has both eligible edge arms, let the challenger be the
lower mirrored mean among the eligible K3500/root11/wide3 and
K3550/root11/wide3 candidates; ties break toward lower K. Admit it only when
its mirrored mean is strictly lower than K3450/root11/wide3 in that edge
screen. With no eligible comparator or challenger, admit no edge candidate.
Let the final tournament's descriptive lowest mean be the leader. Admit that
leader and every other finalist whose 99% interval for `leader minus other`
has an upper endpoint greater than or equal to zero. If the tournament
selected a pipeline, this rule admits only that pipeline. Append the edge
challenger if one qualifies, then deduplicate exact candidate identities.
Keep the original tournament order and place a new edge challenger last. If
fewer than two distinct pipelines remain, record `NO_OPEN_COMPARISON` and stop.

## Confirmatory measurement

Freeze exact candidate IDs, manifests, source hashes, base hashes, index
sizes, target and sidecar hashes, a 1,500-second per-arm cap, 330-minute job
cap and 14.5 GiB RSS guard before running. Execute the admitted complete
pipelines in list and reverse-list order on the same cloud allocation and all
4,000 fresh points. Use fresh processes with empty indexes and logs, identical
source and resource limits. Archive source, binary, host, cgroup,
calibration, usage, full target output and base dumps for every arm. Account
for setup and final relation linear algebra separately from each target's
online interval. Preserve failed, timed-out and over-cap arms.

Independently replay each subgroup point, base digest, four-point relation,
scalar, target-online phase sum and source/host receipt. A comparison is
eligible only if every arm passes under its frozen RSS guard with zero new
throttled periods and OOM kills. Require before and after SHA-256 calibration
triples to match, each triple's max/min wall time at most 1.15, and after/before
median ratio between 1/1.15 and 1.15. Use the target-wise average of the two
repetitions for each pipeline, then 80 contiguous blocks of 50 public points
and 10,000 block bootstrap resamples with seed `20261015601 + 16*i + j` for
the ordered pair of candidate positions `i,j`. The lowest mirrored mean wins
only when the upper endpoint of its two-sided 99% interval against every
other admitted pipeline is below zero. Otherwise record `BOUNDED_UNKNOWN`
with every measured interval. This is a secondary factor-base selection panel;
a separate Tier 1 same-point IC/rho result is required for the primary
one-target online comparison.
