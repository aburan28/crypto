# Durable compute for the a=0 full-subgroup native Jacobian campaign

- **Status:** investigating
- **Kind:** measured work-rate and persistence/resource gap; not a proved obstruction
- **Affected work:** [full-subgroup E/Jacobian comparison](../research/binary_genus2_speed_20261014/PROTOCOL.md)
- **Next action owner:** current native campaign
- **Last updated:** 2026-10-14

## Problem

Completed full-order timing for the pinned a=0 native genus-two image remains
outstanding: its prime order is 2417851639230796216685689, and the observed
checked-Cantor work rate implies a months-scale heuristic on one continuously
occupied core, while the available current worker has no durable-retention
guarantee. A running checkpointed search is not a completed solve.

## Evidence and uncertainty

- **Observed:** the public-job campaign is executing both GF(2^83) subgroups,
  retaining full collision tables and lane state. The a=0 five-round prefix
  median is 10,677.876657510944 steps/s for checked Cantor plus image
  conversions; source and raw records are in the linked study.
- **Derived, heuristic:** sqrt(pi*r/(2*166)) is about 151,258,823,772 walk
  steps. Dividing it by that recorded prefix rate gives about 14,165,628
  seconds before scheduling and other overhead. This is neither a measured
  solve duration nor a proved lower bound or reliable ETA.
- **Available:** executable subgroup-preserving transfer, native arithmetic,
  atomic checkpoints, restart replay, independent scalar verification and
  an active single-worker search. The inverse-image-to-E route is separately
  implemented and should be assessed as such.
- **Missing:** completed a=0 native-Jacobian solves, guaranteed long-lived
  compute/checkpoint retention, complete common-cost pricing, and independent
  host replay. No new paid infrastructure has been provisioned.
- **Not established:** impossibility of recovery, security, a genus-two
  index-calculus cost, or that optimized native formulas would have this rate.

## What would unblock it

- [ ] Preserve and resume the exact public-job checkpoints on compute with
  sufficient retention, or independently validate a separately named improved
  native arithmetic/solver implementation on the same full-order workload.
- [ ] Complete and independently verify the requested full-subgroup targets,
  with setup, transfer, native search, collision, memory, checkpoint and
  verification costs retained. Do not replace them with intervals or prefixes.

These conditions establish the missing executions. Admissible isolation,
common pricing and independent-host replay remain separate comparison gates.

## Next discriminating step

Finish the first a=1 native full solve as a correctness and work-rate control,
retain a=0 progress, and keep the uncapped current-worker campaign running
with atomic checkpoints. Publish each later snapshot in a new directory.
An actual a=0 verified collision establishes a completed target; continued
progress is a censored observation, not a negative result. An infrastructure
loss calls for resumption from the last checkpoint, not a silent restart.
Recheck when an actual solve, a different measured native method, or durable
compute access becomes available. The original workload remains open.

## Resolution and history

Unresolved. The named blocker documents the required continuation; recording
or merging it does not complete the campaign.
