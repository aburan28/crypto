# Prime61 exhaustive horizontal class exceeds the current in-memory walker

- **Status:** open
- **Kind:** resource limit and implementation gap; full class-count certification gap
- **Affected work:** [prime61 experiment](../research/prime61_isogeny_20261016/README.md), [follow-up card](../backlog/XP-20261016-prime61-horizontal-class/card.json)
- **Next action owner:** unassigned native isogeny-walker lane
- **Last updated:** 2026-10-16

## Problem

The requested exhaustive horizontal-class experiment at `p=2^61-1` cannot
be executed with the existing all-in-memory walker in this 16 GiB
container. The source is `a=p-3,b=111`, prime order
`2305843010818082053`, with exact G,Q and measured EC1 in the linked report.
Fifteen outgoing prime-degree kernel edges through 61 were tested; that
one-hop panel does not establish exhaustive class coverage.

## Evidence and uncertainty

- **Observed:** native walk verification passes 16 curves and 15 edges;
  192 capped streams yield zero relations and solves. PARI full
  `bnfcertify` times out at 300 seconds (exit 124).
- **Derived:** discriminant `-6649310858224389603` is fundamental, so all
  base-field orders have conductor one and no vertical base-field edges
  exist. The [exact reduced-form certificate](../research/prime61_isogeny_20261016/results/class-order-certificate.txt)
  proves one class has order 68,162,026; the standard ordinary CM class
  action therefore implies at least this many vertices.
- **Resource bound:** native [header sizes](../research/prime61_isogeny_20261016/results/min-class-node-edge-memory.json)
  are 240 bytes per Node and 152 per Edge. Even a spanning record requires
  26,719,514,040 bytes before heap allocations, versus the 17,179,869,184
  byte limit. This applies to the existing retained-node/edge backend.
- **Conditional:** three native class-group routines agree on 545,296,208,
  with group invariants `[68162026,2,2,2]`; full certification is unfinished.
  That larger count is not used to prove the resource blocker.
- **Missing capability:** external-memory or streaming traversal with a
  verifiable coverage certificate and resumable per-vertex workloads, or
  an adequately provisioned implementation after exact sizing.
- **Not established:** exhaustive reachability for the chosen degree set,
  a certified full class count, any full-class ranking, full-cost gain, or
  mathematical impossibility of enumeration with a different backend.

## What would unblock it

- [ ] Produce a certified complete class structure or independent exhaustive
  coverage certificate for this exact discriminant and trace.
- [ ] Implement bounded-memory, checkpointed enumeration with exact model
  deduplication, ordered kernel routes and generator/target transport.
- [ ] Freeze resources and protocol, then retain a verified workload row
  for every class vertex, including caps and failures, with independent replay.

These conditions would discharge the exhaustive-class interpretation.
They would not by themselves reproduce the screenshot's unlocated
candidate search or establish an operation-cost advantage.

## Next discriminating step

Implement and validate an external-memory walker on small classes whose
complete native enumeration is already known; compare vertex and edge
sets and all transported targets, then measure retained memory on a
larger checkpointed prefix. Freeze a new resource protocol before an
exhaustive run. A passing prefix enables sizing, not full coverage; a
set mismatch blocks rollout. A timeout remains a bounded failure. Recheck
when that backend or a certified full class structure becomes available.

## Resolution and history

2026-10-16: unresolved. The first node-only lower bound did not exceed
16 GiB; adding the mandatory connecting Edge headers does. Both records
are preserved. No exhaustive traversal or paid remote run was started.
