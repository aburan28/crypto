# Upper conductor-gap floor construction and complete cost receipts

- **Status:** investigating
- **Kind:** resource limit; implementation gap; proof gap (separate issues below)
- **Affected work:** [protocol](../research/conductor_gap_cm_orbits_20261015/PROTOCOL.md), [evidence](../research/conductor_gap_cm_orbits_20261015/REPORT.md)
- **Next action owner:** conductor-gap research workstream
- **Last updated:** 2026-10-15

## Problem

The native study cannot yet certify the prepared 95/96-bit floors, price every
native solver operation in one common cost unit, or supply a complete Lean CM
classification proof.

The rational-torsion family has explicit independently checked floors through
91 field bits. Its 95/96-bit inputs have rigorous prime and kernel certificates,
but the first large constructor was stopped after the 91-bit result when shared
memory reached 15,982,964,736 of 17,179,869,184 bytes. This is a resource limit,
not evidence that the floor does not exist.

The older near-limit candidate with prime conductor 129149576694539 is a
separate instance. Its Frobenius scalar does not provide the rational torsion
point required by this square-root Vélu implementation. A new rational-torsion
family endpoint does not construct that exact old candidate's floor.

## Evidence and uncertainty

Observed: all 100 old direct codomains match the fast moment constructor;
15 new floors have independent PARI order/prime/stratum receipts. See
`results/pair-certificates.jsonl`, `results/floor-certificates-v2.jsonl`,
`results/larger-construction-stop.json` in the study.

The native relation backend reports field/group work and queries, while modular
exponentiation, table/hash work, allocation and elimination still lack a complete
common price. A zero charged-GAE lower bound is not zero work. No eligible
IC/rho speedup is asserted. The five online phases include independent replay.

Lean proves conditional map-iteration equality and a norm-numerator bound.
The CM norm identification, endomorphism-order lattice, actual rational-map
homomorphism and order classification are unformalized bridges, listed in
`lean/README.md`. No impossibility or general ECDLP hardness gap is established.

## What would unblock it

- [ ] Complete a 95/96-bit rational-family codomain under recorded memory/time
  limits and independently certify its exact model, prime group order and floor
  conductor. Preserve failed or censored attempts.
- [ ] Construct and certify the floor for the original non-rational-kernel
  96-bit candidate using a method whose kernel prerequisites actually hold.
- [ ] Count and price every native component, then repeat complete same-point
  IC/rho windows with matching resource envelopes and earned isolation.
- [ ] Formalize the stated CM bridges in a pinned proof environment and attach
  checked proof receipts; arithmetic tests alone do not discharge them.

## Next discriminating step

Retry the version-2 precomputed-remainder-tree constructor after the primary
pair-table campaign releases memory. Use the original prepared inputs and a
recorded per-process memory cap. Success requires independent PARI checks;
a cap is a bounded negative result and leaves the exact floor requirement open.

For the older candidate, test the proposed construction's rational-kernel
prerequisite before spending on point-based evaluation. Recheck this record
when a new floor, a complete native cost ledger, or a CM proof dependency lands.

## Resolution and history

2026-10-15: first native implementation, 100-reference equality checks, 15 new
floor certificates, 8,000 orbit controls and conditional Lean proofs recorded.
Full end-to-end and replication campaigns remain governed by the study ledger.
