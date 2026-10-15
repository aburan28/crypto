# Upper conductor-gap floor construction and complete cost receipts

- **Status:** investigating
- **Kind:** resource limit; implementation gap; proof gap (separate issues below)
- **Affected work:** [protocol](../research/conductor_gap_cm_orbits_20261015/PROTOCOL.md), [evidence](../research/conductor_gap_cm_orbits_20261015/REPORT.md)
- **Next action owner:** conductor-gap research workstream
- **Last updated:** 2026-10-15

## Problem

One 96-bit rational-family floor is now independently certified. Other upper
models, the older non-rational-kernel candidate, complete common pricing and
the Lean CM classification bridges remain open.

The rational-torsion family has explicit independently checked floors through
96 field bits. Its other upper inputs have rigorous prime and kernel certificates,
but the first large constructor was stopped after the 91-bit result when shared
memory reached 15,982,964,736 of 17,179,869,184 bytes. This is a resource limit,
not evidence that the floor does not exist.

The older near-limit candidate with prime conductor 129149576694539 is a
separate instance. Its Frobenius scalar does not provide the rational torsion
point required by this square-root Vélu implementation. A new rational-torsion
family endpoint does not construct that exact old candidate's floor.

## Evidence and uncertainty

Observed: all 100 old direct codomains match the fast moment constructor;
16 additional floors have independent PARI order/prime/stratum receipts. See
`results/pair-certificates.jsonl`, `results/floor-certificates-v2.jsonl`,
`results/larger-construction-stop.json` in the study.
The separate `upper-96-floor-certificate.jsonl` certifies prime conductor ratio
125000000302579 over p=78125000378223625457769953627; the torsion certificate
provides independent point-order and Weil-pairing witnesses.

The native relation backend reports field/group work and queries, while modular
exponentiation, table/hash work, allocation and elimination still lack a complete
common price. A zero charged-GAE lower bound is not zero work. No eligible
IC/rho speedup is asserted. The five online phases include independent replay.

Lean proves conditional map-iteration equality, a norm-numerator bound and
scalar-orbit recurrence equality.
The CM norm identification, endomorphism-order lattice, actual rational-map
homomorphism and order classification are unformalized bridges, listed in
`lean/README.md`. No impossibility or general ECDLP hardness gap is established.

## What would unblock it

- [x] Complete a 95/96-bit rational-family codomain under recorded memory/time
  limits and independently certify its exact model, prime group order and floor
  conductor. Preserve failed or censored attempts.
- [ ] Construct and certify the floor for the original non-rational-kernel
  96-bit candidate using a method whose kernel prerequisites actually hold.
- [ ] Count and price every native component, then repeat complete same-point
  IC/rho windows with matching resource envelopes and earned isolation.
- [ ] Formalize the stated CM bridges in a pinned proof environment and attach
  checked proof receipts; arithmetic tests alone do not discharge them.

## Next discriminating step

The bounded-batch constructor matched all 100 reference codomains. Its first
95/96-bit attempts hit the recorded 1,800-second or 6-GiB limits. A separate
11-GiB, 3,600-second retry reached 5,242,880/5,590,170 giant points and then
timed out. A separate precomputed-tree retry of the exact input uses a 14-GiB,
3,600-second cap after solver tables exited. Concurrent test-check compilation
then caused a shared OOM kill (exit 137, 939.36 wall seconds). A new algebraically
equivalent bounded-batch variant has coefficient/reference checks and a launch
guard against local compiler/solver overlap. Its separate retry completed
under an 11-GiB, 5,400-second cap after the combined-tree checks passed.
It used 2,697.60 wall seconds and a maximum 9,362,772 KiB RSS, with no new
OOM kill. Independent PARI prime/order/conductor and torsion-module checks pass.
See the versioned receipts and `native/run_96_v4_solo.sh`.
Fast point transport and dual-map certification on this pair are the next
construction step. Other upper models retain their individual open status.

For the older candidate, test the proposed construction's rational-kernel
prerequisite before spending on point-based evaluation. Recheck this record
when a new floor, a complete native cost ledger, or a CM proof dependency lands.

## Resolution and history

2026-10-15: first native implementation, 100-reference equality checks, 15 new
floor certificates, 8,000 orbit controls and conditional Lean proofs recorded.
The subsequent full campaign verified all 500 primary and 265 replication
attempts. All 47 session audits passed. Independent replay checked 473 IC
records and 797 target recoveries, including controls. Orbit controls now total
44,000 trials, including a pinned local repeat and 20 fresh-process batches.
The primary floor/crater CPU ratio is 0.980 [0.866, 1.111]; endpoint execution
waves and shared-host activity confound timing. The new 96-bit floor and 137th
torsion certificate complete this particular upper construction. Complete
common pricing, other upper floors and the stated CM proof bridges remain open
independently. All completed large-gap ECDLP timing comparisons remain 53-bit.
