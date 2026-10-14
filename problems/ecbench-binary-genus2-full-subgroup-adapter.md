# ecbench adapter for full binary genus-two image solves

- **Status:** open
- **Kind:** implementation and measurement integration gap
- **Affected work:** [full-subgroup comparison](../research/binary_genus2_speed_20261014/PROTOCOL.md)
- **Next action owner:** current campaign implementation
- **Last updated:** 2026-10-14

## Problem

The existing ecbench runner cannot yet seal and independently replay the full
elliptic-versus-genus-two comparison at both pinned GF(2^83) subgroups: its
wide Koblitz constructor assumes cofactor 4 or 2, whereas a=1 needs cofactor
1128547018, and its group adapter does not represent Mumford divisors.

## Evidence and uncertainty

- **Observed:** `src/cryptanalysis/koblitz_wide.rs::WideKoblitz::new` sets
  cofactor 4/2 and rejects a composite remaining order. `ecbench/methods.rs`
  dispatches the existing wide strong-rho path; neither its workload nor its
  counted-group adapter exposes the constructed Jacobian image.
- **Available:** `tools/binary-genus2-cover/src/campaign.rs` uses the existing
  signed-Frobenius implementation with the pinned prime orders and validates
  its eigenvalue on the pinned generator, so the full-size arithmetic itself
  is usable. Native records and restart replay are separate from ecbench
  sealed sessions. They must not be represented as ecbench audit receipts.
- **Missing:** registered width-generic image workload/method adapters,
  complete common-unit cost calibration, independent-host replay and an
  admissible isolation certificate. Native E additions and checked Cantor
  additions are different work units.
- **Not established:** any speedup, optimum native Jacobian arithmetic cost,
  full solve duration for an unfinished search, or obstruction to the cover.

## What would unblock it

- [ ] Add explicit, validated prime-subgroup selection to the wide workload
  without changing prior registered workloads' replay behavior.
- [ ] Register the image adapter and deterministic resumable method as new
  method IDs, retaining setup, mapping, lookup, native arithmetic and external
  verification work; audit sealed public-target sessions with replay.
- [ ] Price both native work vectors in a measured common unit and repeat
  matched controls under admissible isolation and on an independent host.

The first two establish harness integration. The final condition is still
needed for a qualified performance claim; none predicts a faster algorithm.

## Next discriminating step

Run the public-job native correctness/restart checks, calibrate the checked
arithmetic paths, and execute full-order public targets using the already
available arithmetic. Preserve checkpoints and censored outcomes. In parallel
with later integration work, keep common priced work and speedup null.
The full subgroup orders and m=83 remain fixed; there is no replacement by
small or interval targets. A successful replay unblocks restart reliability,
not the remaining ecbench integration or independent-host gate.

## Resolution and history

Unresolved. Recording this gap does not finish the requested campaign.
