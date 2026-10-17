# Photographed-volcano experiment checkpoints are unavailable

- **Status:** open
- **Kind:** access/tool limit; unknown cause of missing files
- **Affected work:** [photographed-curve campaign recovery](../research/photo_volcano_recovery_20261017/README.md)
- **Next action owner:** unassigned; archive custodian or restored runtime
- **Last updated:** 2026-10-17

## Problem

The existing matched ECDLP campaign cannot resume or be independently audited
because its graph, original random fixtures, native adapters, raw measurements
and checkpoints are absent from the currently accessible storage.

The selected source is `icv1-fp57-t354478198-afcea726` over
`F_72057980207847001`, with subgroup order `18014494963342201`.
The original requirement is ten generators times ten paired random scalar
values on every curve of the selected 2-isogeny component, measured with rho,
BSGS and full IC recovery. Prior metadata reports 46 curves and 90 directed
maps across five levels; that topology certificate is also unavailable here.
Regenerating fixtures would start a different experiment and would not recover
or preserve the earlier campaign.

## Evidence and uncertainty

- **Observed:** the original workspace directory, reference snapshots,
  extracted temporary archive copy and matching archives were absent from the
  checked workspace, temporary and shared-download locations. A successful
  Docker state read returned no containers or volumes. No new solver was run.
  [Validation scope](../research/photo_volcano_recovery_20261017/VALIDATION.md)
  identifies the checks and limitations.
- **Derived or proved:** no mathematical obstruction established. Public source
  parameters remain retrievable from the pinned shortlist; random fixtures,
  generated maps and prior measurements do not follow from that source alone.
- **Suspected explanation:** unknown. Missing files do not establish when or
  why they became unavailable or that no copy exists elsewhere.
- **Missing capability:** access to the two original full archives or an
  equivalent exact-byte backup and its verification receipts.
- **What this evidence does not establish:** the prior numerical results have
  not been re-audited; the full campaign is not complete; IC recovery is not
  mathematically impossible. A prior zero-yield diagnostic is not full IC.

## What would unblock it

- [ ] Obtain both archives named in the [manifest](../research/photo_volcano_recovery_20261017/recovery_manifest.json)
  and match their expected byte counts and SHA-256 digests.
- [ ] Verify both extracted `SHA256SUMS` manifests and exact fixture, curve,
  solver-source, seed and checkpoint identities without altering old records.
- [ ] Make the recovered originals durable using the authorized immutable
  store, commit a real locator, and verify readback against the same digests.
- [ ] Validate the saved resume path on the current host, retaining source,
  binary and replay receipts before starting additional measured cases.

These conditions resolve checkpoint access and continuity only. They do not
complete the remaining rho/BSGS cases, full-cost timing protocol, full IC
recovery or its viable algorithm configuration, independent audits, reports,
heatmaps, or the user's full workload.

## Next discriminating step

- **Action:** restore the two original archives and perform the documented
  archive and extracted-manifest verification before resuming any solver.
- **Inputs and controls:** exact expected sizes/digests in the recovery
  manifest; original fixture bytes and preserved failed attempts are controls.
- **Execution limits and stop condition:** this step reads and verifies only
  the supplied archives. It launches no solver, paid resource, regeneration or
  background campaign. Stop on a digest, manifest or checkpoint mismatch.
- **Outcomes and next actions:** exact matches permit resume-path validation;
  a mismatch remains failed evidence and must be investigated without replacing
  expected hashes; unavailable archives leave the same access blocker open.
- **Recheck trigger:** either original archive or an equivalent verifiable
  backup becomes accessible.

## Resolution and history

2026-10-17: recorded the observed storage gap and prior manifest metadata.
No original raw evidence, solver implementation or checkpoint was recovered.
Unresolved. Merging this record does not resolve the experiment or promote any
reported historical measurement.

2026-10-17 restart update: a fresh replacement campaign was authorized and its
[new frozen inputs and native sources](../research/photo_volcano_restart_20261017/README.md)
are preserved separately. Missing old checkpoints no longer block that new
execution. This record remains open for historical evidence recovery; it is
not a claim that the old bytes were restored or the full workload completed.
