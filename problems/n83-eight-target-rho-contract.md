# Degree-83 eight-target rho workload rejected by the frozen worker

- **Status:** open
- **Kind:** implementation gap and missing matched measurement
- **Affected work:** [bounded study](../research/n83_basis_setup_20261015/README.md), PR1899
- **Next action owner:** unassigned; a separately authorized prospective correction
- **Last updated:** 2026-10-16

## Problem

The predeclared eight-target shared-setup IC/rho comparison cannot be established:
the frozen rho worker accepts only one or32 points and rejected all three
8-point cohort inputs. This is the exact legacy a=1/tail149 degree83 curve,
r8569786107849059, Fold/reference-u64/256lanes/DP6/escape-off/capfactor2000,
CPU3/oneOSthread/12GiB/2100s. The cold singleton comparison is independently valid.

## Evidence and uncertainty

- **Observed:** all three cohort rho processes return2 before any answer;
  stderr says “a singleton or exactly32 independent public points is required”.
  [Frozen source](../examples/n83_basis_rho.rs) lines107–108 and
  [raw records](../research/n83_basis_setup_20261015/records/cohort-campaign/)
  establish the input-contract mismatch. All six cohort IC rank/log audits and
  the final rho A/A pair pass; the stage receipt fails closed (8/11 processes,
  50/74 answers). The frozen complete-stage analyzer correctly refuses admission.
- **Derived or proved:** no batch IC/rho ratio or interval is identifiable from
  these records. A millisecond failed launch is not a rho solve cost.
- **Suspected explanation:** no unknown runtime cause remains for these failures;
  the inherited32-target input guard was not reconciled with the new8-target plan.
- **Missing capability:** measured, independently verified eight-target rho cohorts.
- **Limits:** this is not a mathematical obstruction or evidence of rho difficulty.
  No old32-target timings, cold rho times, projected times or new replacement runs
  fill this gap. C4 training allocations are a separate resource-policy failure.

## What would unblock it

- [ ] Separately version and freeze a corrected8-point workload contract, with
  native positive8-point and negative duplicate/malformed-input tests before timing.
- [ ] In a separately authorized prospective campaign, complete matched fresh
  cohort triplets under the full common envelope and independently audit every
  answer/transcript. Preserve this study and all failed attempts unchanged.

That would establish a new local batch comparison; it would not retroactively
complete these three failed cells or establish independent-host/L2/security claims.

## Next discriminating step

Review the worker/controller workload contract and add the missing8-point native
acceptance regression in a new version. No correction or new solve is executed
as part of this publication. A new campaign would retain the original resource
limits unless prospectively changed by its protocol. If acceptance fails, stop
before timing; if measurements fail, preserve them and withhold a winner.
Recheck on a separately approved correction study.

## Resolution and history

2026-10-16: failure found during post-campaign independent audit. No raw directory,
worker, protocol or audit was rewritten and no solve retried. Unresolved.
