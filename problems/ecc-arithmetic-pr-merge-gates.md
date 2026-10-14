# ECC arithmetic PR merge gates

- **Status:** `investigating`
- **Kind:** validation failures, stale generated artifacts, curve identity defect, and an unverified CI repair
- **Affected work:** [ECC arithmetic PR #1798](https://github.com/aburan28/crypto/pull/1798), [platform measurement problem](ecc-arithmetic-platform-measurement.md), [CI repair PR #1812](https://github.com/aburan28/crypto/pull/1812)
- **Next action owner:** authoring session for #1798, with independent review for curve identity and frozen evidence corrections
- **Last updated:** 2026-10-14

## Problem

PR #1798 cannot meet the repository's merge gate while its current head lacks a completed, passing check suite. The completed suite for head `9e747fafd14c7eb8725d53bbd28c7e8c500f7402` had twelve failed jobs; the later head `333fcbb68` contains a focused rustfmt correction and requires its own complete check suite. The failures concern several different components, so one successful arithmetic test or one repaired workflow cannot establish merge readiness.

The requested ECC work remains scoped by [the platform measurement problem](ecc-arithmetic-platform-measurement.md). This record concerns delivery and validation of the implemented slice; passing CI would not establish a measured whole-curve speedup or complete the wider backend portfolio.

## Evidence and uncertainty

- **Observed:** The previous head was approved and mergeable in GitHub, but its check rollup contained failures in `rustfmt`, `clippy`, `ecbench`, curve cover certificates, curve traits, Jacobian transfers, perf history, IC integration, and four isogeny release targets ([PR #1798 checks](https://github.com/aburan28/crypto/pull/1798/checks), provenance `internal`). The changed-file rustfmt log identified `src/binary_ecc/f2m.rs:1616`; moving one inline polynomial comment made `rustfmt --edition 2021 --check src/binary_ecc/f2m.rs` pass locally, and commit `333fcbb68` was pushed.
- **Observed:** The perf-history parser read `perf-history/docs/perf/history` on PR runs although the workflow skipped checking out that branch ([failed job](https://github.com/aburan28/crypto/actions/runs/38501436964/job/115649076254), provenance `internal`). The branch contains that directory. PR #1812 checks out the branch for those runs; its CI result is pending.
- **Observed:** The four isogeny matrix jobs failed during Cargo's package fingerprint walk with `Failed to update the excludes stack` ([release run](https://github.com/aburan28/crypto/actions/runs/38501437080), provenance `internal`). PR #1812 moves target output outside the package and disables the nested target cache; whether this clears the failure is unknown until that matrix runs.
- **Observed:** The Jacobian transfer check said its committed JSON is stale against `registry.json` ([failed job](https://github.com/aburan28/crypto/actions/runs/38501437123/job/115649062572), provenance `internal`). A prebuilt local generator emitted obsolete schema `/1` while current source declares `/2`; its output was rejected and not committed.
- **Observed:** Curve traits rejected `icv1-fp21-tm700-56b6783a` because its stored order does not equal an exhaustive count ([failed job](https://github.com/aburan28/crypto/actions/runs/38501437222/job/115649063324), provenance `internal`). The cited source declares `y^2 = x^3 + a x^2 + b x`; the retained registry row instead describes `y^2=x^3+a*x+b`. A native normalization cache and correct normalized rows also exist. The registry's append-only carry-forward currently retains the older, mismodeled row, so a one-field order edit would not repair the identity.
- **Observed:** Curve cover reports a stale linkage graph ([failed job](https://github.com/aburan28/crypto/actions/runs/38501437188/job/115649063256), provenance `internal`); ecbench's negative-control expectation differs from its current validation errors ([failed job](https://github.com/aburan28/crypto/actions/runs/38501437186/job/115651162990), provenance `internal`); IC integration reports three failures and twelve errors in frozen evidence tests ([failed job](https://github.com/aburan28/crypto/actions/runs/38501437056/job/115649063132), provenance `internal`). Clippy reports 27 errors in untouched cryptanalysis files ([failed job](https://github.com/aburan28/crypto/actions/runs/38501437170/job/115651248848), provenance `internal`).
- **Derived or proved:** The equation forms in the source and retained registry row differ. This is an identity mismatch for that row, not evidence about the security or speed of either curve.
- **Suspected explanation:** The isogeny fingerprint error may involve cached target files inside the package, but this has not been isolated. Other failures have distinct observed symptoms and are not attributed to the ECC arithmetic patch merely because they ran on its PR.
- **Missing capability:** A validated correction for the carried-forward curve identities and dependent catalogs, reviewed repairs for frozen evidence and baseline lint failures, and complete passing checks on the final PR head.
- **What this evidence does not establish:** It does not show the arithmetic kernels are incorrect, that the CI repairs work, or that the remaining failed checks can safely be skipped. Earlier-head results are not final-head validation.

## What would unblock it

- [ ] Identify and correct the explicit-model import/carry-forward error without assigning the old identity to the corrected model; regenerate and independently check the registry, traits, transfer report, cover graph, and every dependent pin affected by the change.
- [ ] Make the ecbench and IC frozen-evidence checks pass with their expected fail-closed semantics preserved, and resolve the Clippy errors in their owning source paths without disabling checks.
- [ ] Validate the perf-history and four isogeny matrix repairs on their own PR head, then bring the accepted changes into #1798 by merging `main` into its branch.
- [ ] On the final #1798 head, obtain completed successful required and applicable checks, no unresolved blocking review thread or conflict, and approval where required; merge with the permitted merge-commit method and confirm the merge commit.

These conditions are necessary for this PR's delivery. They would establish merge readiness, not a qualified arithmetic performance advantage or completion of the wider curve and platform request.

## Next discriminating step

- **Action or proof obligation:** Review #1812's perf-history and isogeny job results, then run a separate identity correction from the pinned native normalization cache against the retained registry rows before regenerating dependent catalogs.
- **Inputs and controls:** Compare `out_seed7.json`, `docs/curves/sources/prime-models.json`, the old and normalized ICV1 preimages, and exhaustive point counts; preserve the old registry bytes and inspect all affected references. Keep the ECC arithmetic diff as the unchanged control.
- **Execution limits and stop condition:** Use the existing finite catalogs and standard CI jobs. Stop promotion on an identity conflict, changed frozen evidence without replay, or a failed current-head job; retain each failed receipt.
- **Outcomes and next actions:** If the normalized identities and catalogs replay, review and merge their scoped correction, then rerun #1798 on a new head. If they conflict, record the exact models and unresolved mapping rather than choosing an order. If CI infrastructure fails, repair that path and rerun its gate without treating it as mathematical evidence.
- **Recheck trigger:** A merged prerequisite PR, an independently reviewed identity correction, or new completed checks on #1798's final head.

## Resolution and history

2026-10-14: Opened from the completed previous-head CI receipts and the focused rustfmt correction. Unresolved. Recording this problem does not satisfy any merge or performance criterion.
