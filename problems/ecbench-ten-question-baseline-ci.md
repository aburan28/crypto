# Unchanged m10 archive file-set assertion blocks the campaign PR

- **Status:** open
- **Kind:** baseline CI/evidence-layout mismatch
- **Affected work:** [PR1853](https://github.com/aburan28/crypto/pull/1853); [campaign validation](../research/ecbench_open_questions_20261014/VALIDATION.md)
- **Next action owner:** archive/workflow owner; campaign author retains the merge gate
- **Last updated:** 2026-10-14

## Problem

The existing leaf-m10 archive verifier fails its exact directory file-set
assertion on both the unchanged baseline and the campaign head.

## Evidence and uncertainty

- Observed: GitHub workflow38641917091, job116119276893 fails at verify_archive.py
  line48 before any new campaign work runs.
- Observed: the same verifier on exact baseline
  f9011e5375155bd9ef7839abe4eea12718e4b401 fails at the same assertion.
  [Baseline log](../research/ecbench_open_questions_20261014/base-leaf-m10-replay.log)
  and exit1 are preserved.
- Observed: the campaign changes neither archive, manifest, verifier nor its
  workflow. The original manifest lists102 files; the directory also contains
  a separately added run_36707051160 subtree outside that list.
- Derived or proved: no campaign source/count difference explains this file-set
  mismatch. No mathematical obstruction or failed new EC result is inferred.
- Missing capability: an archive-owner-approved versioned file-set contract
  that verifies both immutable artifacts without deleting or ignoring evidence.
- What this does not establish: that pending CI/review passes, or that the PR
  may merge despite the red gate.

## What would unblock it

- [ ] Preserve every original archived byte and independently verify the added
  subtree under its own receipt/manifest.
- [ ] Correct/version the directory scope contract and pass the applicable
  archive workflow on the exact current PR head.
- [ ] Complete all other applicable CI/review gates before merge.

## Next discriminating step

Archive owner checks how the second run belongs in the immutable archive
contract. Do not rewrite frozen evidence, delete the subtree, or disable the
check to make CI green. Recheck on a corrected contract/new workflow activity.

## Resolution and history

2026-10-14: baseline failure reproduced; PR remains unmerged. Unresolved.
