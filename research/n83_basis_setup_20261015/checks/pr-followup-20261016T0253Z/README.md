# Publication CI follow-up — 16 October2026,02:53UTC check

PR1899 publication head14e56ddb0902161eaa289425c5131a58fefd5243 was ready/open/
mergeable with17 observed workflows:14successes and3failures. Cursor review
COMMENTED at02:04:37UTC: PR too large/incomplete automated analysis, no approval,
human review required; no reviewer assigned/requested. No merge attempted.
This is CI/presentation work only, not another experiment. Scientific completion
watch stays disabled; missing cohort comparison remains unresolved.

## Change-caused repairs

1. [ic-leaderboard run38838216432](https://github.com/aburan28/crypto/actions/runs/38838216432)
   identified five missing input paths. Two are new native basis results. Both
   new paths are added to BOTH push and pull_request lists, preserving the path
   coverage guard. The remaining three paths were already inputs of the parent
   builder and absent from both parent filters; see base-failure proof below.
2. [Dashboard run38838216453](https://github.com/aburan28/crypto/actions/runs/38838216453)
   correctly rejected stale overview metadata after the new timeline/panel.
   Existing renderer now explicitly retains BOTH the historical matched N83 panel
   and the new basis panel. Their published section bytes are unchanged. Marker
   retention keeps the native-value presenter/overview renderer composable.
   Generated page/data refreshed; no evidence pin updated or guard removed.

Native scientific workers/controllers/v1postprocessing, protocols, fixtures,
raw records, old reports, final native JSON, PDFs, hashmanifest and evidence ZIP
remain byte-identical. PUBLICATION_ARTIFACTS_UNCHANGED.json verifies the final
publication artifacts against14e56ddb0. New follow-up evidence is add-only here;
the old archive/receipt is not repackaged or edited.

## Unrelated inherited failures, preserved

[Leaf archive run38838216481](https://github.com/aburan28/crypto/actions/runs/38838216481)
failed verify_archive.py line48: manifest path set differs from all files under
the archive root. The base59e9c695b4df32bcdb0702b83cff86db1a761e4c has the EXACT
same archive Git tree, manifest, verifier and workflow. Its manifest declares102
paths, all present, while the Git tree contains207 files, including105 additional
paths (largely another nested run archive). Full path/object proof is
BASE_FAILURE_EVIDENCE.json. None of this archive's evidence is deleted/resealed
or patched by the basis publication follow-up.

The same base builder already reads the three n83_matched_20261014 statistics/
replay/supplement inputs absent from both base path filters. They remain uncovered
and the GLOBAL path guard still FAILS on them, rather than being disabled or
falsely reported green. Base workflow snapshot and independently checked path
coverage are retained. The child only fixes its two newly introduced paths.
These inherited issues remain CI/merge gates for parent/dependency resolution.
No claim is made that a base GitHub workflow ran: equality of inputs and the
base coverage check establish that the defects precede this child publication.

## Validation and runtime limits

-9 existing path-check unit tests pass. Both new native sources are covered for
 push and pull_request. Global check now correctly reports only3inherited paths.
-Node renderer --check passes all frozen evidence pins/status/historical ledger.
-Composition regression compares both N83 section bytes to14e56ddb0, verifies
 unique IDs and exact idempotence of native-value presenter + overview renderer.
 First Node child-process attempt hit the managed shell's EPERM; the same test
 passes with the normal network-enabled execution permission, with guards intact.
-48 existing site tests pass after regeneration; generated leaderboard/browser
 and panel-index checks pass. No Rust/scientific source changed; original
 post-campaign4089/0 release-library result remains bound in published evidence.
-Local browser check could not validate the page: container Chromium first lacks
 a usable SUID sandbox, and after normal container launch setup the managed
 browser denies file:// navigation (“Your organization doesn't allow you to view
 this site”). Diagnostic screenshot/log retained. No enterprise policy is
 changed or bypassed and no browser test/expectation is weakened. GitHub CI's
 unchanged browser check remains required; local renderer/site passes do NOT
 imply browser pass.

Required human review including Cursor remains unsatisfied. No reviewer request,
merge, experimental rerun, extrapolated batch value or ledger promotion occurs.
New-head CI must be checked at the next activity cadence, not foreground-polled.
