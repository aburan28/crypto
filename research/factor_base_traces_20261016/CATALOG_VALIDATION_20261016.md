# Historical unpublished catalog-validation receipt

Recovered on 2026-10-17. This receipt belongs to the older, **unpublished**
`d78eed09` + `d2225365` local integration snapshot. Another session subsequently
published main merges and catalog maintenance through `0678100c`. This snapshot
is retained as historical validation and is not pushed over that later work or
claimed as validation of its final source.

This is PR #1975 catalog maintenance, not an experiment, native timing study,
new algorithm implementation or end-to-end performance result. Historical
CPython generation was not rerun.

The maintained input head is `d78eed09a48bec50b8fab3cafc9c0203fb7d5b99`;
the additional main parent is `d22253655de460de5fbf4774959c97d8a03de682`.
The original Copilot request was still unanswered after nearly ten hours.
The [continued fallback was documented before edits][takeover]. A read-only
merge-tree check found only additive store-index entries in conflict.
The existing native index helper preserved all **23** entries from both
parents, refusing differing values for any shared name. No remote object
was uploaded, replaced or removed.

## Change-caused CI failures

The registry has 2,160 models after the 48 historical registrations, while
the old trait and Jacobian catalogs each still had 2,112 records. Their
staleness checks failed, although the relevant unit tests passed. The
leaderboard roster and registry digest were also stale. The failed CI logs
and local reproductions are retained separately from successful checks.

- Native `curve_traits build --budget 262144` refreshes the derived traits.
- Native `jacobian_transfers build --budget 200000` refreshes its catalog
  with the existing cost model and budget. Every previous catalog record
  compares unchanged; the extra 48 records are metadata from existing native
  formulas and invariants, not new measured speedups or ranking claims.
- Existing metadata-only `scripts/update_curve_standards.py` refreshes the
  leaderboard roster and rendered displays. Every frozen measurement field
  is unchanged. Its check and the full leaderboard render check pass.

The native catalog builders' implementation did not change in the new main
parent. No new Python cryptographic implementation or measurement producer
was introduced. Future measured comparisons still require symmetric native
operation pricing bands; a verified engineering improvement counts without
requiring a win over rho.

## Preservation and validation

All 24 non-store creation payload blobs retain their original hashes. Both
trace bundles, archived column order, recorded provenance and signed lookup
caches are unchanged. The original report, diagrams, PDF and failed authoring
logs remain unchanged. The trace-flow diagram has no changed semantics and
the existing research plots remain frozen. The affected leaderboard displays
are regenerated together; no new measurements require benchmark plots.

The #1963 card remains draft, unrun, zero-funded and unleased. No protocol,
checkpoint, runnable command or campaign was promoted.

The first library compile encountered a tracked `hwcost` module outside the
sparse checkout. That compile was interrupted (status 130) after its E0583
diagnostic; the initial log is retained. The existing module was materialized
unchanged before retrying. No test assertions, limits or ignore annotations
were changed.

The recovered pipeline completed all of its sequential commands. The release
library command passed **4,168 tests, zero failed, 117 existing ignores** in a
single invocation with networking enabled. Production native build, all-targets
compilation, strict Clippy, four trace tests and verification of all 384
source-backed trace/cache records passed. Both full native catalog staleness
checks passed with the original budgets; the metadata display/render checks,
naming check (zero problems) and exact preservation audit passed. These are
correctness/integration results on this historical source, not new native
performance measurements or final-head CI results.

Exact commands and compressed original logs are in
`validation/catalog-20261016/`. `log-manifest.json` binds each retained log to
its original byte length and SHA-256. It includes the failed CI jobs, failed
staleness reproductions, initial sparse-source compiler diagnostic and every
successful retry. Existing authoring and `merge-20261016` receipts/logs remain
byte-identical.

Conductor is unavailable; its workflow file remains unavailable (404), not a
passing check. Rust CodeQL on the input head was cancelled at its existing
25-minute PR limit during query evaluation. That is an execution blocker,
not a pass, and no timeout, coverage or review gate was weakened. All applicable
CI and required reviews must pass on the resulting exact head before merge.

The skipped-job conditions were inspected: release publishing requires a
push or manual dispatch on `main`; performance-history recording excludes
PR/merge-group events; long Rust CodeQL is manual-dispatch-only. Their ordinary
PR validation jobs remain applicable. A cancelled normal Rust CodeQL analysis
is never substituted by one of these skipped jobs.

[takeover]: https://github.com/aburan28/crypto/pull/1975#issuecomment-6220383469
