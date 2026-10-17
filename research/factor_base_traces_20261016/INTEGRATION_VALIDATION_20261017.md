# Integration maintenance: 2026-10-17

This is archival/catalog integration maintenance for PR #1975, not a research
campaign or performance result. No historical CPython generator was rerun,
no new point logarithms were obtained and no end-to-end gain is claimed.

## Inputs and resolution

The feature input is `0678100c4e0660f69ef707724d07b22ed3a4b58e`, including
the other session's two main merges and eight generated-catalog repairs.
The pinned main input is `db69b9bd58976dfdadab7623b9c75165b23ff21d`.
Main was merged into that feature head, without rebase or force push.

The original Copilot request remained unanswered after more than 21 hours,
with no visible active resolver. The continued fallback was documented
[before edits][takeover]. A read-only merge-tree check and the actual merge
both found conflicts only in three generated views: `docs/browser/data.json`,
`docs/ic-leaderboard.html` and `docs/ic/LEADERBOARD.md`.

The leaderboard conflict combines the historical-model registry digest with
main's new citations. The existing metadata-only display updater renders the
merged canonical JSON. The existing browser builder regenerates its index
from the tracked frozen sources. No frozen row was selected by a merge side.

## Preservation

The read-only `validation/merge-20261017/check_preservation.py` audit verifies:

- All **55** non-store payload blobs from maintained head `d78eed09` remain
  byte-identical, including all **24** non-store creation payloads.
- All **27** historical evaluation blobs present at creation remain unchanged.
- Both parents' complete model, alias and store records survive unchanged:
  **2,160 models, 5,289 aliases and 23 store entries**.
- Every non-roster leaderboard field equals pinned main's canonical input.
- Both parents' browser identities, exact sessions/rounds and yields survive.
  The combined view has 2,160 curves, 199 methods, 373 factor bases,
  271 candidates, 30 rounds, 111 sessions and 13,402 yields.
- The #1963 card is unchanged: draft, zero-funded, no frozen protocol,
  checkpoint or executable command. No lease or campaign was introduced.

No remote store object was uploaded, replaced, deleted or newly verified.
The original checksum-verified owner-store bundle remains
`factor-base-traces/historical-20261015/registered-20261016.zip`, 2,462,662
bytes, SHA-256
`923579b932103882f85b4475118a5edc7404fa4138823c1bdd59a633ca6f9c9c`.
The maintenance does not renew that historical live-store verification.

The original explanatory diagram and three-page PDF retain exact blobs.
The live browser and all leaderboard renditions are regenerated together.
Other canonical dashboards match pinned main; no new measured result requires
a benchmark point, graph or date change. Future comparisons still require
symmetric native-operation pricing bands and matched references; verified
engineering improvements count without requiring a win over rho.

## Recovered earlier validation

The interrupted older `d78eed09` + `d2225365` snapshot completed **4,168
library passes, zero failures and 117 existing ignores**, plus its native
build, all-targets check, strict Clippy, four trace tests and all 384 record
verifications. Its failed attempts and successful follow-ups are retained
unchanged in `CATALOG_VALIDATION_20261016.md` and
`validation/catalog-20261016/`. The complete historical source/log snapshot is
also retained at [archive commit `a1fdb3de`][archive], on
`codex/archive-pr1975-catalog-validation-20261016`. Its published Git tree
matches the local archive. It was not pushed over the later shared-branch
commits, and its tests do not validate this integration's final sources.

**Provenance correction:** the 4,168-pass count and successful command logs
are observed results. They do not establish that the earlier unit binary
embedded the exact catalog bytes later archived at `a1fdb3de`: catalog
generation overlapped its compilation, and the shared test-binary path was
subsequently reused. The archive binds the recovered working tree and logs;
its exact compiled-fixture equivalence is unverified. Earlier shorthand
describing it as a fully validated source snapshot was too strong. The final
run below pins both the source commit and final fixture before compilation.

## Validation and limits

Validation on source snapshot `552beb731455b9b349732f5b55246bd51c641fde` completed successfully:

- `cargo test --release --lib -- --test-threads=1`: **4,236 passed, zero failed, 117 existing ignores**, in one command with networking enabled.
- Production native build and `cargo check --all-targets` passed.
- Strict release Clippy for `factor_base_trace`, all four native trace tests, and verification of all 384 source-backed traces/caches passed.
- Native trait and Jacobian freshness checks passed. Cover verification passed for all 2,160 models after restoring unchanged off-cone report inputs; the initial missing-input failure is retained.

The final unit binary embeds the exact final catalog bytes (catalog SHA-256 `45e928f3290f1e3f76eb887e233658492d91fff29c6212c658b68f2b0fe198c2`). Its SHA-256 is `1bacf134f82ae128d48ca2006eb3fae20360460147b6e9c19f67c1e319344b2e`. The frozen manifest covers 1,374 tracked native sources and embedded/config inputs. Initial and final production binary hashes differ; no bit-identical-build claim is made.

Completed local checks also include
the preservation audit, naming (zero problems), both leaderboard/browser
freshness checks, five identity tests, three normalization tests, three
archived-label tests and nine workflow-path tests. The workflow filters cover
all 216 files the leaderboard builder reads.

The first full native trait check identified a stale generated catalog after
the merge. Full native regeneration with the existing 262,144-iteration
factoring budget produced 2,336 records (registered models plus the existing
isogeny-source models). Exactly one record description changes: the
`endomorphism.levels` method for `icv1-fp24-t925-8af643e2` says “the crater”
instead of “above the floor”, matching the pinned main implementation. Its
numerical values do not change, and every pinned-main record stays unchanged.
This is generated wording consistency, not a newly discovered transfer or
performance result. The file is also a compiled fixture in the library's query tests.
The library suite was started before that regeneration was complete; its
preliminary run was intentionally stopped after 1,836 passes, returning Cargo
status 101 because the unit process received SIGTERM. Its complete partial
log is retained, and it is **not a passing full-suite result**. The final
catalog fixture and merged native code are pinned in local source commit
`552beb731455b9b349732f5b55246bd51c641fde`, with feature head `0678100c`
and pinned main as its two parents. The final full library run used those inputs.
No assertion, test limit or ignore annotation was changed.

The initial new-worktree setup and incomplete sparse input attempts are
retained as operational failures. The first browser output omitted a candidate
because the materialization helper treated `**/` as requiring a directory;
the parent-preservation audit caught it. The helper now uses the repository's
recursive glob matcher, all tracked inputs are restored read-only, and both
parents' identities are retained. The initial compile-input materializer
also encountered a newly added off-cone source; restoring the maintained
source directories allowed the unchanged repository helper to materialize
all 95 compile-time inputs. No test assertion, limit, ignore annotation or
tracked input was changed.

Conductor remains unavailable (command not found). The pinned workflow file
returns 404. Neither is reported as passing. On input head `0678100c`, normal
Rust CodeQL was cancelled at its 25-minute query limit after 36/37 evaluations;
its neutral summary is not a completed scan. Current-head CI, merge state and
required review policy must all pass before an expected-head-guarded merge.
No check, coverage, resource limit or review gate was weakened.

## Concurrent branch update

Before publication, the shared feature branch advanced independently to
`1a83e379ade9c2cd4c52e426cf4f4777d863e6c8`, merging main
`b59187e5549adb02a43c977f4bea9979a46bf3a3`. GitHub reported it mergeable
with fresh CI queued and a current-head Cursor approval. This completed
local snapshot is archived separately, without overwriting that newer head.
The 4,236-pass result validates source `552beb73`, not the later shared
branch sources. Current-head GitHub checks remain required.

Raw logs, failed/partial attempts, commands and manifests are retained in
`validation/merge-20261017/runs.tar.xz` with a per-member checksum manifest
and deterministic repacking instructions.

[takeover]: https://github.com/aburan28/crypto/pull/1975#issuecomment-6229839449
[archive]: https://github.com/aburan28/crypto/tree/a1fdb3de8f385bb0c5df55a37f1fb9e894b5c474/research/factor_base_traces_20261016
