# PR 1850 main-merge validation — 2026-10-15

This is publication validation of the existing conductor-gap instruction
census. It adds no new solve, measurement, target, claim or hardness result.
The original evidence and its earlier failed library run remain unchanged in
`../ecbench_conductor_gap_ir_20261014/`.

## Revision and reconciliation

The existing Copilot request, comment 6177073533, ended with GitHub timeline
event `copilot_work_finished_failure` at 2026-10-14T19:55:29Z. The authoring
session documented the AGENTS fallback in PR comment 6181895530 and merged
main without rebasing or force-pushing.

Validated merge: `efb72d1954b2c687bdf79fa83a34aaffdc265ecb`, tree
`d6380d909c5628595b9cc16eec11351b79b5eb55`; parents:

- PR head: `cab08ef0a3644b2ca36261b90b4908dc140db86e`.
- Main snapshot: `34b2e238ab6cbbc9b94574a9aedee907834949f2`.

The benchmark index retains the exact union of both parents' rows: 19
entries. Conflicted leaderboard/browser outputs were regenerated with
`scripts/build_ic_leaderboard.py` and `scripts/build_lab_browser.py`;
`scripts/site/panel_index.py` regenerated the panel index. Browser totals
retain main's 102 sessions, 183 methods, 363 factor bases and 12,442 yields.
The four live instruction-census timeline points, embedded data and table
are dated 2026-10-15. Frozen report and study dates were not changed.

All 257 tracked files across the two original study trees are unchanged.
The paired `*-before.tree` and `*-after.tree` files are identical Git
object inventories. The measurement module and native example also retain
their pre-merge bytes; `unchanged-source.sha256` records their digests.
Fresh archive extraction and native Rust analysis reproduced both original
analysis files byte-for-byte; see `native-replay.sha256`.

## Completed validation

The wrapper began at 2026-10-15T01:03:25Z and completed at
2026-10-15T02:07:45Z. These are correctness-test execution times, not new
algorithm performance observations.

| Check | Observed result | Evidence |
|---|---|---|
| Full release library, one test thread | 4,102 passed, 0 failed, 117 ignored | `library-tests.log` |
| Native driver example | 2 passed, 0 failed | `driver-tests.log` |
| Process-exclusive strict-session test | 1 passed, 0 failed, 0 ignored | `strict-session.log` |
| Release library `cargo check` | exit 0 | `cargo-check.log` |
| Site | 48 passed | `site-tests.log` |
| Rho status | 65 passed | `rho-status-tests.log` |
| Curve names, full and diff lint | 0 problems | `curve-names-final.log` |
| Browser freshness | current | `browser-freshness.log` |
| Leaderboard and panel-index freshness | both exit 0; 303 panels current | official tools, checked before publication |
| Scoped Rust formatting and workspace formatting check | both exit 0 | rustfmt / cargo fmt checks before publication |

The full library run includes unchanged passing executions of both JV-cover
tests that failed in the earlier concurrent run. This does not erase or
reclassify that earlier failure. The strict-session test was also executed
separately with `--exact --include-ignored`, and actually ran one test.

`validate.sh` and the `.command` files retain the complete commands. Build
settings were release opt-level 3, default features, three compilation jobs,
256 codegen units and one library test thread. This build supplied no
experimental instruction counts. The strict-session command used the same
settings and native test binary through `tools/isolated_bench.py busy`.

The initial curve-name attempt stopped because a tracked skill file was
absent from the sparse checkout. Its failed log and status are retained;
after restoring the tracked prose/pin inputs, the full and diff lint passed.
No frozen label or test was changed to satisfy the check.

The auxiliary autofix run 38640843389, job 116115931909, failed before
checkout at `Require ANTHROPIC_API_KEY` because its repository service secret
was not configured. This external service setup failure is reported on the
PR, not patched around. Conductor remains unavailable in the workspace;
the required scope check was attempted and exited 127.

GitHub CI and review requirements still apply to the published PR head.
These local results alone do not certify the final GitHub merge result.
