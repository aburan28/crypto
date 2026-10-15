# PR1853: catalog and publication CI repair

Scope: follow the already delivered research PR through CI/review. No scientific
campaign, target solve, parameter search or performance comparison is rerun.
The inspected head is 763ab2b06c1be26b45a2e349ddbba773435e5e19, with exact base
a80b15bbc5c1ea4f2db2b294a9bcddc31f06fdd7. The original campaign, protocols,
measured code snapshots, records, reports, PDFs and replays remain unchanged.

## Findings and repair

- Native traits and Jacobian-transfer freshness checks fail on the old head.
  Their registry-derived files omit the two already registered conductor-gap
  models. Rebuild with the original budgets262144 and200000 and unchanged
  native implementations. All2099 old per-curve records remain unchanged; only
  icv1-fp19-t241-6012dab8 and icv1-fp19-t241-b992f548 are appended. The current
  catalog has2101 curves. The transfer catalog has5152 formula rows, with the
  two appended curves classified boundary under the existing cost model.
  No formula parameter, solver price or scientific verdict is changed.
- The base's48 site tests pass. The original current-head run reports panel
  summary/index and internal-link failures, then times out at180 seconds in
  the browser freshness comparison. The old Python producer cannot reproduce
  the native-augmented index. This is a change-caused publication defect.
- Add summary metadata and absolute, immutable evidence links to the two
  generated diagnostic sections; regenerate their panel index. Generate the
  progress point table from the same canonical JSON as its data embed,
  retaining the recorded values, units and incomplete-cost labels.
- Dispatch native browser freshness to its locked Rust producer, retaining
  every cross-reference assertion and the historical producer path for
  historical documents. Both PR and publishing workflow jobs install the
  pinned Rust toolchain and fetch the native document's exact baseline.
  No assertion, workflow, time budget or negative control is disabled.

## Validation

- Full release library:4115 pass,0 fail,117 existing ignored in one invocation.
- Both native catalog check commands pass; the generated Jacobian Markdown
  table matches exactly. Preservation comparison and the raw catalog patch
  show no changes to any old per-curve record.
  The Markdown table is byte-identical: it lists formula-advantage/marginal
  rows, and both added curves have boundary status.
- All48 site tests pass on the repaired tree. Two intermediate retries expose
  tracked inputs absent from the sparse checkout: candidate-reference files
  and one round directory. Materialize those exact tracked inputs; retain the
  failed logs and all original assertions.
- All65 status-workflow tests pass, including equality of the PR and
  pre-deployment step lists. The initial local socket-permission failure is
  retained; the retry permits the suite's local mock HTTP server.
- The touched publisher example test, warnings-as-errors clippy, rustfmt,
  native freshness, display/name checks and source-pinned dashboard check pass.
  The workflow parser's four tests pass as well.
- Desktop,390px/320px mobile, dark-theme, fragment-link, no-JavaScript and
  no-external-request browser checks pass in documented inline mode.
  The new desktop/mobile diagnostic screenshots were visually inspected.
- Original PDF hashes still pass. No original study PDF, graph, report, session,
  specification, library implementation or Cargo dependency changed.

There is no new graphable experimental result. The historical research charts
and benchmark copies remain immutable. Current catalog/table and website views
are refreshed; the research completion and claim gates do not change.

## Commands and evidence

Build environment uses /workspace/toolchains/cargo/bin in PATH,
RUSTUP_HOME=/workspace/toolchains/rustup, CARGO_HOME=/workspace/toolchains/cargo,
and CARGO_BUILD_JOBS=2.

    cargo +1.98.0 build --release --bin curve_traits --bin jacobian_transfers
    target/release/curve_traits check
    target/release/jacobian_transfers check
    target/release/curve_traits build --budget 262144
    target/release/jacobian_transfers build --budget 200000
    target/release/jacobian_transfers table --write-doc docs/curves/JACOBIAN_TRANSFERS.md
    cargo +1.98.0 test --release --lib
    cargo +1.98.0 test --release --example ecbench_open_questions_publish
    cargo +1.98.0 clippy --release --example ecbench_open_questions_publish -- -D warnings
    cargo +1.98.0 build --release --locked --manifest-path research/ecbench_open_questions_20261014/native-publish/Cargo.toml
    research/ecbench_open_questions_20261014/native-publish/target/release/ecbench-ten-question-publish
    python3 research/ecbench_leads_20261015/render_catalog.py
    python3 scripts/site/panel_index.py
    python3 scripts/site/test_build.py -v
    python3 scripts/rho_status/test_rho_status.py
    python3 scripts/test_workflows_parse.py
    node tools/render_ic_dashboard.mjs --check

Presentation Python lays out recorded documents or runs existing site/control
tests; it does not analyze the new native sessions. Native freshness delegates
to the same strict Rust check already used in catalog CI. Scientific outcomes
and seals are preserved, as recorded in frozen-preservation.log.

ci-state-763ab2b06.json is the inspected remote state. Six failed/cancelled
GitHub job logs are retained in ci-logs/ as deterministic gzip files with
checksums. catalog.patch.gz and repair.patch.gz retain the exact changes.
Conductor remains unavailable, exit127; both requested scopes are logged.

## Remaining merge gates

The original head's leaf-m10 assertion again fails at the known exact archive
file-set contract; see the existing baseline-CI problem. No archive, manifest,
verifier or archive workflow is edited.

Rust CodeQL and IC autolab jobs are cancelled during analysis/regressions.
Their logs do not establish a change-caused defect or a successful check.
The separate cancelled site job led to the publication repairs above. No
cancelled result is treated as passing; the new head needs its own checks.

The current Cursor review explicitly says Not approved and requires human
review because Bugbot did not complete its analysis (usage limit/size).
A successful approval-agent job is not an approval. There are no assigned
human reviewers or unresolved inline threads at this inspection. Leave the
PR open until all applicable current-head CI and required review permit merge.
