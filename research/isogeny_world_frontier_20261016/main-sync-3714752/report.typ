#set page(paper: "a4", margin: 18mm)
#set text(font: "Libertinus Serif", size: 10.5pt)
#set par(justify: true)
= Current-main merge: local verification complete
16 October 2026 · source integration record

This report records integration and test observations. The separately sealed P192 degree 50341 scientific replay continues; no mathematical admission or speed claim follows.

#image("flow.svg", width: 100%)

*Source.* Merge 3b7b24b1620f92d9d9c108e18f1688041a14ad2a joins published repair a5c1c39d31 with current main 3714752ee6. Its only conflict was the release workflow draft setting. The resolution retains the fully prefixed version and draft=true until binary upload, download and original-manifest checks succeed. Other upstream bytes are preserved; pushed history was not rewritten.

*Fresh root gate.* The required release-library test command passed at the merged source: 4148 passed, 0 failed, 117 ignored, exit 0. Recorded execution 325.17s plus separate compile 5m52s are shared-host operational durations, with no speed comparison. The isolated test binary hash is e6b0366e…ec7399; the complete hash, command, toolchain and raw log are in README.md and verification-summary.json.

*Publication helpers.* Nine local transport-fixture cases and ShellCheck pass after restoring the tracked sparse helper paths. Both helper sources match the tested repair; the final workflow hash is recorded in workflow-final-check.txt. No root-library source changed during shell materialization.

*Retained failures.* An initial bootstrap stopped on preexisting main SVG/PDF whitespace before test execution. Its complete log is compressed beside this report; the manual workflow resolution passes its scoped whitespace check. First shell checks exited 127/2 because the merge reapplied sparse checkout and removed tracked helpers from disk. Their Git blobs remained present. Initial receipts and the successful materialization/retry are both preserved; no source failure was hidden.

*CI wiring.* A dedicated release-assets job runs the existing shell checks on PRs and main; official publication depends on that job. This workflow-only addition leaves the tested library source unchanged.

*Remaining.* Current-head CI and live repaired official publication remain outstanding. v0.1.2372 still has failed checksum admission. The scientific process uses its own immutable executable and input; these library checks neither restart it nor admit its map. The diagram describes gate flow, with no comparative performance data and no changed canonical mathematical graph.
