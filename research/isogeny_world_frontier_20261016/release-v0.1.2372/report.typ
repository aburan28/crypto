#set page(paper: "a4", margin: 17mm)
#set text(font: "Libertinus Serif", size: 10pt)
#set par(justify: true)
= Official isogeny release: failed byte admission
#text(size: 10pt)[v0.1.2372 · run38873645089 · 16 October 2026]

This report records a release engineering observation. The separately sealed P192/50341 scientific profile continues; no mathematical frontier or speed claim changes.

#image("flow.svg", width: 100%)

*Observed.* All 12 ReleaseMe jobs report success, and tag v0.1.2372 resolves to built source 3714752ee6a9686e0fce79cb308423d9fd09cde9. All four downloaded isogeny archives fail the published checksum manifest; their hashes match GitHub asset digests. No failed-checksum archive was extracted or executed.

*Exact cause.* The original macOS CI archive matches the manifest (31916732…de258f). UTF8 decoding/re-encoding reproduces the public corrupt asset hash (da8b8e69…051b7a) exactly. The pinned ReleaseMe upload source reads files as UTF8. Full hashes, original source and raw headers are retained beside this report.

*Native control.* After checksum validation, the original CI macOS executable reports source 3714752, dirty=false. Public P192 degree 13 construction and fresh independent replay both exit0; both eigenline maps pass kernel, codomain, exact map and subgroup checks, with 20 scalar transports per map. This validates the pre-publication payload, with no performance comparison.

*Repair.* ReleaseMe creates a draft. A shell helper validates the local manifest, uploads files as bytes, downloads into a fresh directory, binds the downloaded manifest to the original and checks all hashes before publishing the draft as latest. Any failure leaves the draft; existing assets refuse overwrite.

*Local checks.* Nine transport-fixture cases pass: high-byte gzip/zip roundtrip, corrupted download, changed manifest, upload/download/publication failure, wrong local checksum, existing output and existing assets. Shell syntax, ShellCheck and workflow YAML wiring pass. No Rust or Python source changed; earlier source gates remain 4148 root-library and 94 tool tests passed.

*Remaining.* The repaired workflow has not executed. Official delivery still requires current-head CI, live downloaded-byte verification and a native control on the released executable. v0.1.2372 remains failed checksum admission. Canonical curves, seasons, scores and mathematical graphs are unchanged. The diagram represents custody paths; there is no performance data series.

#text(size: 8pt)[Primary source: github.com/dev-build-deploy/release-me/blob/cede87fd8a4829ef571dff6bde807ebbfcaeb3e9/src/assets.ts. Complete evidence and commands: README.md, verification-summary.json and raw receipts.]
