Official ReleaseMe v0.1.2372: binary corruption diagnosis and repair
===================================================================

This is a release engineering observation, not an isogeny frontier result or a
speed comparison. The native scientific profile RUN-ISOU-93b2b3 uses separately
sealed source e82 / binary b3cbd6 and continues independently. No mathematical status
is changed by this report.

Observed publication failure
----------------------------

ReleaseMe workflow 38873645089 reports success for all 12 jobs, including the four
native isogeny targets and release. Tag v0.1.2372 resolves to built source
3714752ee6a9686e0fce79cb308423d9fd09cde9. Release was published
2026-10-16T08:54:06Z. Exact run, tag and asset metadata are retained beside this
report. This is job success, not verified binary delivery.

All four downloaded isogeny archives fail the published SHA256SUMS. Each
observed download hash equals GitHub's release asset digest, ruling out a
local-only download discrepancy for these exact bytes. The macOS archive
starts 1f ef bf bd 08 rather than the gzip magic 1f8b08. No failed-checksum
download was extracted or executed. Existing published bytes were not replaced.
verification-summary.json gives all exact expected/observed digests and sizes.

The original macOS CI artifact from the same immutable run has SHA256
31916732d06a215331654e428c901efd4af3de6a42284001697995e21ade258f,
matching the original manifest. Decoding that correct archive as UTF8 and
encoding it back to bytes reproduces the public asset SHA256 exactly:
da8b8e6946a7bdfb067e4daa7d11993be6523bb110f731e6c3343c2df5051b7a.
This is a deterministic metadata-byte reproduction, not cryptographic search.
The pinned ReleaseMe source at cede87fd8a4829ef571dff6bde807ebbfcaeb3e9 reads
upload files with encoding utf8 (release-me-pinned-assets.ts). Primary source:
https://github.com/dev-build-deploy/release-me/blob/cede87fd8a4829ef571dff6bde807ebbfcaeb3e9/src/assets.ts

Verified pre-publication installation control
--------------------------------------------

After checksum and gzip validation, the original CI archive was extracted in
scratch outside the checkout. Its standalone executable SHA256 is
618766b9cd1f1dc9410a46fb4c40d64cbd04f7b6cefb7a847d6663294566140a;
version reports source 3714752..., dirty=false and the expected standards hash.
Direct native `search --curve p192 --ell 13` and a separate `verify --input`
both exit0. Both eigenline maps pass kernel, codomain, exact rational-map and
subgroup checks, with 20 public scalar transports per map. Full raw controls
are retained in original-ci-control.tar.gz with per-file hashes. These checks
validate the original CI payload only; they do not make the corrupt release
assets installable. No performance number or large-degree improvement follows.

Repair and validation
----------------------

ReleaseMe still creates the official release and uses the fully prefixed tag
(the prefix fix already landed on main). It now creates a draft without passing
binary files to its UTF8 uploader. publish_release_assets.sh checks the local
manifest, uploads with gh's byte-preserving file transport, downloads into a
fresh directory, compares the downloaded manifest to the original local one,
and checks every downloaded digest before publishing the draft as latest.
Existing assets and existing verification directories refuse rather than overwrite.
A transfer, checksum, or publish failure fails the job and leaves the draft.

Nine local shell transport-fixture cases pass: exact high-byte gzip/zip payloads,
corrupted download, replaced remote manifest, failed upload, failed download,
failed publication, bad local checksum, existing output and existing assets.
They exercise actual bytes via a mocked local gh transport; they are not a live
GitHub publication test. Shell syntax, ShellCheck and workflow YAML wiring pass.
Initial eight-case logs are preserved before the draft-publication failure case
was added; final logs contain all nine. No Rust/Python/example source changed.
The unchanged root-library gate remains4148 passed, 0 failed, 117 ignored, and
the dedicated-tool source gate 94 passed, 0 failed, 1 exercised helper definition
ignored. Their earlier full receipts remain in this study's root/replay reports.

The replacement workflow has not run at this report checkpoint. Acceptance of
future release delivery requires the current-head CI jobs, successful live
byte-verification step, exact tag binding and a native control on a downloaded
released executable. v0.1.2372 remains recorded as failed artifact admission.
Do not use its corrupted archives as an installation route. A later good release
will have its own version and evidence; this record is not rewritten.

The diagram shows the failed and proposed paths because this is one custody
incident, with no comparable performance series. Canonical mathematical curves,
challenge seasons, scores and IC/rho graphs are unchanged. Editable SVG/Typst and
an inspected one-page PDF accompany this report.
