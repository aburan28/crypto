Corrected clean native executable and complete retake
====================================================

Source: e82d80feaf82d3730576ccd9137911b336f9dab8, public PR1976.
The checkout was asserted clean before building; version.json reports that exact
source with git_dirty=false. The executable is
/Volumes/SSD990-2/crypto/isogeny-degree-challenge-builds-20261015/compact-final-e82d80feaf/isogeny,
36087488 bytes, SHA256
b3cbd6ac5d457a0f6d6398d1fcefc408f627e60be52c81c9009406486b61587b.
The production copy is read-only. All full native tests invoke a byte-identical
compiled target; compiled-installed-comparison.exit-code.txt is0. This build
uses no Cargo at scientific runtime.

The clean full suite passed94 tests, zero failed, with one standalone helper
ignored and exercised explicitly by its parent interruption control. Build and
control exits are0. The focused native writer-error control was retaken with
nocapture so its expected error/cleanup stderr is retained; exit0. Source/core
arithmetic is unchanged from the successful precommit94-test suite; the change
here is committed clean build provenance.

raw-controls.tar.gz is the complete raw build/control package,110195 bytes,
SHA256 bf0b366d4ce08e7c23ff1287ed8789fa83b36068a758fb6a75ae321d37f3273e.
It contains every compact/trace/RSS/interruption fixture record, logs, source,
version, executable hashes and a per-file raw-controls.sha256 manifest. Extract
the archive to inspect individual records; no raw control was replaced by a
summary. Synthetic RSS and no-cryptographic-work fixtures remain labelled.

Both objects are sealed under fresh write-once names in the configured verified
owner store. binary-store-push.json and controls-store-push.json bind their exact
hashes/sizes. store-remote-verification.txt confirms cache/remote verification.
store/index.json holds logical names:
- evidence/isogeny-world-frontier-20261016/isogeny-compact-trace-e82d80feaf-macos-arm64
- evidence/isogeny-world-frontier-20261016/compact-trace-clean-controls-e82d80feaf

The native version's standards_sha256 is
22b2ffef99720f2b2c7396a267412d2241a868e038a835c4fb1fd12209b102bf;
this raw registry digest is distinct from a derived challenge-catalogue digest.
The clean source/build/proof is ready for Coordinator input binding and final
protocol freeze. No50341 profile, regenerated map, larger expanded map or M8191
scientific target was executed by this retake. The old50341 null score remains.
