Current-main merge and required verification
============================================

This is a source integration and test record, not a scientific result or speed
comparison. The frozen P192/50341 process continues independently with its
sealed e82 executable and unchanged input. It was not restarted for this merge.

Merge 3b7b24b1620f92d9d9c108e18f1688041a14ad2a has parents
a5c1c39d3122be85c4fad62eb980ba25dfe428a0 (published release repair) and
3714752ee6a9686e0fce79cb308423d9fd09cde9 (current main). The sole conflict was
.github/workflows/release.yml: retain draft=true, the fully prefixed version,
and the byte-preserving upload/download verification before publication. Other
upstream main bytes are preserved. No pushed history was rewritten.

The fresh required command in root-release-command.txt passed at this merge:
4148 passed, 0 failed, 117 ignored; recorded test execution 325.17s on the shared
M4 Pro host, plus separate compilation 5m52s. These are operational gate durations,
not paired performance observations. Full stdout/stderr and exit 0 are retained.
Compiled test executable SHA256
e6b0366e26794fe703741ba309cd3c7cb168150c827fa2a23ba068038aec7399,
rustc 1.93.1 (01f6ddf75 2026-02-11), Cargo 1.93.1. The cached target was copied
into an isolated directory before the build; active scientific executable bytes
were not overwritten or reused as a Cargo output.

Nine release transport-fixture cases and ShellCheck pass after materialization.
The two helper source hashes equal the tested repair-source.sha256 in the sibling
release-v0.1.2372 report. No extra root rerun was needed to materialize shell and
report files, which are outside the root-library source dependency closure.

Failed operational attempts are retained:

- Initial bootstrap stopped before commit/test because a full staged whitespace
  check encountered existing upstream SVG/PDF whitespace. Raw output is preserved
  compressed. The single manually resolved workflow path passes its scoped check;
  upstream research assets were not altered to silence their whitespace.
- The first post-merge shell checks exited 127 and 2 because sparse checkout
  excluded the two tracked helper files after the merge. Their Git blobs existed
  throughout. Sparse paths were added, the source hashes rechecked and only the
  affected shell checks rerun. Initial and retry receipts are both preserved.

This completes local integration gates. Current-head CI and live repaired official
publication remain separate requirements. v0.1.2372 is still recorded as failed
binary checksum admission; no replaced public assets or new frontier claim follows.
The companion editable diagram and inspected PDF describe the gate sequence;
there is no scientific performance series or changed canonical mathematical graph.

The helper transport tests now have a dedicated release-assets CI job on PRs as
well as main; the official release depends on it. This wires existing local
checks into remote admission without changing the root-library test closure.
