# Validation record

The first supervisor test compile failed because the included repository SHA256
tests require the `hex` crate. Added `hex` as a development dependency; this does
not change the launch executable's dependencies or map acceptance.

Original diagnostic:

```text
error[E0433]: failed to resolve: use of unresolved module or unlinked crate `hex`
  --> ../../src/hash/sha256.rs:342:9
342 | hex::decode(hex_str).unwrap()
error: could not compile `additional-prime-isogeny-campaign` (bin "isogeny-campaign" test)
exit status: 101
```

The updated CLI passes 63 tests plus 3 field-reference tests. The frozen running
supervisor passes 14 tests. The later archival utility passes 16 tests, including
rejection of alternate decimal encodings and unreduced coefficients, and strict
Clippy with `-D warnings`. Changed CLI and research source passes rustfmt 1.98.

The installed clean binary exactly replays two saved steps of a 1050-step
request and reports `complete: false`. The saved native mutation controls pass:
changed source/model, changed mapped generator with recomputed hashes, and empty
chain are rejected. Positive two-step replay passes. Extracting the archived
runtime outside the checkout also passes a complete two-map P384/19 construction
and independent replay. These control maps are excluded from campaign yield.

The former root release library process and scratch log were lost before a
completion receipt could be retained. A fresh required gate now runs with its
log under `runs/restart-20261011/validation/`. Historical merged-tree validation
does not replace this check. Saved native test/control logs remain under
`evidence/validation/` and `evidence/controls/`.

The first commit included the CLI changes but sparse-checkout staging omitted
the new study directory. The following protocol commit explicitly stages the
complete study and lockfile before any measured campaign execution. The clean
producer revision is recorded in LAUNCH.json.

The first mutation-control command referenced an absent cached executable and
returned 127 before replay began. Its diagnostic is retained. Rebuilding the
saved native control program corrected the command, and the controls passed.
The first archive smoke command used a checkout-relative archive path after
changing to `/private/tmp`; extraction and invocation returned 127. Empty stdout
and stderr were saved with `wrong-path` filenames before the corrected absolute
path command passed. Neither failed invocation contributes maps or timing data.

The initial PDF has two pages, opens successfully, and was rendered and visually
reviewed. Its counts are a frozen snapshot: 6 P192 and 19 P224 repeated steps,
0 additional maps. Live status is separate and the 2000-map requirement stays open.
