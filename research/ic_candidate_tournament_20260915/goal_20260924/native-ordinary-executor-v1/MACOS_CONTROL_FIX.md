# macOS temporary-path control fix

At source head `77488a0183cb8deb93fba5b499b1e549e29dfb71`, Linux worker
controls passed. The macOS disclosed replay failed one of 11 worker controls:
`helper_requires_new_claim_parent_durable_start_and_exact_public_query`.
[The original log](first-macos-ci-failure.txt.gz) is retained from
[run 37239177599](https://github.com/aburan28/crypto/actions/runs/37239177599).
The helper rejected the fixture's claim path; this is a failed CI gate, not a
scientific solver failure. No workflow was restarted to hide it.

The fixture returned an uncanonicalized temporary directory. On macOS, TMPDIR
can contain a `/var` alias while the helper canonicalizes it to `/private/var`.
Production execution records already store canonical paths. Canonicalize the
test fixture after creating it; preserve the strict production claim check.
No source-runtime gate, algorithm, mathematical outcome or old registration
was relaxed or changed.

All 11 worker/transport controls passed with an explicit symlinked TMPDIR:

```sh
TMPDIR=/private/tmp/ic-ordinary-contract-alias-20261004 \
CARGO_TARGET_DIR=/Volumes/SSD990/llm/tmp/ic-native-preparation-build-cache-20261004 \
CARGO_INCREMENTAL=0 CARGO_PROFILE_DEV_DEBUG=0 CARGO_PROFILE_TEST_DEBUG=0 \
/Volumes/SSD990/llm/tmp/ic-native-busy-20261002 busy -- \
  cargo test --locked --offline --bin prepared_ordinary_worker -- --test-threads=1
```

The alias points to this task's new directory
`/private/tmp/ic-ordinary-contract-alias-target-20261004`. The exact test output
is retained in [alias-control-tests-v1.txt.gz](alias-control-tests-v1.txt.gz).
Formatting and diff checks passed. Final-head Linux/macOS CI remains required.
These are source/transport controls, not a new scientific query or runtime
admission. Old target controls and confirmation sets stay closed.
