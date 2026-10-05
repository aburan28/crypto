# F5 one-target registration: published, not dispatched

This is the retained **pre-dispatch** registration record. The subsequent
sole execution and original frozen audit are in `result-v1/RESULT.md`.

The scientific target capsule was frozen from committed source
`5e4a4d8ea19bae3afa2398b7d9a310e069cae151` and completed full offline
source/build verification. Its external registration seal is
`fcc320ea595781d94ab9cfc9461759c8f639fc58505489c1c8e4715e28be7d2e`.
The exact public point `[61889,74818]`, seed 2026100501, eight-query limit,
900-second worker limit and original consumed F5 preparation are frozen in the
capsule. It is scientific (`validation_only=false`) but has made **zero**
target worker calls. The earlier build-only parser failure is retained in
`FREEZE_FAILURE.md`; it consumed no target registration.

`publication-v1/` contains the complete 6,951-member target capsule archive,
four original sidecars and all build receipts. The archive SHA-256 is
`8a23af98b6afc104ccfb9eaf0c550746b6b6832faafdfc67fa06ee9aa0681f42`.
The source inventory SHA-256 is
`05b7672b5e7f66c59c6fd7ecaf36bd713efa11b90cd48d3a200794d517b17473`;
the compiled worker and checker hashes are
`1b46e7157bb790a56f3857bfb32d242757bb8264b0d24fefa10dfdaa110acac9`
and `27b4494a0cbd8b55bae7e326fb38b66492f16846ce6d8cb34cfdafd9b4bcc147`.
Publish-time replay and a separate second data-only replay returned byte-identical
`PASS_DATA_ONLY_TARGET_SCIENTIFIC_REGISTRATION_CUSTODY` receipts, with all
6,951 archived regular files verified and zero archived binaries executed.
The retained receipt is `publication-replay.json`.

The original capsule remains at
`/Volumes/SSD990/llm/tmp/ic-native-f5-one-target-registration-20261005-capsule-v1`.
The pending sole execution must use that capsule's `immutable/bin/icprog`, this
publication through `--publication`, a new execution directory, and the seal
above. The controller must verify the published archive *before* it consumes
the claim. Do not retry or refill this registration on failure.

This registration proves data custody, not a DLP recovery. Target attempts,
online wall time, correctness, candidate/workload/run IDs, rho pairing,
speedup and fresh-target qualification are all **unknown or null** at this
stage. CPU timing from this ordinary Mac will remain exploratory even after
a verified target audit.
