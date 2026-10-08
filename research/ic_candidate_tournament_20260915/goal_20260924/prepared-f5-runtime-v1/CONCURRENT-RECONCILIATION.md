# Accepted registrar preserved

PR #1096 added `prepared_f5_runtime_v1.py` and its original panel/asset schema
while this retained-input and frozen-audit implementation was in review.
The latter is now `prepared_f5_runtime_v2.py` with a separate panel schema.
The accepted version-one runtime remains unchanged. Neither registrar has a
new production invocation frozen or executed by this follow-up.

The first local combined prepared suite after reconciliation ran 42 tests in
36.850 seconds, with one test failure. The incoming version-one wrong-platform
test constructed a macOS fixture and assumed the host was Linux; on macOS it
therefore reached the explicitly mocked meter, which raised `worker ran`.
No native process or IC query executed. The test now fixes the observed host
to Linux/x86-64 before checking rejection of the macOS fixture. The runtime
platform guard and accepted panel/asset schema are unchanged.

The earlier PR head `eec8b19b39584097b2e6d1bce734f4ef5bf9efe4`
passed all 11 applicable CI checks, including 423 harness and 16 boundary
controls and the strict prepared-target/mixed-adapter path. That result
does not establish validation of the reconciled head; its local checks and CI
must pass independently before merge. All measured campaign jobs were skipped.

The corrected combined local suite passed all 42 tests in 31.259 seconds.
The earlier head's [Linux integration](https://github.com/aburan28/crypto/actions/runs/36795770120)
and [scaled producer controls](https://github.com/aburan28/crypto/actions/runs/36795769926)
are retained separately. The latter passed four isolated native unit controls
in 28.76 seconds and the mixed-adapter replay: 90 verified native/profile pairs,
12 expected incomplete profiles, 102 trial slots and 162 execution keys.
Those counts are wiring/correctness evidence, not a new qualification or win.
