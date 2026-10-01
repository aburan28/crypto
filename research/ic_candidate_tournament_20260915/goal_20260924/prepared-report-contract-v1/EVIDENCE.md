# Prepared report repair and source-witness evidence

The prepared F5 native report now carries `query_schema_version: 1`. Both its
complete and incomplete Rust controls assert that contract. The autolab CI
workflow explicitly runs all four ignored prepared Rust unit controls and
passes the rebuilt release CLI's actual reports to the independent Python
mathematical auditor. Previous Python controls repurposed a cold report, which
already had the header; they could not catch the real prepared serializer bug.

[PROTOCOL.md](PROTOCOL.md) was committed at `6840c014a` before the controls.
The Native CLI tests use only the disclosed n17 point, accepted preparation and
existing repeatable correctness seed `2026093032`. These do not reopen either
consumed registration from PR #1122 or sample a fresh point.

The local ARM64 macOS controlled build passed both actual CLI/auditor cases:

| Repeatable fixture | Native termination | Mathematical replay | Attempts retained | Missing-header negative control |
| --- | --- | --- | --- | --- |
| one-attempt cap | incomplete, exit 2 | incomplete, no scalar | 1 | rejected |
| eight-attempt cap | complete, exit 0 | independently verified scalar | 3, including two failures | rejected |

Source and executable gates passed before/after both invocations. The build
manifest SHA-256 is `ba8bbb6db71515344f6bf57046770cb936463804822a3f602d1860fb4d8a385f`;
worker SHA-256 is `fc168e3fb380a8354573316707cec66f2ee0d7ac789795f8def49f5b059e5df5`.
These are controlled native-build/interface checks, not a scientific `IC1` run
with complete preexecution Python/interpreter attestation. Raw timer fields
remain fixture diagnostics. No online speedup, fresh qualification or promotion
is established. Registry dependency file hashes are in the build manifest;
full registry source bytes still need the later versioned retained-input gate.
The old native input v1 source pin and runtime v2 are unchanged and cannot admit
this corrected worker.

[native-controls-macos-arm64/evidence.tar.gz](native-controls-macos-arm64/evidence.tar.gz)
retains the actual build, root source archive, worker, raw complete/incomplete
reports, process receipts, preparation, jobs, mathematical audits and summary.
Its [receipt](native-controls-macos-arm64/receipt.json) identifies every byte,
mode and digest. [publish_controls.py](publish_controls.py) captured the terminal
evidence without executing a solver. Extract into a new directory, verify the
archive and every member against that receipt, and replay each raw report through
`audit_native_target` using its retained job and preparation. Verify native build
binding against the retained worker without executing it on another platform.
CI is configured to execute the same CLI controls against a separately rebuilt
Linux worker. Require the linked PR's final exact-head CI success before treating
those controls as accepted. This does not establish performance or physical CPU
coverage.

The [SAT source check](sat-source-check/result.json) independently substitutes
the retained geometric witness `[29,51,2]` for query `[62577,27783]` into the
unchanged source exports from PR #1122's sealed original archive. All 50 ANF
equations, 2,364 CNF clauses and 50 XOR rows pass; all 716 auxiliary AND variables
have explicit source definitions. The complete 767-bit assignment lifts back to
the exact query point. Deliberately flipping a primary source bit or an auxiliary
AND bit rejects in the respective checker. No native process or search ran.

This proves that the particular SAT query that exhausted its conflict budget
has a source-valid lifting solution. It does not prove the encoder globally,
change the original budget-inconclusive result, recover the original IC target
through its native SAT solver, estimate natural yield or qualify a speedup.
The source files and assignment are retained beside the result for independent
replay. `check_retained_sat_source.py` is restricted to this sealed toy invocation
and original published witness; it generates no new query or target.

Next work is a versioned corrected-worker asset/admission adapter and newly
frozen complete runtime controls, plus bounded SAT solver/encoding diagnostics
on the source-valid retained case. Natural collection and failed-attempt evidence
stay attached to the original preparation certificates. Fresh exposure/reference,
comparison-host build, calibration/resource/order and paired-target gates remain
pending. All three old confirmations and both prepared registrations stay closed;
the full goal is not complete.
