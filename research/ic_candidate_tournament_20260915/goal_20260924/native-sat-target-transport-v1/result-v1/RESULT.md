# Marked CryptoMiniSat disclosed transport control, original result

The predeclared nine roles ran once on three previously disclosed n17 source
instances, after the one-use natural SAT preparation had passed its original
audit. The controller used the source published in commit `d66bdb93d`, the
accepted CryptoMiniSat binary with SHA-256
`6c509f09622f103d8a3ad90afc151e1c4031275c052d7f8af481465b4f27f2af`,
and a distinct post-buffer-marked binary with SHA-256
`128f074850bbd9d973bd50be26c6ffc161e97c7a365b9234f8af171545fe7263`.
The pre-execution [registration](registration.json) has SHA-256
`426859d7f6d0527358f4dad01a948005d280c814dae6c2fddbe1c03291d8c94e`.
Each role had one CMS thread, a 1,000,000-conflict cap and a 120,000 ms wall
watchdog. No role was rerun or refilled.

| Disclosed public point | Accepted file | Accepted stdin | Marked prepared stdin |
| --- | --- | --- | --- |
| `[40991,73355]` | proved source UNSAT | proved source UNSAT | proved source UNSAT |
| `[73003,104622]` | proved source UNSAT | proved source UNSAT | proved source UNSAT |
| `[59775,2910]` | source-valid SAT model | source-valid SAT model | source-valid SAT model |

The [independent original audit](original-audit.json) returned
`PASS_DISCLOSED_CMS_TRANSPORT_PARITY` for 9/9 roles. It replayed the raw
stdout/stderr and exit codes; checked the marked READY line before any stdin
byte, exact CNF delivery, source and binary hashes, watchdogs and process
drains; verified each SAT model against the frozen ANF/CNF and a full-point
lift; and independently checked both negative three-sum classes. The auditor
made zero native solver calls. Its 1,076-byte JSON has SHA-256
`1cecc6161179445317644a92a4ecc9e13fa714f3d2dc015ea5862818f79e3da5`.

The [original execution archive](execution-original.tar.gz) retains all 49
files from the untouched controller output (53 tar entries including
directories), including each raw role, READY receipt, input copy and PID
ledger. Its 112,666 bytes have SHA-256
`15a065debdd8d5755a2e003687225f4c52f8d78dc776cf615e88958767c3b6fe`.
The original output remains at
`/Volumes/SSD990/llm/tmp/ic-marked-cms-control-execution-20261006-v1`.
The earlier 60-second file-mode timeout remains a separate failed historical
control; this result does not rewrite it.

This passes **disclosed-input transport parity only**. These points were
known before the marked binary was built and are not a natural-yield sample
or a fresh target. `source_bound_target_admitted`,
`fresh_paired_qualification` and `online_speedup` remain false, false and
null. A SAT target capsule must independently bind the exact exporter and
marked solver build to a fresh postpublication public point, then pass its
one-use source/role audit. Only a same-point, isolated-host F5/SAT/incumbent/
rho experiment can support a speed claim.
