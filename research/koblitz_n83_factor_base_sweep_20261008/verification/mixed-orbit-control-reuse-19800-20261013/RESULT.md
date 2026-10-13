# Known-answer solver-access control after the 19,800-second ledger: source only

The first cross-stream K600/K600 natural-target call imported two disjoint
signed-Frobenius bases, then reached the 5,000-conflict cap on its first sign
pattern without a model. The earlier nested K600/K2048 and K1182/K2048 calls
also reached their caps. These observations leave the solver's access to the
known five-point witness unanswered. The prepared control uses the pinned
K1182/K2048 fixture and the previously built Linux worker to test that
specific solver-access question; it is not a measurement of natural relation
yield or complete index-calculus runtime.

The historical [17,400-second control gate](../mixed-orbit-control-reuse-17400-20261012/RESULT.md)
is preserved unchanged. This new `reuse_control_supervisor.py` requires
authorization schema `n83.additional-local-budget/v7`, purpose
`one_known_witness_control_reuse_k600_worker_after_cross_stream`,
`prior_total_charged_seconds="19800"`, and the exact hashes of both the
K600/K2048 audit and the [cross-stream audit](../cross-stream-first-gate-source-20261013/audit.json).
The new audit SHA-256 is
`dfc96f11e4a399a42559847dd39e94b753f37911a54d806e0e54e830e7735478`.
The K600-built Linux ELF remains pinned at SHA-256
`26dcb2b49388fcb89853eabbda01db71d148620cd465bd5458099e0606dcf606`.
The supervisor checks the historical K600 17,400-second audit and the newer
2,400-second cross-stream charge separately; it cannot mistake the older
build receipt for evidence of the current cumulative ledger.

Read-only S3 HEAD checks found the K1182 object at 5,956,114 compressed bytes
and the K2048 object at 10,313,647 bytes, matching their committed manifests.
The second object's ETag is multipart, so HEAD metadata is only an availability
and length check; a future launch must download both objects and match their
pinned compressed SHA-256 values before import.

The proposed additional conservative charge is 1,200 local seconds: 120 for
read-only `verify-control`, 900 for one `control-run` SAT call, and 180 for
cleanup and audit. A fully launched run would bring the ledger to 21,000
seconds. The worker watchdog is 840 seconds, with one CPU, 4 GiB, zero swap,
read-only inputs and no network. Both retained S3 archives, the fixture,
binary, source and prior audit hashes must match before the supervisor creates
an output directory. A separate fresh user approval reference and fresh
output path are mandatory. No authorization record has been issued for this
gate.

The separate `audit.py` is now included among the frozen supervisor inputs.
It independently checks Git source blobs, authorization and archive hashes,
the exact offline Docker commands and cleanup, worker counters, the complete
1,200-second conservative charge, and any candidate against both archived
point domains and a Python group sum. It records raw-file SHA-256s and uses
`UNVALIDATED` for incomplete or inconsistent receipts. A conflict cap remains
`VALIDATED_CENSORED`, while only a group-replayed solver candidate can receive
`VALIDATED_CONTROL_ACCESS`; neither status is a natural-yield or runtime win.
The auditor exits successfully only for that verified control-access result;
censored and failed results still retain their raw hashes and explicit status.
The auditor has not yet processed a retained control run.

Seven focused source-only tests passed. They reject the stale v5 grant before
launch, changed prior audit bytes, a changed conflict counter, network access,
and a success claim without archive/group replay; the constructed fixture
also replays in the independent Python group arithmetic. `--help` and Python
compilation pass. No container, archive import, or SAT search was launched
here. The current conservative ledger remains 19,800 seconds and the
selected base remains undetermined.
