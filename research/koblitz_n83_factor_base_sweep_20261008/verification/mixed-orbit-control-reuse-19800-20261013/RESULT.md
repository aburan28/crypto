# Known-answer solver-access control after the 19,800-second ledger

The authorized v7 control imported the retained K1182/K2048 bases and passed
the read-only fixture verification. Its single known-answer SAT call then
exhausted 5,000 conflicts on the first sign pattern in 471.165 seconds of
worker wall time. It returned zero models and no candidate. The independent
auditor classified the result `VALIDATED_CENSORED` (`UNKNOWN_search_cap`),
not a refutation of the known five-point witness. This control measures solver
access to that witness only; it supplies no natural relation yield, rank,
verified column logarithms, or complete index-calculus runtime.

The [retained evidence](evidence-01/run/audit.json) contains the authorization,
frozen input inventory, raw verify and control outputs, Docker command and
container states, budget, and independent audit with SHA-256s of every run
file. The supervisor ran from clean source commit
`bc9c162d1e3efaaded0e5f9c07d215302eecdf83` and reused the Linux ELF at
SHA-256 `26dcb2b49388fcb89853eabbda01db71d148620cd465bd5458099e0606dcf606`.
The pinned image ID was
`sha256:7c4ae649a84014c467d79319bbf17ce2632ae8b8be123ac2fb2ea5be46823f31`.
Both retained S3 archive downloads matched their compressed SHA-256s before
launch; the K1182 and K2048 objects were 5,956,114 and 10,313,647 bytes.
The container had one CPU, 4 GiB, zero swap, no network, read-only inputs,
and an 840-second worker watchdog. Its measured cgroup peak was 442,339,328
bytes; base import took 7.143 seconds. The SAT model had 101,636 peak
variables, 714,072 peak clauses, 4,980 peak XORs, and 465,738 domain clauses.
The supervisor wall time was 486.084 seconds; both containers exited cleanly
and their cleanup was verified.

The user approved up to 1,200 additional local seconds after the audited
19,800-second ledger and permitted a Docker Desktop restart. Docker's API was
unavailable before launch; it restarted successfully and the exact pinned
image identity was rechecked. The frozen authorization is preserved in
`evidence-01/authorization.json`. The conservative full charge is 120 seconds
for verification, 900 for the one search, and 180 for cleanup and audit:
**1,200 seconds**, bringing the cumulative local ledger to **21,000 seconds**
with zero remaining under this approval. Measured wall time is separately
reported above and does not reduce the conservative charge. The historical
[17,400-second control gate](../mixed-orbit-control-reuse-17400-20261012/RESULT.md)
and the [cross-stream audit](../cross-stream-first-gate-source-20261013/audit.json)
remain unchanged. The v7 supervisor pins both prior audit hashes separately.

`audit.py` checks the Git source blobs, authorization and compressed archive
hashes, offline Docker command, cleanup, worker counters, budget and any
candidate against the archived point domains and a Python group sum. It exits
with code 1 for this valid but censored result; its audit JSON is the evidence
classification. The source-only focused tests passed before the run. No
factor-base winner is selected: complete, matched one-target cold runs and
independent replay remain necessary for a total-runtime comparison.
