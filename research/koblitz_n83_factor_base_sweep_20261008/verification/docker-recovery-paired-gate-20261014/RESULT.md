# Controlled Docker recovery for the frozen K600 pair

Status: `PRODUCER_FAILURE_DOCKER_DESKTOP_RESTART`. The user approved one
controlled Docker Desktop restart before the frozen K600 Gray/hash paired
gate. The daemon socket returned HTTP 503 with `Docker Desktop is unable to
start` before the attempt. `docker desktop restart --timeout 180` exited 1
after 180.293 measured seconds because Desktop's existing backend and UI
processes had not stopped by the deadline. The socket still returned the
same 503 response immediately and after twelve five-second probes.

The exact command, UTC start/end times, exit status, socket responses, and
raw-stream SHA-256 digests are in `restart.json`; `restart.stdout` and
`restart.stderr` are the byte-preserved CLI streams. Their hashes were
rechecked after copying them into this directory. The earlier VM log excerpt
in `../k600-gray-hash-paired-protocol-20261014/DOCKER_READINESS.md` records
disk I/O failure, an aborted ext4 journal, and a read-only remount. This
bounded restart did not establish a healthy daemon; no direct VM disk edit,
forced process kill, image deletion, or filesystem repair was attempted.

The pinned image `sha256:7c4ae649a84014c467d79319bbf17ce2632ae8b8be123ac2fb2ea5be46823f31`
could not be rechecked through the daemon. No run-specific grant was created,
no paired N83 solver arm was started, and no factor-base runtime comparison
resulted. The separate native import receipts are under
`../k600-gray-hash-native-import-20261014/`; they do not substitute for the
frozen Linux relation-search stages. The recovery failure is infrastructure
evidence, not an observed relation-yield result.
