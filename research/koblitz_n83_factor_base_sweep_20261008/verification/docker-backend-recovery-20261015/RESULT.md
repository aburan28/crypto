# Bounded Docker backend recovery before the N83 Linux gate

Status: `PRODUCER_FAILURE_DOCKER_DAEMON_UNAVAILABLE`. Docker Desktop's Unix
socket still returned HTTP 503 with `Docker Desktop is unable to start` on
2026-10-15. The earlier controlled Desktop restart and UI quit/start are
retained under `../docker-recovery-paired-gate-20261014/`. The VM log there
records write I/O errors, an aborted ext4 journal and a read-only remount.

The committed [recovery script](recover.py) first required the exact backend
parent PID `22317`, process path, parent PID 1, and the 503 socket response.
Two no-signal tests reject a changed process identity or a healthy socket.
It sent one SIGTERM to that parent and waited 30 seconds. The three observed
backend/build child PIDs exited; the parent remained in state `S`. The socket
then refused connections. `result.json` records the preflight, signal time,
process snapshots, socket responses and 31.108-second measured stage. Because
the parent had not exited, the script withheld its start command.

A separate Docker Desktop status probe then reported that Desktop was not
running. One bounded `docker desktop start --timeout 180` exited 0 in 0.327
measured seconds but printed `Docker Desktop is already running`. The socket
still refused connections immediately and after twelve five-second probes.
`post_start.json`, `post_start.stdout`, and `post_start.stderr` retain the
exact command outcome and hashes. The raw copied streams and both JSON
receipts were checked after copying into this directory.

No SIGKILL, direct VM disk edit, image deletion, or filesystem repair was
performed. The daemon never passed the health gate, so the pinned Linux
image and frozen 4 GiB/zero-swap cgroup could not be checked. No N83 relation
search or paired arm was launched. This is an infrastructure failure and
provides no relation-yield or factor-base runtime comparison.
