# Interrupted counted main attempt

The first execution of `specs/counted-main.json` stopped when its controlling
agent session was interrupted. The runner process was gone and the shared
benchmark lock had no holder when inspected on 2026-10-09. Its manifest still
says `running`, with no finish time or records digest; the JSONL has 208 of
the 288 planned executions. The directory was moved intact to
[`interrupted-counted-main-208`](interrupted-counted-main-208), outside the
CI path for completed sessions. No result from this attempt enters a
comparison. The frozen spec was started again in a new directory.

| file | SHA-256 |
|---|---|
| `interrupted-counted-main-208/session.json` | `9ab45974bc7cdfa7be6714360d009cddc975e46d59e51cad7efc1a4fa25b019a` |
| `interrupted-counted-main-208/records.jsonl` | `ed56203f059315dfaf91d49a8810095879aa0c9c3c5c4dc07c5e17e84aacda38` |

The incomplete manifest is preserved as written. It cannot pass a complete
session audit, and its partial records must not be silently pooled with a
finished run.

The restart was stopped deliberately as soon as the pilot's missed
per-call conflict-budget gate was found. Its SIGTERM handler wrote an
`interrupted` manifest and killed the active child; 80 of 288 planned
records were sealed, 79 verified and the last recorded a signal exit.
Its directory is
[`interrupted-counted-main-restart1-80`](interrupted-counted-main-restart1-80).
The audit with `--allow-interrupted --replay-all` reports zero problems and
31 of 31 measured deterministic replays identical. Its receipt is
[`interrupted-counted-main-restart1-80.audit.json`](interrupted-counted-main-restart1-80.audit.json),
SHA-256 `8fc96ed55d36b8deaaa52020510b0fc4bea55c4e411738d412af861079f2e4e4`.

| file | SHA-256 |
|---|---|
| `interrupted-counted-main-restart1-80/session.json` | `0b8aac5e70be6d396ff19aaf725ff288891fe34cd7757f5631bbc9a9d555e1ed` |
| `interrupted-counted-main-restart1-80/records.jsonl` | `c6d408e5e5afc3f906612af55b780049120b11d6ff8797dcd06dba12036831a0` |
