# Preregistered whole-process cost control

**Status: frozen before V4 meter implementation or execution.** The
[V3 native timer run](RUN_V3.jsonl) completed on the corrected public
target order and passed [separate arithmetic replay](REPLAY_V3.jsonl).
The V3 raw SHA-256 is
`ac85539489c417981adcc41ab5f1cca00539f55269c8ec83ebe4d7231a683e10`;
its replay SHA-256 is
`47d6ee15964c390968fa9636068e4e0b7a9cb4b519a548622c17053dbff85d1f`.
Its native timer begins inside `batch-run` before JSON reading and
curve construction, but excludes solver process startup and final
stdout serialization. The V4 control measures one entire solver child
process per arm and round, using the exact same answer-free input and
the same hash-pinned solver binary as V3.

The hypothesis is that process overhead may change the local CPU
ratios at `N=1,4,16,60`. The reference remains same-prefix direct
source BSGS and negation rho. There is no selected winning prefix.
Use three rounds and the same seven-arm rotation, 60 distinct targets,
1,220 relation trials per target, 32 rho retries, 512 MiB RSS ceiling
and 120 s wall ceiling. Keep the exact V3 public/checker files, source
and solver binary SHA-256 values. A separately compiled Rust meter
launches each solver process with public input only, records child user
and system CPU from `getrusage(RUSAGE_CHILDREN)`, wall time from before
spawn through child exit and output completion, exit status and child
peak RSS. Preserve individual stdout, stderr, meter JSON and exit
files; refuse to overwrite any V4 output. The script and meter
source/binary hashes must be committed in a V4 spec before execution.

Run the existing separate native arithmetic replay over all 21 V4
solver records, requiring 1,260 scalar checks and every relation
witness. A failed child, timeout, missing meter receipt or replay
failure censors the corresponding complete-prefix cost. One table
reports child **user+system CPU microseconds** for every arm and
prefix, with paired ratios to BSGS and rho; show process wall and RSS
separately. The meter observes only the complete `N=60` child, so
prefix `N=1,4,16` values remain V3 native timer diagnostics unless
separately launched prefix-limited children are implemented and
hash-pinned. Do not label those prefix values whole-process costs.
This shared-host process measurement does not supply an isolated wall
result, second-host replay, or portable counted `S`.
