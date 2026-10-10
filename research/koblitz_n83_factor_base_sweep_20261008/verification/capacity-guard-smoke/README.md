# Capacity supervisor synthetic guard check

`fixture.c` is a public synthetic Linux worker, built with Zig 0.16.0 for
`x86_64-linux-musl` (SHA-256 of the final executable:
`9cb6e06e0325769aea8ca1f6b442f2c8fbe7cdd4d3b1c418cec65139aec5568e`).
It does not import a factor base or construct an S4 model. K=64 writes a
structurally valid worker receipt using the container's `memory.max`; K=256
sleeps past a 0.5-second wall cap; K=600 touches 256 MiB under a 128 MiB
memory cap. The final-source fixture was used for all three retained cases.
The supervisor pinned the local container image to
`sha256:0dd364ba7e10242f07755449e3a3d0e35f9efd987952737b90def6709ab0c5ce`.

`pass/outer.json` records `PASS_model_construction_only` for the synthetic
worker, `timeout/outer.json` records `UNKNOWN_wall_cap` with the container
removed, and `memory/outer.json` records `UNKNOWN_resource_or_worker_exit`
with worker exit code 137. A separate `malformed_fixture.c` writes a JSON
array where an object is required; `malformed/outer.json` records
`PRODUCER_FAILURE_worker_receipt` and preserves the parser reason. Every case
used Docker's `--memory 128m`,
`--memory-swap 128m`, `--network none`, and a locally cached image; a direct
container check read `memory.max = 134217728` and `memory.swap.max = 0`.
The `PASS` label applies only to this guard fixture, not to N83 model
construction or an index-calculus runtime.

The original exact config, outer receipt, worker receipt when available,
stdout and stderr are retained under each named directory. The supervisor
requires a new output directory and refuses a non-Linux binary; its native
worker independently checks the cgroup ceiling before importing retained
points. The first local attempt to use macOS `RLIMIT_AS` failed with `EINVAL`,
so the committed implementation uses the verified Docker cgroup path.
