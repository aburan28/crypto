# Bounded construction gate for a larger N83 factor base

`v2_object_supervisor.py` wraps the existing one-object exporter with a Docker
cgroup memory cap, equal memory/swap limits (zero swap), one CPU, no network,
and an outer process-wall deadline. It writes `config.json`, separate process
logs and `outer.json` into a fresh directory outside the repository. The
worker writes its object and manifest under `object/`. On timeout it kills the
named container and checks that it is gone. Worker exit 137 is censored as an
unknown resource/worker exit; it is not asserted to be an OOM event.

The v2 exporter streams JSONL through a plain-byte BLAKE3 hasher into a gzip
file named `.partial.jsonl.gz`. After gzip finishes, it hashes the compressed
file and renames it to the content-addressed object name. It still holds the
representatives and point coordinates needed for the sorted point-set hash;
peak memory for a declared larger size has not been measured. An interrupted
worker may leave the partial file, but no completed manifest names it and the
uploader requires a matching replay receipt. A small public fixture checks
that streamed and buffered serialization have identical plain and compressed
bytes.

This gate checks a clean source commit, executable SHA-256, the one-object
manifest's declared curve/policy/size/seed and destination, and basic file
presence. `PASS_construction_only` does **not** certify the compressed BLAKE3
digest, each point, or subgroup closure. The separate `v2-replay-one` command
does that with generic multi-limb arithmetic. Only its matching PASS receipt
allows the existing `upload` command to publish the object and verify a S3
download byte hash. The outer receipt keeps total-runtime and winner fields
null. The pilot's 3,600-second allocation is exhausted, so this command has
not been run on a larger N83 object.

After a new construction and replay budget is assigned, build a Linux ELF
exporter from a clean committed checkout, retain its SHA-256 and exact build
command, and choose a locally present immutable Docker image. For example,
substitute absolute paths and approved caps in:

```sh
python3 research/koblitz_n83_factor_base_sweep_20261008/v2_object_supervisor.py \
  --checkout /Volumes/SSD990/crypto/worktrees/codex-n83-factor-base-sweep-20261008 \
  --repository-root /Volumes/SSD990/crypto \
  --binary /absolute/path/to/linux/exporter \
  --output-dir /private/tmp/n83-v2-k1182-unique \
  --curve-a 0 --policy public_x_hash --columns 1182 --seed 2026100801 \
  --wall-seconds 600 --memory-mib 4096
```

The Docker image must contain the worker's dynamic libraries and `git`, and
the repository root must be accessible through Docker file sharing. The
worker uses its compiled `CARGO_MANIFEST_DIR`; its emitted source commit must
match the supervisor's frozen commit. Build provenance for the full binary is
still a separate receipt requirement. Run `v2-replay-one` with its own bounded
wall and memory guard before `upload`. Never treat a construction or replay
receipt as a measured relation yield, matrix rank, individual log, or
single-target index-calculus runtime.
