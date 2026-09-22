# Campaign seed ownership

The seed layout is `run_id << 48 | lane_index << 16 | restart_counter`.
Run IDs therefore identify disjoint seed ranges only while the lane index
fits in 32 bits and the restart counter fits in 16 bits. The client refuses
seed-counter exhaustion instead of carrying into the adjacent lane. Unique
seeds prevent duplicate starts; independently seeded walks can still meet,
which is expected in collision collection.

Storage directories, provider names, GPU types and process IDs do not create
new seed ranges. Every curve-131 production launcher must consult the same
campaign bucket's `seed-registry/curve-131/` records before starting its client.
Moving to another bucket requires a campaign migration, not a new independent
registry for the same collection.

`aws/seed_registry.py` enforces:

* A run ID has a permanent stream binding. A different S3 prefix, slot backend
  or Modal run cannot claim it. Old records are retained after shutdown.
* Acquisition uses S3 conditional writes; competing starts have one winner.
* An active owner never expires automatically. An unreachable process might
  still be computing. Recovery must establish that it stopped before clearing
  its owner, even if its provider lease or nominal deadline has expired.
* A used run requires a full checkpoint for that curve and run ID, at least as
  advanced as the recorded checkpoint floor. A 40-byte status header is not
  a resumable checkpoint. The client separately validates complete contents
  and build geometry before walking.
* Missing registry access, missing checkpoints, conflicting identities and
  ambiguous ownership prevent a launch. If a supervisor raises while its child
  is still alive, it retains ownership.

The registry is at the bucket root, outside each producer's storage prefix.
The legacy S3 worker, the prefixed cloud runner and the standalone Modal GPU
and CPU collectors use it. Modal GPU IDs are 8000–8999 and CPU sidecars are
9000–9999; GPU-only calls cannot take sidecar IDs. `next_run_id` consults both
the volume and permanent registry, including a CPU partner's prior use. It is
a suggestion; the conditional acquisition at launch remains authoritative.

The Modal driver preserves its call IDs and stops for recovery after a failed
status read, cancellation, missing checkpoint or unsuccessful client exit.
It never treats an uncertain call result as permission to start another GPU.

## Initialization and deployment

`seed_audit.py --bucket BUCKET --out audit.json` inventories root and prefixed
slot authorities, checkpoints, and historical upload directories. Unmapped
directories are decoded record by record instead of inferring IDs from their
names. Review active owners and historical conflicts before initializing with
`--initialize`; readiness is written last, and repeated initialization cannot
overwrite a binding, clear a live owner or rewind its checkpoint floor.

The 2026-09-21 repair imported 225 historical IDs and 11 active owners. The
eleven active assignments were distinct. Nine older upload directories needed
full decoding: 27,524,849 records, 880,795,168 bytes. Their otherwise untracked
IDs are retained and quarantined. The prior run 1–4 Modal/AWS conflicts remain
recorded, with the original slot bindings preserved.

A live control in a separate `checks/seed-registry/` namespace raced eight
claimants against actual S3 conditional writes: exactly one acquired the ID,
seven were refused, and conditional release plus checkpointed resume passed.
It allocated no production run IDs and started no GPU clients.

The shared bucket's published worker received only the seed-guard wrapper;
its prior implementation was saved under a content-addressed backup key.
Older bootstraps that download only `worker.py` can fetch the hash-pinned guard
before starting a client. Updated bootstrap scripts install it directly.

Both Modal launchers were deployed using existing production images. The
`ECC_SEED_BASE_IMAGE` deployment option layers the supervisor guard over that
image, preserving the client binaries and checkpoint formats. No replacement
GPU campaign was launched as part of this repair. Already-running processes
keep their existing code; their imported owner records protect their IDs from
new guarded launches. An older image or manually invoked binary that bypasses
the guarded launcher is not covered by this mechanism and must not be used
for campaign collection.

## Recovery after a crash

Do not delete the identity record or reset its progress. Read its active owner,
identify the exact old process/container, and establish that it has stopped.
Retain the newest compatible full checkpoint and corpus. Only then release
that exact owner using `SeedRegistry.release(run_id, owner, checkpoint_path)`;
the conditional update preserves the binding and maximum progress. If no
matching checkpoint survives, retire the ID and allocate a new one. Never
reconstruct its old seed range from iteration zero.

Imported `legacy:` owners deliberately require this same review when those
pre-guard processes finish. The tradeoff is explicit: an uncertain worker
stops automatic recovery rather than authorizing potentially duplicated work.

This is an operational ownership guarantee for cooperating guarded launchers,
not a claim that different seeds produce different curve points, or a change
to the mathematical collision-search complexity.
