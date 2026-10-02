# GPU-wide exact-v3 hint queue: frozen first implementation

Source parent: `60fb2ae47fcca9352d8bae2b280fce77371e57f0`.

The active goal is 26 billion complete scalar updates/s on one RTX PRO 6000.
The exact table-v3 path currently measures about 5.10 B/s. The reference's
block queue packs rare hinted owners into one or a few resolver warps, while
every other warp in that block waits at a barrier. This experiment separates
that cold pass from the hot Montgomery arithmetic and packs hinted owners
across the full GPU. It preserves the exact v3 map and raw-cycle anchor.

## Candidate and accounting

`TABLE_GLOBAL_HINTS=1` selects four ordered stream operations per scalar step:

1. Clear one device queue counter.
2. Select raw tags and perform the existing DP/guard checks for every slot.
   Each worker reserves one contiguous queue range for its pending slots.
3. A GPU-wide resolver grid calls the unchanged `tableResolveHintSlot` on each
   queued `(slot,worker)` exactly once.
4. Run one unchanged prefix/inversion/reverse update from the resolved tags.

The queue capacity is the full walk population, so it cannot overflow: each
slot emits at most one hint. Keys encode `slot*workers+worker` in a 32-bit word;
setup rejects a population above UINT32_MAX. Kernel boundaries provide the
global ordering that a block barrier cannot. Queue state is transient and is
rebuilt before every update; it is not checkpointed.

The first prototype uses ordinary stream launches, not a CUDA graph. All
counter resets, queue traffic, shared-table copies, host launch overhead,
synchronization and report processing remain inside the complete-update rate.
The experiment predicts no rate from its sparse operation count. A graph may
be a later separately registered scheduling follow-up if launch cost dominates.

## Frozen comparison

One NVIDIA RTX PRO 6000 Blackwell Server Edition, sm_120, CUDA 13.3.73.
Both arms use the current exact-v3 B16/T512/min1 preset with fast2, shared
tables, compact state, square table and polynomial inversion mode2. The control
uses block hints. Candidate disables block hints and selects GPU-wide hints;
the resolver uses 128-thread blocks and one block per SM as its initial grid.
Both use 96,256 workers for correctness and 385,024 for timing.

Before timing, require native field/queue controls, matching source/flags,
successful native CUDA builds, odd launch lengths, 300/300 reference replay,
zero drops, deterministic spread replay with nonzero trail coverage and
identical sorted complete v3 corpora. Checkpoint continuation must match in
both directions. Queue reset, empty/full queues, partial blocks and repeated
launch boundaries receive explicit controls.

Timing starts only after every gate passes. Exclude one warmup per arm, then
five alternating A/A pairs and five alternating A/B pairs, each completing
`385024*16*1024*32` scalar updates. All logs must bind exact final work and
feature/resource markers. A/A symmetric drift must stay below 1%.

## Decision and stop

Goal achieved only with a correctness-verified candidate median at least
26,000 M complete scalar updates/s on that one GPU. An engineering candidate
qualifies at geometric-mean paired ratio at least1.10 and every pair above1.0
under the noise gate. Otherwise preserve the negative result without selecting
the mode. Any correctness/identity/overflow failure suppresses timing. This
first bounded panel ends after the five pairs; no search or solver is run.

## Native producer details frozen before dispatch

`gpujob.sh` runs native packed field and strict v3 replay controls, then builds
both arms and the production queue device control with the same preset flags.
The device control uses 64 public fixture points, including the established
scalar1184 inverse-edge pair. It compares production selected/resolved histories
against scalar-reference tags, checks every queue owner exactly once and both
allocation canaries, and reuses allocations through empty/full/empty/sparse/
endpoint/full/empty patterns at1,127,128,511,512,513 and1537 workers. Dead flags
suppress report emission only in this queue control; full-walk controls retain
the ordinary DP and reseed behavior.

Whole-walk controls use95 steps for7 launches, run-id7 and DP weight48. The
96,256-worker corpora require300 replayed reports per arm and independent
evenly spread300-record replay with at least299 nonzero trails. Separate
511/513-worker device runs require the same300 replay and corpus identity.
For checkpoint continuation, each arm produces a513-worker prefix of4 launches,
then both arms resume copies of each prefix for3 launches. All four final
checkpoint payloads and sorted continuation corpora must be byte-identical.

Each arm's one excluded warmup completes the same32 launches as a ranked row.
The fixed timing schedule is five A/A pairs followed by five A/B pairs, with
odd pairs control-first and even pairs candidate-first. Native summary rejects
missing/reordered rows and invalid rates/digests; producer separately requires
the exact201,863,462,912 update count, zero drops and feature/resource markers
in every timed log. Binaries and helper executables are retained in the raw
archive. Any failed build/control/sample exits nonzero and remains an attempt
artifact; no failed attempt may be used as a throughput result.

Require exactly one visible GPU and a clean40-hex source revision. Runtime
markers bind the complete frozen backend/feature identity and188 resolver
blocks. Each binary's first full verification log records its hot/select/resolve
resource tuple; every subsequent log must match that tuple. Local stack is
recorded and charged, including the unchanged cold anchor's legitimate local
array, rather than admitted through an assumption of zero local bytes. The
native summary requires the exact preflight marker and records its SHA256.
The producer manifest covers binaries, corpora, checkpoints, logs and source
manifests, excluding only itself and the wrapper's still-growing job.log and
later-added exit-code; an independent post-run manifest/archive binds those.
