# Checkpoint-compatible sigma deployment

The live campaign has about 48.1 GB of checkpoints and 6.1 GB of
distinguished-point objects. Its walk and frozen geometry remain:

```text
WALK_TABLE=0 BATCH=16 THREADS=256 MINBLOCKS=2
```

Only `PACKED_INLINE_POLY=3` is selected. Six alternating pairs improve from
14.414191 to **14.854760 B/s** (1.03056×), with every pair favoring the
candidate. The candidate replayed 300/300 reports with zero dropped.

The broader table-walk optimization bundle regressed to 14.387 B/s on the
sigma walk and is not deployed. Table walk, batch 17, 512-thread geometry and
halving all remain excluded because they would orphan or mix the live corpus.

`build.sh` accepts `INLINE_POLY=3` so this exact compatible build can be staged
without moving the live pointer. `comparison.log` and the build/verification
logs retain the evidence.

## Deployment

Activated `bin/f3ab725a533b777e` as kernelVersion 3. All eight active GPU
walkers adopted it and report an aggregate 40.52 B/s. No walk, geometry,
checkpoint or corpus field changed.

One legacy CPU worker cannot report kernelVersion adoption, but the staged and
live host binary hashes are identical. It was deliberately left running
because its legacy slot has no committed checkpoint pointer; terminating it
solely to clear the rollout counter could discard progress.
