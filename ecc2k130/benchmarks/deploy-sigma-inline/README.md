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
