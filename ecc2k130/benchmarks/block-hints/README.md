# Block-compacted v3 cycle-hint screen

This bounded RTX PRO 6000 experiment keeps the selected B16/T512
split-forward, reconverged-hints and exact-fast2 kernel fixed. The candidate
adds `TABLE_BLOCK_HINTS=1 TABLE_HINT_QUEUE=512`.

Every thread first records its hinted slots exactly as the control does. The
block then atomically assigns those `(owner thread, slot)` pairs to a 512-entry
shared queue. When the queue fits, consecutive CUDA threads resolve consecutive
entries, concentrating sparse cold probes into dense warps. A full queue falls
back to the measured per-lane mask schedule. All owners wait at a block barrier
before reading their resolved tags and building the original prefix chain.
Queue order cannot affect a result because slots are independent and each
history word has one writer.

The experiment requires both arms to replay 300/300 reports with zero drops,
identify the exact B16 geometry and feature markers, and produce identical
header-aware table-v3 corpora before timing. Two warmups are excluded. A
control-candidate-control screen must clear the faster control by 0.5% before
three alternating 64-launch pairs are admitted. Instrumented or overflowing
rows are not substituted for this fixed queue-512 candidate.

```sh
modal run --detach modal_job.py \
  --job benchmarks/block-hints/gpujob.sh \
  --out /tmp/ecc2k-block-hints --gpu RTX-PRO-6000
```
