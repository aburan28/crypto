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

## Measured result

Both arms replayed 300/300 reports with zero drops and produced identical
sorted multisets of 1,480,278 unique v3 records, SHA-256
`2fc6f84d4b1969240418cf6f8e779df9538b1d8f01eb6ff4c6a00be7324fae01`.
An independent spread replay added 300 records per arm, with 298 nonzero trails
and maxima of 572/575 table steps. Both kernels used 128 registers, a 400-byte
frame and zero spills. The candidate adds 1,040 static shared bytes beside the
48,732-byte dynamic table.

The screen measured 4.134526 / 5.470673 / 4.144668 B/s. Three alternating
64-launch confirmation pairs were:

| Pair | Per-lane control B/s | Block queue B/s | Ratio |
|---|---:|---:|---:|
| 1 | 3.558739 | 5.019737 | 1.410538 |
| 2 | 3.559925 | 5.015825 | 1.408969 |
| 3 | 3.563098 | 5.019275 | 1.408683 |

The candidate median is **5.019275 B/s**, versus **3.559925 B/s** control;
the median paired ratio is **1.408969**. [`result.json`](result.json) retains
the full run, and [`independent-audit.json`](independent-audit.json) checks the
raw logs, spread replay, queue semantics and overflow/partial-block model. The
61,214,341-byte archive has SHA-256
`2acc1465a5248c1ea0816a4cd5f68f2035fa4e824d7656f710a349aebfa75ad9`.

The first frozen attempt failed before timing with CUDA `invalid argument`:
the static queue reduced the default dynamic-shared allowance. Commit
`0aad485b` opts into the unchanged table allocation. That failed archive is
retained separately with SHA-256
`432f66c05940040801d0241dfa451b842e0cdfeef426e6c872b867457f83c730`.

```sh
modal run --detach modal_job.py \
  --job benchmarks/block-hints/gpujob.sh \
  --out /tmp/ecc2k-block-hints --gpu RTX-PRO-6000
```
