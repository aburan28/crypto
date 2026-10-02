# Current-source sigma inline-polynomial A/B

This is a bounded RTX PRO 6000 engineering comparison of the checkpoint-compatible
ECC2K-130 sigma walk at current `origin/main`.  It changes one build-time bitmask:

| Arm | `PACKED_INLINE_POLY` | Meaning |
|---|---:|---|
| control | 0 | single and paired polynomial helpers remain out of line |
| inline3 | 3 | both helper classes are inlined |

Every other compiled option is identical.  The fixed shape is B16/T256/minBlocks2,
385,024 workers, TILE256 compact state, weighted-prefix 2, shared sigma, native
carryless multiplication, direct reduction, generated products, polynomial state,
paired products, denominator cache, the permutation sigma network, witness off and
the table walk off.  Both arms use CUDA 13.3.1's sm_120 compiler in one Modal
container and run on one RTX PRO 6000.  Run id 7 is fixed throughout.

## Frozen question and admission rule

The measured unit is billions of complete scalar sigma-walk updates per second.
Each ranked row performs exactly

```
385024 workers * 16 slots * 1024 steps * 32 launches
  = 201863462912 complete scalar updates.
```

This is a stage throughput diagnostic.  It is not a complete ECDLP speedup, a
search, a solver run, or evidence of a cryptanalytic advance.  The standing
engineering objective is 26 B complete scalar updates/s on one RTX PRO 6000;
the result reports the candidate's ratio to that objective without converting it
to a full-solve claim.

`inline3` is eligible for promotion only when all of the following hold:

1. Both binaries and their dedicated arithmetic, compact-storage and shared-sigma
   probes build from the same source manifest and pass with the exact expected
   coverage counts.
2. Small same-geometry checkpoints are byte-identical across uninterrupted,
   same-binary resume and cross-binary resume paths in both directions.  This
   establishes compatibility for the build-time scheduling change without writing
   or importing any campaign checkpoint.
3. Each arm completes a bounded DP-weight-42 correctness run at 385,024 workers,
   reference-replays at least the first 300 reports with zero mismatches and zero
   drops, and emits a nonempty corpus.  The repository's shared format-aware
   `benchmarks/dp_identity.py` must report identical sorted record multisets.  This
   is the specified fallback because the client has no GPU deterministic-spread
   replay entrypoint.
4. Exactly two full-work warmups, one per arm, are excluded.  Five full-work pairs
   then alternate order as AB, BA, AB, BA, AB.  Every row must report the exact
   scalar work above, the fixed feature markers and zero drops.
5. All five candidate/control pair ratios exceed 1.0 and a two-sided 95% paired
   confidence interval over log ratios has a lower bound above 1.0.  A failed
   gate suppresses timing; a failed timing row, mixed source/binary identity, or
   interval touching 1 is a no-go.  The absolute 26 B/s objective is reported
   independently.

No optional early screen changes the five-pair trial count.  Failures and partial
logs remain in the Modal result archive.  Full SASS, extracted cubins, ptxas logs,
runtime resources, source hashes, executable hashes, toolchain identity, GPU
identity and raw timing logs are retained so the derived result can be replayed.
Current main does not emit a runtime marker for `PACKED_INLINE_POLY`; the gate
therefore requires exactly one matching `-DECC_PACKED_INLINE_POLY=0` or `=3`
definition in the corresponding retained nvcc client command and rejects the
opposite definition.  Runtime resources, binary hashes and SASS bind every
sample to that executable.

## Measured result

Attempt 1 stopped before timing because the producer expected a runtime inline
marker that current main does not print.  Its arithmetic, storage, shared-sigma,
checkpoint and replay/corpus observations passed, but the attempt remains a
[`PRODUCER_FAILURE`](attempt-1-producer-failure.json) with no scientific
decision.  The additive repair changed only that producer gate to the exact
retained nvcc-definition check described above.

Attempt 2 ran on one RTX PRO 6000 (UUID
`GPU-85cc08d8-b26d-7c50-1b69-72ea66d6787a`, driver 580.95.05) from source
commit `9f2086cc791bdf1e42fe9f83e0e6d0b520815c37`, with CUDA 13.3.73 native
sm_120 code.  Both arms passed the arithmetic, compact-storage and shared-sigma
GPU gates.  Bidirectional cross-binary resume produced the same final
checkpoint as each uninterrupted run.  Each arm emitted 29,779 headerless v1
records, replayed 300/300 against the reference with zero mismatches and zero
drops, and the shared format-aware comparator produced the same sorted-record
SHA-256, `17e29b9695e175f306b0c76136b7331601596c98b4a9e1c27311f03c6b4a635d`.

The two excluded full-work warmups were 14.839194 B/s control and 15.134070
B/s inline3.  The five ranked pairs were:

| Pair | Order | Control B/s | inline3 B/s | inline3 / control |
|---:|:---:|---:|---:|---:|
| 1 | AB | 14.635580 | 15.007083 | 1.025383552 |
| 2 | BA | 14.552702 | 14.972113 | 1.028820146 |
| 3 | AB | 14.534256 | 14.941312 | 1.028006662 |
| 4 | BA | 14.517639 | 14.936047 | 1.028820664 |
| 5 | AB | 14.510774 | 14.931758 | 1.029011823 |

The control and inline3 medians are **14.534256 B/s** and **14.941312 B/s**,
a 2.800666% ratio of arm medians.  The median paired ratio is **1.028820146**;
the geometric mean paired ratio is 1.028007672 with a two-sided 95% log-ratio
interval of **[1.026123156, 1.029895649]**.  All five pairs favor inline3, so
the frozen admission rule says **promote `PACKED_INLINE_POLY=3` for this exact
checkpoint-compatible sigma preset**.

This does not meet the separate 26 B/s objective: 14.941312 B/s is 57.466585%
of the objective, 11.058688 B/s short.  It is bounded kernel-engineering
evidence on one allocation, not a full-ECDLP result or cryptanalytic advance.
The walk resources moved from 104 to 91 registers; both arms have zero walk
stack/local bytes, 1,792 function shared bytes plus the 1,024-byte driver
reservation, and a 2,816-byte compiled shared extent.  Full cubins and SASS are
retained in the raw archive and hash-bound in the
[`independent audit`](independent-audit.json).

The canonical derived artifact is [`result.json`](result.json).  The raw
68-file archive is stored in the `ecc2k130-jobs` Modal volume at
`594e18cedb0048d6aafba591eba51b11/results.tgz`: 12,187,638 bytes, SHA-256
`bba5b4bc46a22a6d297114fb96933e8b143526cace83b58908414525df7b0fd2`.
Fetch it with:

```sh
modal run modal_job.py \
  --fetch-token 594e18cedb0048d6aafba591eba51b11 \
  --out /tmp/ecc2k-sigma-inline-refresh-audit
```

## Reproduction

From `ecc2k130/` on the committed protocol revision:

```sh
modal run --detach modal_job.py \
  --job benchmarks/sigma-inline-refresh/gpujob.sh \
  --out /tmp/ecc2k-sigma-inline-refresh-20261001 \
  --gpu RTX-PRO-6000

bash benchmarks/sigma-inline-refresh/audit.sh \
  /tmp/ecc2k-sigma-inline-refresh-20261001/results
```

The job script invokes only bounded `--bench`, bounded `--verify 300` and small
checkpoint-continuation commands.  It never invokes `--load`, `--replay`, a search
entrypoint, `mergeCorpus` or `solveCorpus`.
