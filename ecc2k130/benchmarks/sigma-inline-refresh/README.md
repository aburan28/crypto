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
