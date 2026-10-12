# Reconverged hints plus exact raw-two-cycle screen

This bounded RTX PRO 6000 experiment holds the measured B16 reconverged-hints
kernel fixed and changes only the exact raw-two-cycle shortcut. Both arms use
`TABLE_SPLIT_FORWARD=1 TABLE_BATCH_HINTS=1 BATCH=16 THREADS=512 MINBLOCKS=1`.
The control has `CYCLE_FAST2=0`; the candidate has `CYCLE_FAST2=1`.

After the first raw step `R1 = R + A_t`, the candidate short-circuits only
when the next raw tag is `t ^ ECC_TAG_EPS` and the second affine denominator
is nonzero. The two table addends are exact inverses, so associativity proves
`R1 + A_(t^eps) = R`. It applies the same distinguished-point stop, both
length-two cyclic eligibility tests, and the same strict `(weight, canonical
x)` anchor ordering as the general v3 probe. All other hints continue through
the unchanged bounded probe, reusing the already-computed raw tag at `R1`.

The frozen admission protocol is:

1. Build both B16 arms from one source tree and record ptxas resources.
2. Use the same run id, geometry, DP threshold, launch count, and capacity.
   Reference-replay the first 300 reports for each arm with zero mismatches or
   drops, parse the table-v3 header, and require identical sorted 32-byte
   record multisets. A failed or different corpus suppresses all timing.
3. Run exactly two excluded warmups and a control-candidate-control screen.
4. Admit three alternating confirmation pairs only if fast2 reaches at least
   `1.005` times the faster screening control.

The timing counter covers complete table-walk updates. Passing is bounded
kernel-engineering evidence; it does not establish the 26 B/s objective.

## Measured result

Both arms replayed 300/300 reports with zero drops and produced identical
sorted multisets of 1,480,278 v3 records, SHA-256
`2fc6f84d4b1969240418cf6f8e779df9538b1d8f01eb6ff4c6a00be7324fae01`.
Both kernels used 128 registers, a 400-byte frame and zero reported spills.
The screen measured 2.931950 / 4.190032 / 2.931331 B/s. Three alternating
64-launch confirmation pairs were:

| Pair | Control B/s | Fast2 B/s | Ratio |
|---|---:|---:|---:|
| 1 | 2.450018 | 3.593206 | 1.466604 |
| 2 | 2.450247 | 3.592601 | 1.466220 |
| 3 | 2.449482 | 3.592794 | 1.466757 |

The fast2 median is **3.592794 B/s** and the median paired ratio is
**1.466604**. [`result.json`](result.json) retains the full summary. The raw
archive is 61,389,365 bytes with SHA-256
`6662cd707f4a0add15ba47fa240c1dde83c3f8ace7db40637cf7fef51eaf2755`.
[`independent-audit.json`](independent-audit.json) recomputes the corpus,
replay, exact work, raw-log hashes, source manifests, resources and timing.
The post-run checker was strengthened to reopen every timed log and validate
those same hashes, work counters, feature markers, final rates and zero drops.
The executables were not retained, so their recorded hashes remain provenance
rather than independently recomputed identities.

Reproduce or audit with:

```sh
modal run --detach modal_job.py \
  --job benchmarks/batch-hints-fast2/gpujob.sh \
  --out /tmp/ecc2k-batch-hints-fast2 --gpu RTX-PRO-6000
python3 benchmarks/batch-hints-fast2/summarize.py \
  /tmp/ecc2k-batch-hints-fast2/results \
  --out /tmp/ecc2k-batch-hints-fast2/result.json
```
