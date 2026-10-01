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
kernel-engineering evidence; it does not establish the 26 B/s objective. The
harness is prepared only and must not be launched without separate admission:

```sh
modal run --detach modal_job.py \
  --job benchmarks/batch-hints-fast2/gpujob.sh \
  --out /tmp/ecc2k-batch-hints-fast2 --gpu RTX-PRO-6000
python3 benchmarks/batch-hints-fast2/summarize.py \
  /tmp/ecc2k-batch-hints-fast2/results \
  --out /tmp/ecc2k-batch-hints-fast2/result.json
```
