# Reconverged v3 cycle-hint screen

This experiment compares two split-forward builds of the same v3 table walk
on one RTX PRO 6000. Both builds select every slot before the prefix product
is live. The control immediately resolves each hinted slot. The candidate
first pushes every raw tag into that slot's history and records a per-lane
slot mask. A warp-reconverged loop then resolves one pending slot per lane at
a time and replaces only the newest history tag before the addend pass.

`twCycleTag` is independent of incoming history. Slots are independent walks,
and no point or history advances between raw selection and resolution. The
candidate therefore changes SIMT scheduling only: it uses the same hint,
bounded exact point probe, anchor rule, selected tag and final history.

The frozen admission protocol is:

1. Build `TABLE_SPLIT_FORWARD=1` with `TABLE_BATCH_HINTS=0/1` from one source.
2. Re-walk 300 reports for each build with zero mismatches or drops. Parse the
   v3 header and require identical sorted 32-byte record multisets.
3. Run exactly two excluded warmups and a control-candidate-control screen.
4. Admit three alternating confirmation pairs only if the candidate reaches
   at least 1.005 times the faster screening control.

The timing uses the existing complete-update counter. Passing this screen is
kernel-engineering evidence; it does not establish the 26 B/s objective.

## Measured result

The candidate passed the correctness gate on one RTX PRO 6000 with CUDA
13.3.73. Both arms replayed 300/300 reports with zero drops and produced the
same sorted multiset of 1,480,278 v3 records, SHA-256
`2fc6f84d4b1969240418cf6f8e779df9538b1d8f01eb6ff4c6a00be7324fae01`.
Both kernels used 128 registers, a 400-byte frame and zero reported spills.

The 32-launch screen measured 1.473520 / 2.931030 / 1.473470 B/s. Three
alternating 64-launch confirmation pairs measured:

| Pair | Split control B/s | Reconverged B/s | Ratio |
|---|---:|---:|---:|
| 1 | 1.052700 | 2.448134 | 2.325576 |
| 2 | 1.052232 | 2.449438 | 2.327850 |
| 3 | 1.050904 | 2.449169 | 2.330535 |

The candidate median is **2.449169 B/s**, versus **1.052232 B/s** for the
control; the median paired ratio is **2.327850**. Every confirmation row
completed 100,931,731,456 scalar updates. This repairs a large v3 SIMT
divergence cost, but it remains below the correctness-current sigma walk and
far below 26 B/s.

[`result.json`](result.json) contains the frozen summary. The raw archive is
61,385,506 bytes with SHA-256
`c143b3af1f48aee44706af84b4bfcd8c099156afcaac3c585f0cf56819095f7b`.
[`independent-audit.json`](independent-audit.json) independently checks the
corpora, replay, exact work, raw-log hashes, timing, resources, source hashes
and schedule equivalence. The binaries were not retained, so their recorded
hashes are provenance rather than independently recomputed identities.

```sh
modal run --detach modal_job.py \
  --job benchmarks/batch-hints/gpujob.sh \
  --out /tmp/ecc2k-batch-hints --gpu RTX-PRO-6000
python3 benchmarks/batch-hints/summarize.py \
  /tmp/ecc2k-batch-hints/results \
  --out /tmp/ecc2k-batch-hints/result.json
```
