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

```sh
modal run --detach modal_job.py \
  --job benchmarks/batch-hints/gpujob.sh \
  --out /tmp/ecc2k-batch-hints --gpu RTX-PRO-6000
python3 benchmarks/batch-hints/summarize.py \
  /tmp/ecc2k-batch-hints/results \
  --out /tmp/ecc2k-batch-hints/result.json
```
