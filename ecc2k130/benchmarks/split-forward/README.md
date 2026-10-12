# Split-forward v3 table-walk screen

This experiment compares one source tree and one v3 table iteration map on one
RTX PRO 6000. The control is `TABLE_SPLIT_FORWARD=0`; the candidate is
`TABLE_SPLIT_FORWARD=1`. The candidate runs reporting and `twSelect` for all
slots before any prefix product is live, then reloads the unchanged polynomial
coordinates, takes the selected tag from `hist & 0xffff`, and builds the same
prefix chain. It changes scheduling and state traffic, not the step function.

The compile-only CUDA 13.3.73 screen that admitted this experiment measured:

| build | registers | stack | spill stores | spill loads |
|---|---:|---:|---:|---:|
| current v3 control | 128 | 448 B | 36 B | 84 B |
| split-forward candidate | 128 | 400 B | 0 B | 0 B |

These resource figures are not throughput evidence. The GPU protocol is frozen
before timing:

1. Build both arms from the same source and retain compiler logs and hashes.
2. For each arm, re-walk 300 reports with zero mismatches/drops and write a
   forced-common-work v3 corpus. Parse the 16-byte `ECC2KDT3` header, require
   version 3 and 32-byte records, then compare sorted record hashes.
3. Run two excluded warmups, control then candidate.
4. Run one bracketed screen: control, candidate, control. The candidate
   qualifies only when its completed-update rate is at least 1.005 times the
   faster control bracket.
5. Stop an unqualified candidate. A qualified candidate receives exactly three
   paired confirmations, alternating order by pair. No result here changes the
   rho operation count or proves the 26 B/s objective.

Run through an existing launcher, explicitly requesting one RTX PRO 6000:

```sh
modal run --detach modal_job.py \
  --job benchmarks/split-forward/gpujob.sh \
  --out /tmp/ecc2k-split-forward --gpu RTX-PRO-6000
python3 benchmarks/split-forward/summarize.py \
  /tmp/ecc2k-split-forward/results \
  --out /tmp/ecc2k-split-forward/result.json
```

`gpujob.sh` invokes the summarizer before deciding whether confirmation is
admitted. The summarizer treats an invalid header, partial record, failed
replay, missing identity marker, malformed rate, or incomplete confirmation as
invalid evidence rather than a zero or a winner.
