# Stricter-reader compatibility check

The frozen executable at source
`5521605f86ae444a9f4a433ca9b2f80d541e8470` read the original corrected BMSS
baseline and final rung again. It accepted their mandatory source/binary and
output identities, known count model, full matrix and terminal PASS receipts.
Its fresh comparison passed all 96 exact coefficient-output matches. This
rechecks historical evidence; it adds no performance samples or new speedup.

```sh
isogeny-perf compare evidence/baseline/counts-valid.jsonl \
  evidence/rung3/counts.jsonl parser-rechecks-20261015/BMSS_COMPARISON.json
```

The original raw files remain inside `FROZEN_LOCAL_RUNS.tar.gz` and unchanged
locally. Their digests are in the fresh comparison. The initial interrupted
`counts-complete.jsonl` has only its header and is excluded. The stronger reader
retains the existing fail-closed treatment of partial files.
