# N83 fresh confirmation result package

The source protocol and exact two-candidate admission are in
`../n83_edge_confirmation_20261015`. GitHub run 38809654754 executes four
mirrored arms on 4,000 fresh public subgroup points. Its source commit is
`a369212178ffd28ec30e0e20eaa4729ec7b81e8d`.

After the run, archive the raw GitHub artifact and compare the ZIP against
GitHub's SHA-256 digest. Extract it outside this checkout. Replay the exact
base dumps, four-point relations, recovered scalars, resource counters,
calibration and timing phases:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 \
  research/n83_edge_confirmation_20261015/verify.py \
  /absolute/path/to/extracted/results
```

Compare the local `analysis.json` SHA-256 with the cloud analysis SHA-256.
Then run `materialize.py` with the ZIP and analysis digests and artifact ID.
The script freezes all four run rows, the 4,000-point workload, the full
paired comparison, and an index of raw receipt hashes. It refuses changed
bytes for any previously sealed file. Only a `SELECTED` result can populate
`../n83_final_tier1_20261015/selection.json` through its `select.py`.
