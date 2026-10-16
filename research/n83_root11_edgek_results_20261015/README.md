# Sealed N83 root11 K3450/K3500/K3550 result

`RESULT.md` reports the paired 1,000-point edge screen. `selection.json`
retains the complete frozen comparison and its bounded decision.
`candidate_K*_root11_wide3.json` names the exact bases and algorithms;
`measurements.json.gz` stores all six run rows and per-target phase records;
`workload.json.gz` freezes the public points and input law. The small
`records_index.json` pins both compressed files and the raw artifact.

To reproduce these records from the archived raw ZIP and independently
replayed extracted directory, run:

```sh
python3 research/n83_root11_edgek_results_20261015/materialize.py \
  --artifact-zip /absolute/path/to/n83-root11-edgek-38794757719.zip \
  /absolute/path/to/n83-root11-edgek-38794757719/results
```

The materializer verifies the ZIP, analysis, source, public fixture,
allocation and arm identities before writing. It refuses to overwrite any
previously sealed record with different contents.
