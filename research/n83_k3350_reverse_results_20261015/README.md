# N83 K3350/K3300 reversed-order replication records

This directory seals the source-pinned reversed-order comparison defined in
[`PROTOCOL.md`](../n83_k3350_reverse_20261015/PROTOCOL.md). Its primary
records are `selection.json`, one candidate manifest per completed arm,
`workload.json.gz`, `measurements.json.gz`, and `records_index.json`.

The source commit is `53339d82a372e09958dd43bb7391e283de8933f5`.
The generated point and known-scalar answer files have SHA-256 digests
`fc8be4c2c946bdbcb78a6845c8e2be3d1a51728c2ed4d812f6546ea5b000055e`
and `7121dcaf478061b87f54f28a0dc8e5af28e2637fd34a6837dd8f6434b6e29e8b`.
The answer sidecar is used only for independent replay, outside the measured
producer inputs.

After the raw Actions artifact is downloaded and its ZIP and `analysis.json`
digests are fixed, regenerate both fixtures from the source-pinned scripts,
run `verify.py` against the unpacked `results/` directory, then seal records:

```sh
python3 research/n83_k3250_sparse_fresh_20261015/generate.py
python3 research/n83_k3350_sparse_fresh_20261015/generate.py
python3 research/n83_k3350_reverse_20261015/verify.py /path/to/results
python3 research/n83_k3350_reverse_results_20261015/materialize.py \
  --artifact-zip /path/to/raw.zip --artifact-sha256 ZIP_SHA256 \
  --analysis-sha256 ANALYSIS_SHA256 /path/to/results
```

The sealer refuses changed source, fixture, raw-artifact, and analysis digests;
it writes each canonical record once and verifies byte equality on replay.
The source-pinned CI replay retains the exact commands and final status.
