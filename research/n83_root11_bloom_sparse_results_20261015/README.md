# N83 combined sparse/Bloom result sealing

The two-arm producer is frozen at commit
`9289963251767f88340e1c8900d18f50f1b40ee0`, IC source SHA-256
`814afa03ae1c386e62e518be3c7a55fb108b09c41b05bef64ca1371787466c39`.
The arms share the canonical 1,000-point workload ID `0d7b05d9fca6`.

After the cloud workflow uploads its raw ZIP, retain its artifact ID and
SHA-256. Extract the ZIP outside the repository, run the frozen
`research/n83_root11_bloom_sparse_20261015/verify.py` against its `results/`
directory, and compare the local `analysis.json` digest with the cloud value.
Then seal the artifact:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 \
  research/n83_root11_bloom_sparse_results_20261015/materialize.py \
  --artifact-zip /absolute/path/raw.zip \
  --artifact-sha256 <raw-zip-sha256> \
  --analysis-sha256 <independent-analysis-sha256> \
  /absolute/path/extracted/results
```

The sealer checks source and input identity, raw hashes, both arm records,
and the frozen paired statistic. It emits exact sparse/Bloom-off and
sparse/Bloom-on candidate manifests, lossless compressed measurement rows,
a selection record and digest index. It preserves a complete result whether
Bloom gives a confirmed additional gain or the paired interval spans zero.
Incomplete arms remain rows with their recorded status and resource caps.
