# Fresh Bloom-on N83 Tier 1 result sealing

The five-pair producer is frozen at commit
`b70a2611fb8bc130c708ac91f26981592b2c1f78`, IC source SHA-256
`814afa03ae1c386e62e518be3c7a55fb108b09c41b05bef64ca1371787466c39`.
The Bloom-on candidate ID is
`IC1N83Ckb1a1fb547800PDP4rootRCguidedLAgaussTDdirectISO0hc3b6db4ff8c0`.
The five public points have SHA-256
`ef616c5254a9d883e962024b9f4bd7f794cfe5576c73a93b099ac5be6fe45668`.

After the cloud workflow uploads its raw ZIP, retain the artifact ID and
SHA-256. Extract it outside the repository, run the frozen
`research/n83_tier1_bloom_reserve_20261015/verify.py` on its `results/`
directory, and compare the resulting `analysis.json` digest to the cloud
value. Then seal the five one-target workload and measurement rows:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 \
  research/n83_tier1_bloom_reserve_results_20261015/materialize.py \
  --artifact-zip /absolute/path/raw.zip \
  --artifact-sha256 <raw-zip-sha256> \
  --analysis-sha256 <independent-analysis-sha256> \
  /absolute/path/extracted/results
```

The sealer checks the executed source, selected manifest, public points,
answer sidecar, artifact and independent analysis. It preserves each pair's
verification and noise status, exact point, separate target online times,
exclusive IC phases, reusable setup and memory. A failed or noisy pair stays
in its row and does not become a Tier 1 speedup.
