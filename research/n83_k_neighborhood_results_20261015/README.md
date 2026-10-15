# N83 root11 K-neighborhood result sealing

The four-arm producer is frozen at commit
`85ad48ca4ab668e850623416b814a350b162beef`, exact IC source SHA-256
`a1eb3b937106bb7cfdee6c54d6422ceb126b628462361c9ff7a127127f230c89`.
It reuses the canonical 4,000-point workload ID `12e14f274f8c`. The
comparison is secondary exploratory data; the selected K3300 factor base
already has a separate fresh-target Tier 1 result.

The completed raw ZIP is GitHub artifact `11946288154`, SHA-256
`f4e1edb15db6e3a4306f5f208e8f0714ef1a1e81d738ed9b30d3d28d3c30649b`;
its independent analysis SHA-256 is
`fe9c87c6730befc67fa8338204e48edc35dd3c66aa1cbd21b55fa29036b8eb88`.
Extract it outside the repository, run the frozen
`research/n83_k_neighborhood_20261015/verify.py` against its `results/`
directory, and compare the local `analysis.json` SHA-256 with the archived
cloud value. Then seal the artifact:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 \
  research/n83_k_neighborhood_results_20261015/materialize.py \
  --artifact-zip /absolute/path/raw.zip \
  --artifact-sha256 f4e1edb15db6e3a4306f5f208e8f0714ef1a1e81d738ed9b30d3d28d3c30649b \
  --analysis-sha256 fe9c87c6730befc67fa8338204e48edc35dd3c66aa1cbd21b55fa29036b8eb88 \
  /absolute/path/extracted/results
```

The sealer checks the source commit, executable and runner hashes, public
point and answer sidecar digests, arm set, independent analysis and raw ZIP.
It writes a canonical candidate manifest for each completed arm, one
measurement row per arm including failures and caps, a lossless compressed
workload and measurement file, a selection record, and a digest index.
Existing records are compared byte for byte and never overwritten with a
different outcome. Preserve the raw ZIP and local replay log separately.
