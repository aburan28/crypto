# N83 root11 Bloom result records

The completed two-arm cloud run is
[38674728464](https://github.com/aburan28/crypto/actions/runs/38674728464),
executed at commit `a9b44ee1454beb2bbe8614b681714df07fd11f3c`.
Artifact `11945754917` has SHA-256
`098c1ea1383c3ddaf78aaf1ac28bc32c80cb2caa5c69ccfb155bb060f5d2bda2`.
Its independent `analysis.json` has SHA-256
`2976f932e8eb2fd2d285b1f40f065c43dd6573fb6879845bd915f9b7ee675a6a`.

Extract the ZIP outside the repository and run the frozen
`research/n83_root11_bloom_20261015/verify.py` on the `results/` directory.
Compare the resulting analysis digest to the value above. Seal or reproduce
the records with:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 \
  research/n83_root11_bloom_results_20261015/materialize.py \
  --artifact-zip /absolute/path/raw.zip \
  /absolute/path/extracted/results
```

The sealer checks the raw artifact, analysis, source and public-point hashes,
then writes exact Bloom-off and Bloom-on IC1 candidate manifests, lossless
measurement rows, a selection record, and a digest index. The 1,000-point
canonical workload is shared with the earlier factor-base refinement panel.
The current-head workflow independently replays the frozen source and checks
the committed result records byte for byte.
