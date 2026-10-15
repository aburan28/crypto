# N83 sparse canonical selector result records

The completed two-arm cloud run is
[38676235613](https://github.com/aburan28/crypto/actions/runs/38676235613),
executed at commit `39f8c5e74b2916333c6be60f2e653d7651eb7308`.
Raw artifact `11947421783` has SHA-256
`be12dbe701fbe33092b2bb104d2e4cd07d321374596b1f8c69b83c332ac137fc`.
Its independent `analysis.json` has SHA-256
`6d555f4d6f8596fda5451b18e9f6cc0451fc63b2f792fd7baf1f33cf1fcf501c`.

Extract the ZIP outside the repository and run the frozen
`research/n83_sparse_canonical_20261015/verify.py` on its `results/`
directory. Compare the analysis digest to the value above. Seal or reproduce
the records with:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 \
  research/n83_sparse_canonical_results_20261015/materialize.py \
  --artifact-zip /absolute/path/raw.zip \
  /absolute/path/extracted/results
```

The sealer checks the raw artifact, analysis, source and public-point hashes,
then writes exact scan and sparse IC1 candidate manifests, lossless
measurement rows, a selection record and a digest index. The canonical
1,000-point workload is shared with the factor-base refinement panel. The
current-head workflow replays the frozen source and checks the committed
result records byte for byte.
