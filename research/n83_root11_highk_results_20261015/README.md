# Sealed N83 root11 high-K screen

`RESULT.md` records the balanced eight-arm K3300–K3450 root11/wide3
measurement and its frozen challenger choice. The raw GitHub Actions
artifact is `12000855288` from run `38773883285`; retain its ZIP and
extracted `results` directory outside the Git tree. The artifact and
analysis digests are pinned in `materialize.py` and `selection.json`.

To reproduce the versioned candidate manifests, eight canonical run rows,
workload ID, and selection from the archived evidence:

```sh
python3 research/n83_root11_highk_results_20261015/materialize.py \
  --artifact-zip /absolute/path/to/n83-root11-highk-38773883285.zip \
  /absolute/path/to/n83-root11-highk-38773883285/results
```

The script rejects any changed output, source or base. Reproduce the
independent mathematical replay first with
`research/n83_root11_highk_20261015/verify.py` on the extracted `results`
directory. The current-head GitHub workflow repeats both steps against the
archived artifact.
