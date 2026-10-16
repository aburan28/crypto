# Reproduce the N83 direct comparison record

The raw ZIP is GitHub Actions artifact `11997750150` from run `38768544528`.
Its SHA-256 is pinned in `materialize.py`. Keep the ZIP unchanged, extract it
to a separate directory, and run the source verifier against `results/`.
The verifier writes `analysis.json`; its SHA-256 must match
`59e2a1cd6bf404ed46a4e82694d5ecde61f191f5e97142ad27282948cae568a7`.
The executed source commit is
`ae6a10b641e17629bd282d8205787ede129dbf98`. Regenerate the ignored
`research/n83_best_e2e_20261015/targets.jsonl` and `answers.jsonl` with
that commit's `generate.py`; the cloud fixture receipt and all four digests
are checked by the materializer.

From a checkout containing this result package and the regenerated fixture:

```sh
python3 research/n83_best_e2e_results_20261015/materialize.py \
  --artifact-zip /absolute/path/raw.zip \
  /absolute/path/extracted/results
```

The script writes each record only when absent or byte-identical to the
sealed version. `RESULT.md` gives the measured selection and its resource
conditions; `records_index.json` binds the compressed measurement and
workload data.
