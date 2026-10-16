# Reproduce the N83 root4 Bloom result

The raw ZIP is GitHub Actions artifact `11999315231` from run
`38772262240`. Its SHA-256 is pinned in `materialize.py`. Keep it unchanged,
extract it to a separate directory, and run the executed source verifier on
`results/`. The expected `analysis.json` SHA-256 is
`15b47de167f99e004bccada699f739ff7d0d7e7f7e7960239d63b4b00183e547`.

The executed source commit is
`afea7179c40c9e00f8701a32c596e47cd8129f62`. Regenerate its ignored
`research/n83_root4_bloom_20261015/targets.jsonl` and `answers.jsonl`
with that commit's `generate.py`. The materializer checks their hashes
against the cloud fixture receipt and source manifest.

```sh
python3 research/n83_root4_bloom_results_20261015/materialize.py \
  --artifact-zip /absolute/path/raw.zip \
  /absolute/path/extracted/results
```

The materializer accepts an existing sealed record only when it is
byte-identical. `RESULT.md` records the paired comparison and resource
envelope; `records_index.json` binds the compressed measurements and workload.
