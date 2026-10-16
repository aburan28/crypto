# Reproduce the N83 root21/root22 challenger record

The raw ZIP is GitHub Actions artifact `11997698000` from run
`38770004923`. Its SHA-256 is pinned in `materialize.py`. Keep it unchanged,
extract it to a separate directory, and run the executed source verifier on
`results/`. The independently reproduced `analysis.json` must have SHA-256
`7d10662fe9987fd1766734188855a57f355b4244f26a815c5e065d7a93f85551`.

The executed source commit is
`d6293e278a037ef34a6de0101844fcfbf97703e2`. Regenerate its ignored
`research/n83_legacy_challenger_20261015/targets.jsonl` and
`answers.jsonl` with that commit's `generate.py`. The materializer checks
their hashes against the cloud fixture receipt and source manifest.

```sh
python3 research/n83_legacy_challenger_results_20261015/materialize.py \
  --artifact-zip /absolute/path/raw.zip \
  /absolute/path/extracted/results
```

The materializer accepts an existing sealed record only when it is
byte-identical. `RESULT.md` records the comparison and resource envelope;
`records_index.json` binds the compressed measurements and workload.
