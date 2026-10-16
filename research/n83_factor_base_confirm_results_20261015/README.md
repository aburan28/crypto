# N83 fresh factor-base confirmation sealer

`materialize.py` binds the frozen 4,000-point confirmation artifact to the
exact N83 curve, executed source commit
`e65f35916dfe4791a60b938b8ab0bdc875d0ae7b`, both complete factor
bases, two canonical IC1 candidate IDs, the predeclared K3100-minus-K3300
contrast, and one measurement row per arm. It keeps every target interval,
exclusive online phase, setup, rank, memory and process status. The output
`confirmation.json` distinguishes a replayed result from a positive
confirmatory interval. The complete workload and per-target measurement
records are losslessly stored in deterministic `workload.json.gz` and
`measurements.json.gz`. `records_index.json` records the SHA-256 of both
archives and their decoded JSON bytes; use `gzip -dc` to read them.

Compare the downloaded artifact ZIP SHA-256 to GitHub's published artifact
digest and independently rerun
`research/n83_factor_base_confirm_20261015/verify.py` on its extracted
`results` directory. Compare the resulting `analysis.json` digest to the
cloud analysis, then run:

```sh
python3 research/n83_factor_base_confirm_results_20261015/materialize.py \
  --artifact-zip /absolute/path/to/raw.zip \
  --artifact-sha256 <64-lowercase-hex> \
  --analysis-sha256 <64-lowercase-hex> \
  /absolute/path/to/extracted/results
```

The sealer checks source, fixture, analysis, artifact and base digests.
Existing output files accept identical content and reject changed content.
