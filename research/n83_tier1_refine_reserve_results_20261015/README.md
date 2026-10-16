# N83 K3300 fresh Tier 1 result sealer

`materialize.py` binds the five frozen public targets, exact selected
K3300/root11 candidate, executed source and binaries, raw cloud artifact,
independent verifier analysis, and Tier 1 host receipt to five canonical
one-target workload and measurement rows. Its source commit is
`d0e4aa13074c9fea910592c6e26b60b5e97fc439`, the cloud execution
commit. Each row retains the same-point rho online time, IC phase sum,
target-independent IC preparation, process wall and RSS, noise status, and
failed or capped arm status.

After downloading the artifact ZIP, compare its SHA-256 with GitHub's
published artifact digest. Extract it, then run the independent verifier
against the extracted `results` directory and compare its `analysis.json`
SHA-256 with the cloud analysis. Seal the record with:

```sh
python3 research/n83_tier1_refine_reserve_results_20261015/materialize.py \
  --artifact-zip /absolute/path/to/raw.zip \
  --artifact-sha256 <64-lowercase-hex> \
  --analysis-sha256 <64-lowercase-hex> \
  /absolute/path/to/extracted/results
```

The sealer rejects changed source, selection, candidate, or target fixtures.
Existing output files accept identical content and reject changed content.
