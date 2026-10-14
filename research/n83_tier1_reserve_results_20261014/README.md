# N83 K2940 fresh-point Tier 1 result sealer

`materialize.py` binds the five independently replayed same-point IC/rho
pairs to five canonical one-target workloads and measurement rows. It pins
the executed source commit `552e379748cd0befe082ca66b2560e7316769b61`,
the selected K2940 root22 candidate, the panel selection and exact reserve
fixture hashes. Each row keeps online timings, exclusive IC phases, setup,
memory, noise status and the corresponding strong-rho receipt.

After downloading the raw GitHub artifact, independently run
`research/n83_tier1_reserve_20261014/verify.py` on its `results` directory.
Compare the downloaded ZIP SHA-256 to GitHub's artifact digest, then run:

```sh
python3 research/n83_tier1_reserve_results_20261014/materialize.py \
  --artifact-zip /absolute/path/to/raw.zip \
  --artifact-sha256 <64-lowercase-hex> \
  --analysis-sha256 <64-lowercase-hex> \
  /absolute/path/to/extracted/results
```

The sealer verifies the raw digests, selection and candidate manifests,
source and target hashes, and all five pair identities. It writes
`workloads.json` and `measurements.json` once. Repeated runs accept identical
content and reject drift. Incomplete pairs retain their status and null
speedup fields.
