# N83 1,000-point refinement result sealer

`materialize.py` binds the independently replayed refinement artifact to the
canonical N83 curve, exact factor-base point sets, candidate IDs, frozen
workload ID and one measurement row per arm. It is pinned to the experiment
source commit `b1a4131558b3121ed562ca830abce972aa1dc580` and the 1,000
public-point and answer-sidecar SHA-256 values in the script. Existing
candidate identities are reused when the complete base and algorithm are
identical; the executed experiment commit remains in each measurement row.

After downloading the raw artifact and independently running
`research/n83_factor_base_refine_20261014/verify.py` against its `results`
directory, compare the downloaded ZIP digest to GitHub's artifact digest.
Then run:

```sh
python3 research/n83_factor_base_refine_results_20261014/materialize.py \
  --artifact-zip /absolute/path/to/raw.zip \
  --artifact-sha256 <64-lowercase-hex> \
  --analysis-sha256 <64-lowercase-hex> \
  /absolute/path/to/extracted/results
```

The sealer checks the two supplied digests, executed source and fixture
hashes, all frozen arms and selection status. It writes `workload.json`,
candidate manifests for verified arms, `measurements.json`, and
`selected_candidate.json` once; repeated runs accept identical content and
reject drift. A failed or over-cap arm remains a row and prevents a complete
selection receipt.
