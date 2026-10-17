# N83 sparse plus Bloom Tier 1 result sealing

The five alternating same-point IC/rho pairs use the combined K3300/root11
candidate frozen in `research/n83_tier1_combined_20261015/selection.json`.
The candidate's factor base has 547,800 subgroup points and 3,300 folded
columns. Its exact root index has 903,852,075 keys. Both sparse canonical
selection and the three-hash Bloom prefilter are enabled.

The producer source is frozen at commit
`92e6f373312435bf89439b8a08bf49aa383989b2`. Once its raw artifact is
archived, rerun the frozen verifier against the extracted `results/`
directory, compare the analysis digest, then seal five separate one-target
workloads and measurement rows:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 \
  research/n83_tier1_combined_results_20261015/materialize.py \
  --artifact-zip /absolute/path/raw.zip \
  --artifact-sha256 RAW_ZIP_SHA256 \
  --analysis-sha256 ANALYSIS_JSON_SHA256 \
  /absolute/path/extracted/results
```

The sealer checks source, binary, candidate, selection, target, and raw
receipt identities. It preserves every pair status, resource cap, noise
measurement, preparation cost, and independently replayed scalar result.
Only pairs satisfying the Tier 1 correctness and noise gate receive online
IC/rho times and ratios.
