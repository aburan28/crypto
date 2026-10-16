# N83 unresolved-finalists same-point Tier 1 result package

The source protocol in `../n83_unresolved_tier1_20261015` fixes one new
public point and both complete finalists before timing. GitHub run
38818664674 executes two alternating IC/rho pairs on the same cloud
allocation. Its exact source commit is
`e1a77618587533364cf2236f09cc41b5bc816429`.

After the run, archive its raw GitHub ZIP and compare the SHA-256 with the
artifact digest. Independently replay both factor bases, four-point
relations, recovered IC and rho scalars, calibration, cgroup and system
counters, and all exclusive IC online phases:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 \
  research/n83_unresolved_tier1_20261015/verify.py \
  /absolute/path/to/extracted/results
```

Compare the local `analysis.json` SHA-256 with the cloud analysis SHA-256.
`materialize.py` then seals one workload record, two paired run records,
the summary and raw receipt hashes. It preserves failed or noisy pairs and
refuses any changed sealed content. The factor-base selection remains
`BOUNDED_UNKNOWN` in every output from this diagnostic.
