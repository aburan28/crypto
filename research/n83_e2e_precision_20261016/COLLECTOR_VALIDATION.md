# N83 provider-artifact collector replay

The provider-artifact collector recovered the complete five-point holdout
artifact from GitHub Actions run `38955665283` and independently reproduced
its analysis. This checks the collection path used for the pending precision
panel against a completed 245-file result set. It does not add a new timing
observation to the factor-base comparison.

| Check | Recorded value |
| --- | --- |
| Executed head | `e44c5022519e190cfa72fad326e1ba399d49893d` |
| Artifact name | `n83-e2e-holdout-38955665283-1` |
| Provider and downloaded ZIP SHA-256 | `d2e73fff94e3dc972c47be646667eef7ab2fa4a20f661f44f02323b25d8ff060` |
| Extracted result files | 245 |
| Archived `analysis.json` SHA-256 | `b436fcad4ad5ae8af8b773c1dcbfb8c043fda12433fa509ae49639cb77f33b89` |
| Independent verifier exit | `0`, empty stderr |
| Independent verifier stdout SHA-256 | `b3efa6240fe6188e4d4f746bcd6b596469b837ac2c935d16a10ec5e31f2aa2b3` |
| Analysis comparison | Parsed JSON equal; `BOUNDED_UNKNOWN`, five points, `PASS_INDEPENDENT_REPLAY` |

The exact raw ZIP is retained at
`/Volumes/SSD990-2/n83-archive-20261016/raw/n83-e2e-holdout-38955665283.zip`.
The collector rehearsal used this command from the repository root:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 research/n83_e2e_precision_20261016/collect.py \
  --run-id 38955665283 \
  --head-sha e44c5022519e190cfa72fad326e1ba399d49893d \
  --artifact-name n83-e2e-holdout-38955665283-1 \
  --destination /private/tmp/n83-collector-rehearsal-38955665283
```

The independent replay used the holdout verifier on the extracted `results`
directory. Its compact stdout and the archived formatted analysis have
different byte hashes; parsing both as JSON gives identical objects. The
precision collector will use the same provider digest, member extraction,
head binding and local replay checks on its own finished artifact.
