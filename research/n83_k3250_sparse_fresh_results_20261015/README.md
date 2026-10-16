# Fresh N83 sparse factor-base confirmation

The sealed study compares the K3250 and K3300 factor bases with the exact
root11 sparse canonical selector. Its 20,000 public target points are disjoint
from earlier panels. The exact target and answer files are retained in
executed source commit `402c3142fb3f1cfe70b9de58d54e51e50b2a40c5`
and recreated byte-for-byte by
`research/n83_k3250_sparse_fresh_20261015/generate.py` for new runs.
The answer sidecar is used for independent replay only. Both arms use the
same target order, 14 GiB memory guard, and 1,800-second child cap. The
decision rule promotes K3250 only when the upper endpoint of the paired,
contiguous 50-target block bootstrap interval is below zero.

`materialize.py` consumes an archived producer ZIP and its extracted
`results` directory. It verifies the ZIP, analysis, source, fixture, base,
and execution identities before writing immutable candidate, workload,
measurement, and selection records. It preserves incomplete arms and their
resource receipts. The reusable factor-base and log preparation is recorded
separately from per-target online time.

```sh
python3 research/n83_k3250_sparse_fresh_results_20261015/materialize.py \
  --artifact-zip /absolute/path/to/producer.zip \
  --artifact-sha256 SHA256_OF_ZIP \
  --analysis-sha256 SHA256_OF_EXTRACTED_RESULTS_ANALYSIS_JSON \
  /absolute/path/to/extracted/results
```
