# Sealed N83 compact4 cross-layout comparison

[RESULT.md](RESULT.md) gives the measured result and evidence. The source
protocol is [../n83_compact4_crosslayout_20261016/PROTOCOL.md](../n83_compact4_crosslayout_20261016/PROTOCOL.md).
`selection.json` retains the frozen three-way decision, `records_index.json`
pins the lossless compressed `workload.json.gz` and `measurements.json.gz`,
and the three candidate manifests identify the exact bases and algorithms.

To reproduce the versioned records from the archived raw artifact:

```sh
python3 research/n83_compact4_crosslayout_results_20261016/materialize.py \
  --artifact-zip /Volumes/SSD990-2/n83-archive-20261016/raw/n83-compact4-crosslayout-38848691554/raw.zip \
  --artifact-sha256 cff3b20efe5189d0ea952cbb0dc2f2e856c71fc7547ac4698ca582241b3f9c2e \
  --analysis-sha256 3e3079efb0b10759b6bc1f43af0ee52081eaf772e49d6b6fe104991ce5ac7921 \
  --artifact-id 12028240172 \
  /Volumes/SSD990-2/n83-archive-20261016/raw/n83-compact4-crosslayout-38848691554/raw/results
```

The [current-head replay workflow](../../.github/workflows/n83-compact4-crosslayout-results.yml)
downloads that artifact from GitHub, checks its SHA-256, regenerates the
public fixtures, independently replays the six bases and 48,000 target
records, and requires `materialize.py` to reproduce every committed record.
