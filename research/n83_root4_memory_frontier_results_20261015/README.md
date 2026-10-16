# N83 root4 memory-frontier result package

The sealed six-arm measurement and its factor-base decision are in [RESULT.md](RESULT.md). K3596 has the lowest descriptive online mean; the K3500-versus-K3596 selector remains `BOUNDED_UNKNOWN` at the frozen 99% rule.

The source protocol in `../n83_root4_memory_frontier_20261015` freezes three exact root4 standard-Bloom configurations, six mirrored arms, one fresh 1,000-point panel, and a 15.25 GiB common RSS guard before execution. GitHub run `38825895649` executes source commit `9c4a60cf5b49ad926e7eb391b2b80bd016893245`.

After the cloud run, archive its raw ZIP and compare SHA-256 with the GitHub artifact digest. Replay every base, public-target scalar, four-point relation, resource receipt and mirrored comparison with the executed source verifier. The replayed `analysis.json` SHA-256 must match the cloud file. Then run `materialize.py` with the archived ZIP, analysis SHA-256, and artifact ID. It writes canonical `IC1` manifests from executed base hashes, a compressed 1,000-point workload, six per-arm measurement records, a selection record, and a receipt index. Existing sealed records are accepted only when byte-identical.

```sh
PYTHONDONTWRITEBYTECODE=1 python3 research/n83_root4_memory_frontier_20261015/verify.py /absolute/path/extracted/results
PYTHONDONTWRITEBYTECODE=1 python3 research/n83_root4_memory_frontier_results_20261015/materialize.py \
  --artifact-zip /absolute/path/raw.zip \
  --artifact-sha256 EXPECTED_GITHUB_ZIP_SHA256 \
  --analysis-sha256 REPLAYED_ANALYSIS_SHA256 \
  --artifact-id GITHUB_ARTIFACT_ID \
  /absolute/path/extracted/results
```

An incomplete or noisy arm remains in the six-row record with its observed status and resource receipt; it cannot select a pipeline. The frozen paired rule is applied only when every arm is independently replayed and eligible.
