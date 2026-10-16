# N83 root4 high-K confirmation result package

The four-arm confirmation uses 8,000 fresh public points per arm and compares
K3500/root4/standard Bloom with K3596/root4/standard Bloom under one 15.25 GiB
observed RSS guard. The source protocol is
[`../n83_root4_highk_confirm_20261016/PROTOCOL.md`](../n83_root4_highk_confirm_20261016/PROTOCOL.md).
GitHub Actions run `38835419672` executed the frozen source commit
`adf9ca9b254ac7de8bb777446e74788d25e369a8`. Artifact `12021760871`
has ZIP SHA-256
`4a752ba2910ef52773de7128d933df9c7e8f7aad36c6eda3a561f921041d0a78`.

The source verifier reconstructs every factor-base orbit, four-point relation,
and scalar answer. It retains every process status, target phase, root-probe
count, RSS and cgroup receipt. `materialize.py` checks the frozen source,
fixture and artifact identities and writes candidate manifests, a compressed
8,000-point workload, four complete measurement rows, the frozen selection,
and a receipt index. Existing records are accepted only when byte-identical.

```sh
PYTHONDONTWRITEBYTECODE=1 python3 research/n83_root4_highk_confirm_20261016/verify.py /absolute/path/extracted/results
PYTHONDONTWRITEBYTECODE=1 python3 research/n83_root4_highk_confirm_results_20261016/materialize.py \
  --artifact-zip /absolute/path/raw.zip \
  --artifact-sha256 4a752ba2910ef52773de7128d933df9c7e8f7aad36c6eda3a561f921041d0a78 \
  --analysis-sha256 eb212870d61ddeaab1e5e18f6b7a95ffc83d16ee9dafc5c3986ad705d208c6f6 \
  --artifact-id 12021760871 \
  /absolute/path/extracted/results
```

The preregistered selector and repeatability diagnostics are separate fields
in `selection.json`. A repeated-process discrepancy is retained in full and
does not silently change the frozen decision rule. The exact-head replay
workflow reproduces the result records from the archived ZIP and checks that
the current commit leaves them byte-identical.
