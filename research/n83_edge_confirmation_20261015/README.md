# N83 fresh-point confirmation

The frozen final tournament left K3425/root4/standard and
K3450/root11/wide3 unresolved at the declared 99% paired boundary. The
separately replayed root11 K3450/K3500/K3550 edge screen admitted no new
challenger under `PROTOCOL.md`. `admit.py` binds those source selections and
their artifact hashes to `admission.json`; it refuses a changed admission.

`generate.py` creates the already committed fresh public-point law and its
answer sidecar. The four arms run in candidate-list and reverse-list order,
using fresh processes and indexes on the same cloud allocation. The exact
producer source SHA-256, candidate manifests, source-screen packages,
resource guard and target fixture are checked before execution. The workflow
archives every arm and independently replays the base, four-point relations,
subgroup scalars, source receipt and online phase sums.

Run the local preflight after restoring the historical public fixtures:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 research/n83_edge_confirmation_20261015/generate.py
PYTHONDONTWRITEBYTECODE=1 python3 research/n83_edge_confirmation_20261015/admit.py
PYTHONDONTWRITEBYTECODE=1 python3 research/n83_edge_confirmation_20261015/run.py --preflight
```

The GitHub workflow `.github/workflows/n83-edge-confirmation.yml` executes
the frozen cloud run. Its raw artifact must be archived, checked against the
GitHub digest, replayed locally, and sealed before promoting any selection.
