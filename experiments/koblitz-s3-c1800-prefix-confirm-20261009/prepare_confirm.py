"""Freeze a fresh confirmation without changing either candidate identity."""

import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SCREEN = ROOT / "experiments/koblitz-s3-c1800-prefix-20261009/SCREEN_FREEZE.json"
RUNNER = HERE / "run_confirm.py"


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def main():
    screen_raw = SCREEN.read_bytes()
    screen = json.loads(screen_raw)
    assert screen["target_labels"] == [f"T{i:03}" for i in range(109, 117)]
    assert screen["controlled_speedup"] is None
    confirmation = dict(screen)
    confirmation.update({
        "kind": "s3_n53_c1200_c1800prefix_paired_confirmation_freeze_v1",
        "screen_freeze_sha256": sha(screen_raw),
        "runner_sha256": sha(RUNNER.read_bytes()),
        "target_labels": [f"T{i:03}" for i in range(117, 125)],
        "confirmation_root_call_geomean_ratio_min": 2,
        "confirmation_root_call_bootstrap95_lower_strictly_gt": 1,
    })
    del confirmation["screen_advance_root_call_geomean_ratio_min"]
    (HERE / "CONFIRM_FREEZE.json").write_text(
        json.dumps(confirmation, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"targets": confirmation["target_labels"],
                      "runner_sha256": confirmation["runner_sha256"]}, sort_keys=True))


if __name__ == "__main__":
    main()
