"""Freeze fresh T141–T148 confirmation for the unchanged striped candidates."""

import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SCREEN = ROOT / "experiments/koblitz-s3-c1800-stripe-screen-20261009/SCREEN_FREEZE.json"
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
TEMPLATE = ROOT / "experiments/koblitz-s3-c1800-20261009/identities/workloads/W25901ccca0ba.json"
RUNNER = HERE / "run_confirm.py"


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()


def main():
    screen_raw = SCREEN.read_bytes()
    screen = json.loads(screen_raw)
    assert screen["target_labels"] == [f"T{i:03}" for i in range(133, 141)]
    assert screen["controlled_speedup"] is None
    parent = json.loads(PARENT.read_bytes())
    targets = {target["label"]: target for target in parent["targets"] if target["n"] == 53}
    template = json.loads(TEMPLATE.read_bytes())
    template.pop("workload_id")
    workload_dir = HERE / "identities/workloads"
    workload_dir.mkdir(parents=True, exist_ok=True)
    for number in range(141, 149):
        label = f"T{number:03}"
        target = targets[label]
        assert sha((PARENT.parent / target["point_file_relative"]).read_bytes()) == target["point_sha256"]
        record = dict(template)
        record["target_index"] = number
        record["target_label"] = label
        record["target_point"] = target["public_point"]
        record["scalar_derivation_sha256"] = target["scalar_derivation_sha256"]
        workload_id = sha(canonical(record))[:12]
        assert workload_id == target["workload_id"]
        record["workload_id"] = workload_id
        (workload_dir / ("W" + workload_id + ".json")).write_bytes(canonical(record) + b"\n")
    confirmation = dict(screen)
    confirmation.update({
        "kind": "s3_n53_c1200_c1800stripe_paired_confirmation_freeze_v1",
        "screen_freeze_sha256": sha(screen_raw),
        "runner_sha256": sha(RUNNER.read_bytes()),
        "target_labels": [f"T{i:03}" for i in range(141, 149)],
        "repetitions_per_arm_per_target": 2,
        "run_order": "ABBA_even_target_BAAB_odd_target",
        "confirmation_root_call_geomean_ratio_min": 2,
        "confirmation_root_call_bootstrap95_lower_strictly_gt": 1,
        "exploratory_online_geomean_ratio_min": 2,
        "exploratory_online_bootstrap95_lower_strictly_gt": 1,
    })
    del confirmation["screen_advance_root_call_geomean_ratio_min"]
    (HERE / "CONFIRM_FREEZE.json").write_text(
        json.dumps(confirmation, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"targets": confirmation["target_labels"],
                      "runner_sha256": confirmation["runner_sha256"]}, sort_keys=True))


if __name__ == "__main__":
    main()
