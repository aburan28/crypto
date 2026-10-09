"""Freeze fresh T133–T140 workloads and the unchanged striped candidates."""

import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DIAGNOSTIC = ROOT / "experiments/koblitz-s3-c1800-stripe-20261009/DIAGNOSTIC_FREEZE.json"
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
TEMPLATE = ROOT / "experiments/koblitz-s3-c1800-20261009/identities/workloads/W25901ccca0ba.json"
SOURCE = ROOT / "examples/s3_partial_prepared_fullrank_c1800stripe_targetspan.rs"
BINARY = ROOT / "target/release/examples/s3_partial_prepared_fullrank_c1800stripe_targetspan"
RUNNER = HERE / "run_paired.py"


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()


def main():
    diagnostic_raw = DIAGNOSTIC.read_bytes()
    diagnostic = json.loads(diagnostic_raw)
    assert diagnostic["source_sha256"] == sha(SOURCE.read_bytes())
    assert diagnostic["binary_sha256"] == sha(BINARY.read_bytes())
    for arm in ("baseline", "expanded"):
        frozen = diagnostic["arms"][arm]
        manifest = json.loads((ROOT / frozen["candidate_manifest"]).read_bytes())
        candidate_id = manifest.pop("candidate_id")
        assert candidate_id == frozen["candidate_id"]
        assert candidate_id.endswith("h" + sha(canonical(manifest))[:12])
    parent = json.loads(PARENT.read_bytes())
    targets = {target["label"]: target for target in parent["targets"] if target["n"] == 53}
    template = json.loads(TEMPLATE.read_bytes())
    template.pop("workload_id")
    workload_dir = HERE / "identities/workloads"
    workload_dir.mkdir(parents=True, exist_ok=True)
    for number in range(133, 141):
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
    plan = {
        "kind": "s3_n53_c1200_c1800stripe_paired_call_screen_freeze_v1",
        "parent_diagnostic_freeze_sha256": sha(diagnostic_raw),
        "parent_freeze_sha256": sha(PARENT.read_bytes()),
        "source_commit": diagnostic["source_commit"],
        "source_sha256": diagnostic["source_sha256"],
        "binary_sha256": diagnostic["binary_sha256"],
        "runner_sha256": sha(RUNNER.read_bytes()),
        "arms": diagnostic["arms"],
        "target_labels": [f"T{i:03}" for i in range(133, 141)],
        "repetitions_per_arm_per_target": 1,
        "run_order": "AB_even_target_BA_odd_target",
        "relation_seed": 20260928,
        "relation_workers": 14,
        "target_workers": 1,
        "timeout_seconds_per_run": 180,
        "memory_envelope_bytes": 34359738368,
        "screen_advance_root_call_geomean_ratio_min": 3,
        "host_exclusive_partition_receipt": None,
        "controlled_speedup": None,
    }
    (HERE / "SCREEN_FREEZE.json").write_text(json.dumps(plan, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"targets": plan["target_labels"],
                      "runner_sha256": plan["runner_sha256"]}, sort_keys=True))


if __name__ == "__main__":
    main()
