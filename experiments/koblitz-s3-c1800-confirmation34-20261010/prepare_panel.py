"""Freeze a fresh 34-target N53 threshold-confirmation design."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import random


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
REFERENCE = ROOT / "experiments/koblitz-s3-c1800-combined-20261010/CONFIRM_FREEZE.json"
ROBUSTNESS = ROOT / "experiments/koblitz-s3-c1800-combined-20261010/robustness.json"
GENERATED = HERE / "GENERATED_TARGETS.json"
OUTPUT = HERE / "PANEL_FREEZE.json"


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def canonical(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode()


def main() -> None:
    if OUTPUT.exists():
        raise FileExistsError(OUTPUT)
    parent_raw = PARENT.read_bytes()
    reference_raw = REFERENCE.read_bytes()
    robustness_raw = ROBUSTNESS.read_bytes()
    generated_raw = GENERATED.read_bytes()
    parent = json.loads(parent_raw)
    reference = json.loads(reference_raw)
    robustness = json.loads(robustness_raw)
    generated = json.loads(generated_raw)
    assert generated["kind"] == "s3_n53_c1800_confirmation34_public_targets_v1"
    assert generated["parent_freeze_sha256"] == sha(parent_raw)
    assert generated["generator_order_verified"]
    assert generated["existing_t200_fixture_replayed"]
    assert robustness["sample_size_planning"]["estimated_distinct_targets_for_80_percent_power"] == 34
    assert robustness["sample_size_planning"]["process_runs_per_target"] == 4
    assert reference["relation_seed"] == parent["relation_seed"] == 20260928
    assert reference["target_workers"] == 1
    assert reference["relation_workers"] == 14
    targets = generated["targets"]
    labels = [f"T{i}" for i in range(213, 247)]
    assert [row["label"] for row in targets] == labels
    assert len({row["fixture_scalar"] for row in targets}) == 34
    assert len({tuple(row["public_point"]) for row in targets}) == 34
    old_targets = [row for row in parent["targets"] if row["n"] == 53]
    old_targets.extend(reference["targets"])
    old_scalars = {row["fixture_scalar"] for row in old_targets}
    old_points = {tuple(row["public_point"]) for row in old_targets}
    assert not old_scalars.intersection(row["fixture_scalar"] for row in targets)
    assert not old_points.intersection(tuple(row["public_point"]) for row in targets)
    for target in targets:
        point_raw = (ROOT / target["point_file_relative"]).read_bytes()
        workload_raw = (ROOT / target["workload_manifest_relative"]).read_bytes()
        assert sha(point_raw) == target["point_sha256"]
        assert json.loads(point_raw) == target["public_point"]
        assert sha(workload_raw) == target["workload_manifest_sha256"]
        workload = json.loads(workload_raw)
        assert workload["workload_id"] == target["workload_id"]
        assert workload["identity_sha256"] == sha(canonical(workload["record"]))
        assert workload["record"]["target_count"] == 1
        assert workload["record"]["target_point"] == target["public_point"]
    arms = {}
    for name in ("baseline", "combined"):
        old = reference["arms"][name]
        manifest_raw = (ROOT / old["candidate_manifest_relative"]).read_bytes()
        source_raw = (ROOT / old["source_relative"]).read_bytes()
        assert sha(manifest_raw) == old["candidate_manifest_sha256"]
        assert sha(source_raw) == old["source_sha256"]
        manifest = json.loads(manifest_raw)
        candidate_id = manifest.pop("candidate_id")
        assert candidate_id == old["candidate_id"]
        assert candidate_id.endswith("h" + sha(canonical(manifest))[:12])
        arms[name] = {key: value for key, value in old.items()
                      if key != "binary_sha256"}
        arms[name]["linux_binary_sha256"] = None
    assert arms["baseline"]["factor_base_points"] == 127200
    assert arms["combined"]["factor_base_points"] == 190800
    shuffled = labels.copy()
    random.Random(20261010).shuffle(shuffled)
    run_order = []
    for index, label in enumerate(shuffled):
        sequence = (("baseline", 1), ("combined", 1), ("combined", 2),
                    ("baseline", 2)) if index % 2 == 0 else (
                    ("combined", 1), ("baseline", 1), ("baseline", 2),
                    ("combined", 2))
        run_order.extend({"label": label, "arm": arm, "repetition": repetition}
                         for arm, repetition in sequence)
    assert len(run_order) == 136
    plan = {
        "kind": "s3_n53_c1800_confirmation34_panel_freeze_v1",
        "parent_freeze_sha256": sha(parent_raw),
        "reference_freeze_sha256": sha(reference_raw),
        "robustness_sha256": sha(robustness_raw),
        "generated_targets_sha256": sha(generated_raw),
        "fixture_generator_sha256": sha((HERE / "make_targets_sage.py").read_bytes()),
        "fixture_sage_runtime_sha256": sha((HERE / "sage_runtime_info.json").read_bytes()),
        "panel_preparation_sha256": sha(Path(__file__).read_bytes()),
        "target_labels": labels,
        "targets": targets,
        "arms": arms,
        "comparison_variable": "factor_base_and_target_extraction_kernel",
        "distinct_target_count": 34,
        "repetitions_per_arm_per_target": 2,
        "run_order_policy": "seed_20261010_shuffled_targets_alternating_ABBA_BAAB",
        "run_order": run_order,
        "relation_seed": 20260928,
        "relation_workers": 14,
        "target_workers": 1,
        "timeout_seconds_per_process": 600,
        "observed_peak_rss_envelope_bytes": 34359738368,
        "primary_metric": "verified_single_target_online_wall_ms",
        "speedup_threshold": "2.0",
        "confirmation_gate": {
            "all_136_processes_verified": True,
            "all_34_targets_independently_replayed": True,
            "host_isolation_and_noise_gates_pass": True,
            "geometric_online_ratio_at_least": "2.0",
            "two_sided_95_percent_target_bootstrap_lower_strictly_above": "2.0",
        },
        "linux_binaries_and_execution_receipt": None,
        "host_exclusive_partition_receipt": None,
        "controlled_speedup": None,
    }
    OUTPUT.write_text(json.dumps(plan, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"panel_sha256": sha(OUTPUT.read_bytes()),
                      "target_count": len(targets), "process_runs": len(run_order)},
                     sort_keys=True))


if __name__ == "__main__":
    main()
