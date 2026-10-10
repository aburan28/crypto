"""Freeze the C2000 target-call screen against archived C1800 runs."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUTPUT = HERE / "SCREEN_FREEZE.json"
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
REFERENCE = ROOT / "experiments/koblitz-s3-c1800-combined-20261010/analysis.json"
REFERENCE_FREEZE = ROOT / "experiments/koblitz-s3-c1800-combined-20261010/CONFIRM_FREEZE.json"
CONTROL = HERE / "local/control/receipt.json"
SOURCE = ROOT / "examples/s3_partial_prepared_fullrank_prefetch_projective_c2000prefix_targetspan.rs"
BINARY = ROOT / "target/release/examples/s3_partial_prepared_fullrank_prefetch_projective_c2000prefix_targetspan"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode()


def main() -> None:
    if OUTPUT.exists():
        raise FileExistsError(OUTPUT)
    control = json.loads(CONTROL.read_bytes())
    assert control["scalar_verified"] and control["full_rank_before_target"]
    assert control["memory_within_envelope"]
    reference = json.loads(REFERENCE.read_bytes())
    assert reference["run_count"] == reference["verified_run_count"] == 32
    assert reference["independent_sage_replay_verified"]
    old_plan = json.loads(REFERENCE_FREEZE.read_bytes())
    labels = [f"T{i}" for i in range(205, 213)]
    assert old_plan["target_labels"] == labels
    candidate_path = next((HERE / "identities/candidates").glob("*.json"))
    candidate = json.loads(candidate_path.read_bytes())
    candidate_id = candidate.pop("candidate_id")
    assert candidate_id.endswith("h" + hashlib.sha256(canonical(candidate)).hexdigest()[:12])
    factor = candidate["factor_base"]
    flags = candidate["implementation"]["algorithm_flags"]
    assert factor["actual_usable_points_B"] == control["factor_base_points"] == 212000
    assert factor["enumerated_set_digest"] == control["factor_base_digest"]
    assert flags["root_table_slots"] == control["root_table_slots"] == 536870912
    assert flags["prepare_rank"] == 2000
    assert control["full_rank_before_target"]
    assert candidate["implementation"]["example_source_sha256"] == sha(SOURCE)
    arm = {
        "candidate_id": candidate_id,
        "candidate_manifest_relative": str(candidate_path.relative_to(ROOT)),
        "candidate_manifest_sha256": sha(candidate_path),
        "source_relative": str(SOURCE.relative_to(ROOT)),
        "source_sha256": sha(SOURCE),
        "binary_sha256": sha(BINARY),
        "columns": 2000,
        "prepare_rank": 2000,
        "factor_base_points": 212000,
        "factor_base_digest": control["factor_base_digest"],
        "root_table_slots": flags["root_table_slots"],
        "root_table_layout": flags["root_table_layout"],
        "root_table_capacity_policy": flags["root_table_capacity_policy"],
        "target_search_schedule": flags["target_search_schedule"],
        "target_search_prefix_columns": flags["target_search_prefix_columns"],
    }
    assert arm["source_sha256"] == control["source_sha256"]
    assert arm["binary_sha256"] == control["binary_sha256"]
    reference_calls = {point["label"]: point["combined_target_s3_calls"]
                       for point in reference["points"]}
    assert set(reference_calls) == set(labels)
    plan = {
        "kind": "s3_n53_c2000_combined_screen_freeze_v1",
        "target_labels": labels,
        "targets": old_plan["targets"],
        "arms": {"combined": arm},
        "reference_candidate_id": old_plan["arms"]["combined"]["candidate_id"],
        "reference_analysis_sha256": sha(REFERENCE),
        "reference_freeze_sha256": sha(REFERENCE_FREEZE),
        "reference_target_s3_calls": reference_calls,
        "parent_freeze_sha256": sha(PARENT),
        "control_freeze_sha256": sha(HERE / "CONTROL_FREEZE.json"),
        "control_receipt_sha256": sha(CONTROL),
        "relation_seed": 20260928,
        "relation_workers": 14,
        "target_workers": 1,
        "repetitions_per_arm_per_target": 1,
        "run_order": labels,
        "timeout_seconds_per_run": 600,
        "memory_envelope_bytes": 34359738368,
        "advance_root_call_geomean_ratio_min": "1.10",
        "single_runner_sha256": sha(HERE / "run_one.py"),
        "screen_runner_sha256": sha(HERE / "run_screen.py"),
        "host_exclusive_partition_receipt": None,
        "controlled_speedup": None,
    }
    OUTPUT.write_text(json.dumps(plan, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"freeze_sha256": sha(OUTPUT), "candidate_id": candidate_id,
                      "target_count": len(labels)}, sort_keys=True))


if __name__ == "__main__":
    main()
