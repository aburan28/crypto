"""Freeze the paired N53 C1200/C1800 combined one-target panel."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
GENERATED = HERE / "GENERATED_TARGETS.json"
CONTROL = HERE / "local/control/receipt.json"
OUTPUT = HERE / "CONFIRM_FREEZE.json"

ARMS = {
    "baseline": {
        "source": "examples/s3_partial_prepared_fullrank_c2800mixed_targetspan.rs",
        "binary": "target/release/examples/s3_partial_prepared_fullrank_c2800mixed_targetspan",
        "candidate": "experiments/koblitz-s3-c2800-mixed-20261009/identities/candidates/IC1N53Ckb1fb127200PDP4rootRCsampleLAgaussTDdirectISO0hf8b850f65ece.json",
        "columns": 1200,
    },
    "combined": {
        "source": "examples/s3_partial_prepared_fullrank_prefetch_projective_c1800prefix_targetspan.rs",
        "binary": "target/release/examples/s3_partial_prepared_fullrank_prefetch_projective_c1800prefix_targetspan",
        "candidate": "experiments/koblitz-s3-c1800-combined-20261010/identities/candidates/IC1N53Ckb1fb190800PDP4rootRCsampleLAgaussTDdirectISO0h1ab23390e700.json",
        "columns": 1800,
    },
}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode()


def main() -> None:
    if OUTPUT.exists():
        raise FileExistsError(OUTPUT)
    parent = json.loads(PARENT.read_bytes())
    assert parent["input_law"].startswith("1 + int_be(SHA256(")
    generated = json.loads(GENERATED.read_bytes())
    targets = generated["targets"]
    assert [row["label"] for row in targets] == [f"T{k}" for k in range(205, 213)]
    assert generated["parent_freeze_sha256"] == sha(PARENT)
    assert generated["generator_order_verified"]
    control = json.loads(CONTROL.read_bytes())
    assert control["scalar_verified"] and control["full_rank_before_target"]
    assert control["memory_within_envelope"]
    assert control["factor_base_points"] == 190800
    assert control["factor_base_digest"] == (
        "0c3597d24ca943d7bb081ee30078596ffc433766356a8042ba8180c2ced78030")
    arms = {}
    for label, item in ARMS.items():
        candidate_path = ROOT / item["candidate"]
        candidate = json.loads(candidate_path.read_bytes())
        candidate_id = candidate.pop("candidate_id")
        assert candidate_id.endswith("h" + hashlib.sha256(canonical(candidate)).hexdigest()[:12])
        assert candidate["field"]["degree_n"] == 53
        assert candidate["factor_base"]["effective_relation_columns"] == item["columns"]
        flags = candidate["implementation"]["algorithm_flags"]
        arm = {
            "candidate_id": candidate_id,
            "candidate_manifest_relative": item["candidate"],
            "candidate_manifest_sha256": sha(candidate_path),
            "source_relative": item["source"],
            "source_sha256": sha(ROOT / item["source"]),
            "binary_sha256": sha(ROOT / item["binary"]),
            "columns": item["columns"],
            "prepare_rank": flags["prepare_rank"],
            "factor_base_points": candidate["factor_base"]["actual_usable_points_B"],
            "factor_base_digest": candidate["factor_base"]["enumerated_set_digest"],
            "root_table_slots": flags["root_table_slots"],
            "root_table_layout": flags["root_table_layout"],
            "root_table_capacity_policy": flags["root_table_capacity_policy"],
            "target_search_schedule": flags["target_search_schedule"],
            "target_search_prefix_columns": flags["target_search_prefix_columns"],
            "target_workers": flags["target_query_workers"],
        }
        assert arm["source_sha256"] == candidate["implementation"]["example_source_sha256"]
        arms[label] = arm
    assert arms["combined"]["source_sha256"] == control["source_sha256"]
    assert arms["combined"]["binary_sha256"] == control["binary_sha256"]
    assert arms["combined"]["factor_base_digest"] == control["factor_base_digest"]
    assert arms["combined"]["root_table_slots"] == control["root_table_slots"]
    for target in targets:
        point = ROOT / target["point_file_relative"]
        workload = ROOT / target["workload_manifest_relative"]
        assert sha(point) == target["point_sha256"]
        assert json.loads(point.read_bytes()) == target["public_point"]
        assert sha(workload) == target["workload_manifest_sha256"]
        record = json.loads(workload.read_bytes())
        assert record["workload_id"] == target["workload_id"]
        assert record["identity_sha256"] == hashlib.sha256(canonical(record["record"])).hexdigest()
        assert record["record"]["target_point"] == target["public_point"]
    plan = {
        "kind": "s3_n53_c1800_combined_confirmation_freeze_v1",
        "arms": arms,
        "target_labels": [row["label"] for row in targets],
        "targets": targets,
        "comparison_variable": "factor_base_1200_vs_1800_columns_and_combined_online_kernel",
        "parent_freeze_sha256": sha(PARENT),
        "generated_targets_sha256": sha(GENERATED),
        "fixture_generator_sha256": sha(HERE / "make_targets_sage.py"),
        "fixture_generator_sage_runtime_sha256": sha(HERE / "sage-runtime-info.json"),
        "control_freeze_sha256": sha(HERE / "CONTROL_FREEZE.json"),
        "control_receipt_sha256": sha(CONTROL),
        "relation_seed": 20260928,
        "relation_workers": 14,
        "target_workers": 1,
        "repetitions_per_arm_per_target": 2,
        "run_order": "ABBA_even_target_BAAB_odd_target",
        "timeout_seconds_per_run": 180,
        "memory_envelope_bytes": 34359738368,
        "single_runner_sha256": sha(HERE / "run_one.py"),
        "paired_runner_sha256": sha(HERE / "run_paired.py"),
        "analyzer_sha256": sha(HERE / "analyze.py"),
        "independent_replay_sha256": sha(HERE / "replay_sage.py"),
        "local_confirmation_geomean_min": "2",
        "local_confirmation_bootstrap_lower_min": "1",
        "robust_twofold_bootstrap_lower_min": "2",
        "host_exclusive_partition_receipt": None,
        "controlled_speedup": None,
    }
    OUTPUT.write_text(json.dumps(plan, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"freeze_sha256": sha(OUTPUT),
                      "target_count": len(targets),
                      "baseline_candidate": arms["baseline"]["candidate_id"],
                      "combined_candidate": arms["combined"]["candidate_id"]},
                     sort_keys=True))


if __name__ == "__main__":
    main()
