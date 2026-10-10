"""Freeze a matched eight-point batch-8/batch-32 N53 screen."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUTPUT = HERE / "PAIR_FREEZE.json"
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
GENERATED = ROOT / "experiments/koblitz-s3-c1800-combined-20261010/GENERATED_TARGETS.json"
REFERENCE_FREEZE = ROOT / "experiments/koblitz-s3-c1800-combined-20261010/CONFIRM_FREEZE.json"
CONTROL = HERE / "local/control/receipt.json"

ARMS = {
    "baseline": {
        "candidate": "experiments/koblitz-s3-c1800-combined-20261010/identities/candidates/IC1N53Ckb1fb190800PDP4rootRCsampleLAgaussTDdirectISO0h1ab23390e700.json",
        "source": "examples/s3_partial_prepared_fullrank_prefetch_projective_c1800prefix_targetspan.rs",
        "binary": "target/release/examples/s3_partial_prepared_fullrank_prefetch_projective_c1800prefix_targetspan",
        "report_inversion_policy": "one normal-basis field inversion per eight regular state pair products",
    },
    "candidate": {
        "candidate": "experiments/koblitz-s3-c1800-batch32-20261010/identities/candidates/IC1N53Ckb1fb190800PDP4rootRCsampleLAgaussTDdirectISO0hd14498b06bee.json",
        "source": "examples/s3_partial_prepared_fullrank_prefetch_projective_c1800prefix_targetspan_batch32.rs",
        "binary": "target/release/examples/s3_partial_prepared_fullrank_prefetch_projective_c1800prefix_targetspan_batch32",
        "report_inversion_policy": "one normal-basis field inversion per 32 regular state pair products",
    },
}


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def canonical(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode()


def main() -> None:
    if OUTPUT.exists():
        raise FileExistsError(OUTPUT)
    generated_raw = GENERATED.read_bytes()
    generated = json.loads(generated_raw)
    targets = generated["targets"]
    labels = [f"T{i}" for i in range(205, 213)]
    assert [row["label"] for row in targets] == labels
    reference = json.loads(REFERENCE_FREEZE.read_bytes())
    assert reference["target_labels"] == labels
    assert targets == reference["targets"]
    control = json.loads(CONTROL.read_bytes())
    assert control["scalar_verified"] and control["full_rank_before_target"]
    assert control["memory_within_envelope"]
    arms = {}
    for name, item in ARMS.items():
        manifest_path = ROOT / item["candidate"]
        raw = manifest_path.read_bytes()
        manifest = json.loads(raw)
        candidate_id = manifest.pop("candidate_id")
        assert candidate_id.endswith("h" + sha(canonical(manifest))[:12])
        source_sha = sha((ROOT / item["source"]).read_bytes())
        assert source_sha == manifest["implementation"]["example_source_sha256"]
        flags = manifest["implementation"]["algorithm_flags"]
        assert flags["target_query_batch_width"] == (8 if name == "baseline" else 32)
        assert manifest["factor_base"]["actual_usable_points_B"] == 190800
        arms[name] = {
            "candidate_id": candidate_id,
            "candidate_manifest_relative": item["candidate"],
            "candidate_manifest_sha256": sha(raw),
            "source_relative": item["source"],
            "source_sha256": source_sha,
            "binary_sha256": sha((ROOT / item["binary"]).read_bytes()),
            "columns": 1800,
            "prepare_rank": 1800,
            "factor_base_points": 190800,
            "factor_base_digest": manifest["factor_base"]["enumerated_set_digest"],
            "root_table_slots": flags["root_table_slots"],
            "root_table_layout": flags["root_table_layout"],
            "root_table_capacity_policy": flags["root_table_capacity_policy"],
            "target_search_schedule": flags["target_search_schedule"],
            "target_search_prefix_columns": flags["target_search_prefix_columns"],
            "target_query_batch_width": flags["target_query_batch_width"],
            "report_target_query_inversion_policy": item["report_inversion_policy"],
            "target_workers": flags["target_query_workers"],
        }
    assert arms["baseline"]["factor_base_digest"] == arms["candidate"]["factor_base_digest"]
    assert arms["baseline"]["root_table_slots"] == arms["candidate"]["root_table_slots"]
    assert arms["candidate"]["source_sha256"] == control["source_sha256"]
    assert arms["candidate"]["binary_sha256"] == control["binary_sha256"]
    for target in targets:
        point_raw = (ROOT / target["point_file_relative"]).read_bytes()
        workload_raw = (ROOT / target["workload_manifest_relative"]).read_bytes()
        assert sha(point_raw) == target["point_sha256"]
        assert sha(workload_raw) == target["workload_manifest_sha256"]
    plan = {
        "kind": "s3_n53_c1800_batch32_pair_freeze_v1",
        "arms": arms,
        "comparison_variable": "target_state_inversion_and_prefetch_batch_width_8_vs_32",
        "target_labels": labels,
        "targets": targets,
        "parent_freeze_sha256": sha(PARENT.read_bytes()),
        "generated_targets_sha256": sha(generated_raw),
        "reference_freeze_sha256": sha(REFERENCE_FREEZE.read_bytes()),
        "control_freeze_sha256": sha((HERE / "CONTROL_FREEZE.json").read_bytes()),
        "control_receipt_sha256": sha(CONTROL.read_bytes()),
        "relation_seed": 20260928,
        "relation_workers": 14,
        "target_workers": 1,
        "repetitions_per_arm_per_target": 1,
        "run_order": "baseline_candidate_even_target_candidate_baseline_odd_target",
        "timeout_seconds_per_run": 180,
        "memory_envelope_bytes": 34359738368,
        "single_runner_sha256": sha((HERE / "run_one.py").read_bytes()),
        "paired_runner_sha256": sha((HERE / "run_paired.py").read_bytes()),
        "analyzer_sha256": sha((HERE / "analyze.py").read_bytes()),
        "advance_online_geomean_ratio_min": "1.10",
        "host_exclusive_partition_receipt": None,
        "controlled_speedup": None,
    }
    OUTPUT.write_text(json.dumps(plan, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"freeze_sha256": sha(OUTPUT.read_bytes()),
                      "target_count": len(labels), "run_count": 16,
                      "candidate_id": arms["candidate"]["candidate_id"]}, sort_keys=True))


if __name__ == "__main__":
    main()
