"""Freeze matched C3200 projective/prefetch development targets and binaries."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUTPUT = HERE / "SCREEN_FREEZE.json"
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
REFERENCE = ROOT / "experiments/koblitz-s3-c3200-screen-20261010/SCREEN_FREEZE.json"
CONTROL = HERE / "local/control/receipt.json"
SAGE_REPLAY = HERE / "independent_sage_replay.json"
PREFETCH_POLICY = (
    "prefetch_all_mixed_quotient_home_slots_after_eight_state_root_batch_then_original_lookup_order")
ARMS = {
    "baseline": {
        "candidate": "experiments/koblitz-s3-c3200-projective-20261010/identities/candidates/IC1N53Ckb1fb339200PDP4rootRCsampleLAgaussTDdirectISO0h0e595bd6143a.json",
        "source": "examples/s3_partial_prepared_fullrank_c3200projective_targetspan.rs",
        "binary": "target/release/examples/s3_partial_prepared_fullrank_c3200projective_targetspan",
        "prefetch_policy": None,
    },
    "candidate": {
        "candidate": "experiments/koblitz-s3-c3200-prefetch-20261010/identities/candidates/IC1N53Ckb1fb339200PDP4rootRCsampleLAgaussTDdirectISO0h35d70065259b.json",
        "source": "examples/s3_partial_prepared_fullrank_c3200projective_prefetch_targetspan.rs",
        "binary": "target/release/examples/s3_partial_prepared_fullrank_c3200projective_prefetch_targetspan",
        "prefetch_policy": PREFETCH_POLICY,
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
    reference_raw = REFERENCE.read_bytes()
    reference = json.loads(reference_raw)
    parent_raw = PARENT.read_bytes()
    control = json.loads(CONTROL.read_bytes())
    replay = json.loads(SAGE_REPLAY.read_bytes())
    labels = [f"T{i}" for i in range(173, 181)]
    assert reference["target_labels"] == labels
    assert len(reference["targets"]) == 8
    assert reference["parent_freeze_sha256"] == sha(parent_raw)
    assert control["scalar_verified"] and control["full_rank_before_target"]
    assert control["memory_within_envelope"] and control["index_layout_ok"]
    assert control["same_first_relation_witness_as_reference"]
    assert replay["generator_order_verified"]
    assert replay["scalar_times_generator_equals_public_target"]
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
        assert flags["target_query_batch_width"] == 8
        assert flags.get("target_lookup_prefetch_policy") == item["prefetch_policy"]
        assert manifest["factor_base"]["actual_usable_points_B"] == 339200
        arms[name] = {
            "candidate_id": candidate_id,
            "candidate_manifest_relative": item["candidate"],
            "candidate_manifest_sha256": sha(raw),
            "source_relative": item["source"],
            "source_sha256": source_sha,
            "binary_sha256": sha((ROOT / item["binary"]).read_bytes()),
            "columns": 3200,
            "prepare_rank": 3200,
            "factor_base_points": 339200,
            "factor_base_digest": manifest["factor_base"]["enumerated_set_digest"],
            "root_table_slots": flags["root_table_slots"],
            "root_table_layout": flags["root_table_layout"],
            "root_table_capacity_policy": flags["root_table_capacity_policy"],
            "target_search_schedule": flags["target_search_schedule"],
            "target_search_prefix_columns": flags["target_search_prefix_columns"],
            "report_target_lookup_prefetch_policy": item["prefetch_policy"],
            "target_query_batch_width": flags["target_query_batch_width"],
        }
    for key in ("factor_base_digest", "factor_base_points", "root_table_slots",
                "root_table_layout", "root_table_capacity_policy", "columns",
                "prepare_rank", "target_search_schedule", "target_search_prefix_columns"):
        assert arms["baseline"][key] == arms["candidate"][key]
    assert arms["candidate"]["source_sha256"] == control["source_sha256"]
    assert arms["candidate"]["binary_sha256"] == control["binary_sha256"]
    for target in reference["targets"]:
        point_file = PARENT.parent / target["point_file_relative"]
        assert sha(point_file.read_bytes()) == target["point_sha256"]
        assert sha((ROOT / target["workload_manifest_relative"]).read_bytes()) == target[
            "workload_manifest_sha256"]
    plan = {
        "kind": "s3_n53_c3200_prefetch_pair_freeze_v1",
        "arms": arms,
        "comparison_variable": "batch_home_slot_prefetch_before_ordered_exact_mixed_quotient_lookups",
        "target_labels": labels,
        "targets": reference["targets"],
        "parent_freeze_sha256": sha(parent_raw),
        "reference_screen_freeze_sha256": sha(reference_raw),
        "control_freeze_sha256": sha((HERE / "CONTROL_FREEZE.json").read_bytes()),
        "control_receipt_sha256": sha(CONTROL.read_bytes()),
        "control_sage_replay_sha256": sha(SAGE_REPLAY.read_bytes()),
        "relation_seed": 20260928,
        "relation_workers": 14,
        "target_workers": 1,
        "repetitions_per_arm_per_target": 1,
        "run_order": "baseline_candidate_even_target_candidate_baseline_odd_target",
        "timeout_seconds_per_run": 600,
        "memory_envelope_bytes": 34359738368,
        "single_runner_sha256": sha((HERE / "run_one.py").read_bytes()),
        "paired_runner_sha256": sha((HERE / "run_paired.py").read_bytes()),
        "analyzer_sha256": sha((HERE / "analyze.py").read_bytes()),
        "sage_replay_source_sha256": sha((HERE / "replay_sage.py").read_bytes()),
        "advance_online_geomean_ratio_min": "2.0",
        "host_exclusive_partition_receipt": None,
        "controlled_speedup": None,
    }
    OUTPUT.write_text(json.dumps(plan, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"freeze_sha256": sha(OUTPUT.read_bytes()),
                      "target_count": len(labels), "run_count": 16,
                      "candidate_id": arms["candidate"]["candidate_id"]}, sort_keys=True))


if __name__ == "__main__":
    main()
