"""Preserve bounded-scatter feasibility and derive exact candidate/workload IDs."""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
LOCAL = HERE / "local"
PLAN = HERE / "FEASIBILITY_FREEZE.json"
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
OLD = ROOT / "experiments/koblitz-s3-c2800-prefix-20261009/identities/candidates/IC1N53Ckb1fb296800PDP4rootRCsampleLAgaussTDdirectISO0h1f41cb1aeaa9.json"
OLD_REPORT = ROOT / "experiments/koblitz-s3-c2800-prefix-20261009/local/feasibility_report.json.gz"
TEMPLATE = ROOT / "experiments/koblitz-s3-c1800-20261009/identities/workloads/W25901ccca0ba.json"
TEMP_RECEIPT = Path("/Volumes/SSD990/llm/tmp/s3-c2800stream-feasibility-receipt.json")
TEMP_REPORT = Path("/Volumes/SSD990/llm/tmp/s3-c2800stream-feasibility-report.json.gz")


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()


def replace_exact(value, changes):
    if isinstance(value, dict):
        return {key: replace_exact(item, changes) for key, item in value.items()}
    if isinstance(value, list):
        return [replace_exact(item, changes) for item in value]
    return changes.get(value, value)


def main():
    plan_raw = PLAN.read_bytes()
    plan = json.loads(plan_raw)
    receipt_raw = TEMP_RECEIPT.read_bytes()
    receipt = json.loads(receipt_raw)
    report_gz = TEMP_REPORT.read_bytes()
    raw = gzip.decompress(report_gz)
    report = json.loads(raw) if raw else None
    assert receipt["freeze_sha256"] == sha(plan_raw)
    assert receipt["source_sha256"] == plan["source_sha256"]
    assert receipt["binary_sha256"] == plan["binary_sha256"]
    assert receipt["raw_report_sha256"] == (sha(raw) if raw else None)
    LOCAL.mkdir(parents=True, exist_ok=True)
    (LOCAL / "feasibility_receipt.json").write_bytes(receipt_raw)
    (LOCAL / "feasibility_report.json.gz").write_bytes(report_gz)
    if not (receipt["exit_code"] == 0 and receipt["error"] is None
            and receipt["scalar_verified"] and receipt["full_rank_before_target"]
            and receipt["memory_within_envelope"] and receipt["root_slots_within_declared_step"]
            and report is not None):
        print(json.dumps({"feasibility_passed": False,
                          "receipt": str(LOCAL / "feasibility_receipt.json")}, sort_keys=True))
        return
    assert report["regular_states"] == plan["expected_regular_states_if_all_regular"]
    assert report["root_table_slots"] <= plan["root_table_slot_step_limit"]
    assert report["target_search_schedule"] == "c1200_embedded_state_prefix_then_complete_c2800"
    assert report["factor_base_points"] == receipt["factor_base_points"]
    old_report = json.loads(gzip.decompress(OLD_REPORT.read_bytes()))
    assert sha(gzip.decompress(OLD_REPORT.read_bytes())) == plan["reference_raw_report_sha256"]
    old_fields = dict(old_report)
    new_fields = dict(report)
    for key in ("timing_ms", "root_table_capacity_policy",
                "root_scatter_policy", "root_scatter_chunk_state_limit"):
        old_fields.pop(key, None)
        new_fields.pop(key, None)
    semantic_equal = old_fields == new_fields
    (LOCAL / "semantic_comparison.json").write_text(json.dumps({
        "kind": "s3_n53_c2800stream_old_vs_bounded_scatter_semantic_control_v1",
        "old_raw_sha256": plan["reference_raw_report_sha256"],
        "new_raw_sha256": receipt["raw_report_sha256"],
        "excluded_fields": ["timing_ms", "root_table_capacity_policy",
                            "root_scatter_policy", "root_scatter_chunk_state_limit"],
        "non_timing_fields_equal": semantic_equal,
        "old_only_keys": sorted(set(old_fields) - set(new_fields)),
        "new_only_keys": sorted(set(new_fields) - set(old_fields)),
        "differing_keys": sorted(key for key in old_fields.keys() & new_fields.keys()
                                 if old_fields[key] != new_fields[key]),
    }, indent=2, sort_keys=True) + "\n")
    if not semantic_equal:
        print(json.dumps({"feasibility_passed": False,
                          "semantic_comparison": str(LOCAL / "semantic_comparison.json")},
                         sort_keys=True))
        return

    old = json.loads(OLD.read_bytes())
    old_id = old.pop("candidate_id")
    assert old_id.endswith("h" + sha(canonical(old))[:12])
    old_source = old["implementation"]["example_source_sha256"]
    old_commit = old["implementation"]["frozen_source_commit"]
    old_policy = old["point_decomposition"]["root_table_capacity_policy"]
    candidate = replace_exact(old, {
        old_source: plan["source_sha256"],
        old_commit: plan["source_commit"],
        old_policy: report["root_table_capacity_policy"],
    })
    for section in (candidate["point_decomposition"],
                    candidate["implementation"]["algorithm_flags"]):
        section["root_scatter_policy"] = report["root_scatter_policy"]
        section["root_scatter_chunk_state_limit"] = report["root_scatter_chunk_state_limit"]
    assert candidate["factor_base"]["actual_usable_points_B"] == report["factor_base_points"]
    assert candidate["factor_base"]["effective_relation_columns"] == 2800
    assert candidate["point_decomposition"]["target_search_schedule"] == report["target_search_schedule"]
    assert candidate["implementation"]["algorithm_flags"]["root_table_slots"] == report["root_table_slots"]
    assert candidate["implementation"]["algorithm_flags"]["index_state_payload_bytes"] == report["index_state_payload_bytes"]
    assert candidate["factor_base"]["enumerated_set_digest"] == report["factor_base_digest"]
    candidate_id = ("IC1N53Ckb1fb" + str(report["factor_base_points"])
                    + "PDP4rootRCsampleLAgaussTDdirectISO0h" + sha(canonical(candidate))[:12])
    candidate["candidate_id"] = candidate_id
    identities = HERE / "identities/candidates"
    identities.mkdir(parents=True, exist_ok=True)
    (identities / (candidate_id + ".json")).write_bytes(canonical(candidate) + b"\n")

    target_template = json.loads(TEMPLATE.read_bytes())
    target_template.pop("workload_id")
    parent = json.loads(PARENT.read_bytes())
    targets = {item["label"]: item for item in parent["targets"] if item["n"] == 53}
    workloads = HERE / "identities/workloads"
    workloads.mkdir(parents=True, exist_ok=True)
    for number in range(157, 165):
        label = f"T{number:03}"
        target = targets[label]
        assert sha((PARENT.parent / target["point_file_relative"]).read_bytes()) == target["point_sha256"]
        record = dict(target_template)
        record["target_index"] = number
        record["target_label"] = label
        record["target_point"] = target["public_point"]
        record["scalar_derivation_sha256"] = target["scalar_derivation_sha256"]
        workload_id = sha(canonical(record))[:12]
        assert workload_id == target["workload_id"]
        record["workload_id"] = workload_id
        (workloads / ("W" + workload_id + ".json")).write_bytes(canonical(record) + b"\n")
    print(json.dumps({"candidate_id": candidate_id,
                      "factor_base_points": report["factor_base_points"],
                      "peak_rss_bytes": receipt["peak_rss_bytes"],
                      "workload_count": 8}, sort_keys=True))


if __name__ == "__main__":
    main()
