"""Freeze the C3200 mixed-quotient prefetch correctness control."""

from __future__ import annotations

import hashlib
import gzip
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUTPUT = HERE / "CONTROL_FREEZE.json"
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
REFERENCE = ROOT / "experiments/koblitz-s3-c3200-projective-20261010/local/verified_t025/receipt.json"
REFERENCE_REPORT = ROOT / "experiments/koblitz-s3-c3200-projective-20261010/local/verified_t025/report.json.gz"
SOURCE = ROOT / "examples/s3_partial_prepared_fullrank_c3200projective_prefetch_targetspan.rs"
BINARY = ROOT / "target/release/examples/s3_partial_prepared_fullrank_c3200projective_prefetch_targetspan"


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def main() -> None:
    if OUTPUT.exists():
        raise FileExistsError(OUTPUT)
    parent_raw = PARENT.read_bytes()
    parent = json.loads(parent_raw)
    reference_raw = REFERENCE.read_bytes()
    reference = json.loads(reference_raw)
    point = next(row for row in parent["targets"]
                 if row["n"] == 53 and row["label"] == "T025")
    reference_report_raw = gzip.decompress(REFERENCE_REPORT.read_bytes())
    reference_report = json.loads(reference_report_raw)
    assert reference["scalar_verified"] and reference["full_rank_before_target"]
    assert reference["factor_base_points"] == 339200
    assert reference["point_sha256"] == point["point_sha256"]
    assert reference_report["factor_base_digest"] == reference["factor_base_digest"]
    assert reference_report["target_relation_indices"]
    plan = {
        "kind": "s3_n53_c3200_prefetch_control_freeze_v1",
        "candidate_id": None,
        "actual_usable_points_B": None,
        "parent_freeze_sha256": sha(parent_raw),
        "reference_control_receipt_sha256": sha(reference_raw),
        "reference_report_gzip_sha256": sha(REFERENCE_REPORT.read_bytes()),
        "reference_report_raw_sha256": sha(reference_report_raw),
        "source_reference_sha256": reference["source_sha256"],
        "source_sha256": sha(SOURCE.read_bytes()),
        "binary_sha256": sha(BINARY.read_bytes()),
        "runner_sha256": sha((HERE / "run_control.py").read_bytes()),
        "curve_id": parent["curves"]["53"]["curve_id"],
        "field_degree_n": 53,
        "orbit_columns": 3200,
        "prepare_rank": 3200,
        "relation_seed": parent["relation_seed"],
        "relation_workers": 14,
        "target_workers": 1,
        "target_query_batch_width": 8,
        "target_lookup_prefetch_policy": "prefetch_all_mixed_quotient_home_slots_after_eight_state_root_batch_then_original_lookup_order",
        "target_search_prefix_columns": 1200,
        "target_search_schedule": "c1200_embedded_state_prefix_then_complete_c3200",
        "target_label": "T025",
        "point_file_relative": "holdout/n53/T025/public_target.json",
        "point_sha256": point["point_sha256"],
        "public_point": point["public_point"],
        "fixture_scalar": point["fixture_scalar"],
        "expected_factor_base_points": 339200,
        "expected_factor_base_digest": reference["factor_base_digest"],
        "expected_root_table_slots": 1073741824,
        "expected_first_relation_indices": reference_report["target_relation_indices"],
        "timeout_seconds": 1200,
        "memory_envelope_bytes": 34359738368,
        "controlled_speedup": None,
        "host_exclusive_partition_receipt": None,
    }
    OUTPUT.write_text(json.dumps(plan, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"freeze_sha256": sha(OUTPUT.read_bytes()),
                      "source_sha256": plan["source_sha256"],
                      "binary_sha256": plan["binary_sha256"]}, sort_keys=True))


if __name__ == "__main__":
    main()
