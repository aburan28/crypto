"""Freeze the C1800 full-rank row-elision control before execution."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUTPUT = HERE / "CONTROL_FREEZE.json"
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
REFERENCE = ROOT / "experiments/koblitz-s3-c1800-combined-20261010/local/control/receipt.json"
SOURCE_PARENT = ROOT / "examples/s3_partial_prepared_fullrank_prefetch_projective_c1800prefix_targetspan.rs"
SOURCE = ROOT / "examples/s3_partial_prepared_fullrank_prefetch_projective_c1800prefix_targetspan_norow.rs"
BINARY = ROOT / "target/release/examples/s3_partial_prepared_fullrank_prefetch_projective_c1800prefix_targetspan_norow"


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
    assert reference["scalar_verified"] and reference["full_rank_before_target"]
    assert reference["factor_base_points"] == 190800
    assert reference["point_sha256"] == point["point_sha256"]
    assert sha(SOURCE_PARENT.read_bytes()) == reference["source_sha256"]
    plan = {
        "kind": "s3_n53_c1800_norow_control_freeze_v1",
        "candidate_id": None,
        "actual_usable_points_B": None,
        "parent_freeze_sha256": sha(parent_raw),
        "reference_control_receipt_sha256": sha(reference_raw),
        "source_parent_sha256": reference["source_sha256"],
        "source_sha256": sha(SOURCE.read_bytes()),
        "binary_sha256": sha(BINARY.read_bytes()),
        "runner_sha256": sha((HERE / "run_control.py").read_bytes()),
        "curve_id": parent["curves"]["53"]["curve_id"],
        "field_degree_n": 53,
        "orbit_columns": 1800,
        "prepare_rank": 1800,
        "relation_seed": parent["relation_seed"],
        "relation_workers": 14,
        "target_workers": 1,
        "target_query_batch_width": 8,
        "target_full_rank_row_policy": "skip_row_when_prepared_base_logs_available",
        "target_search_prefix_columns": 1200,
        "target_search_schedule": "c1200_embedded_state_prefix_then_complete_c1800",
        "target_label": "T025",
        "point_file_relative": "holdout/n53/T025/public_target.json",
        "point_sha256": point["point_sha256"],
        "public_point": point["public_point"],
        "fixture_scalar": point["fixture_scalar"],
        "expected_factor_base_points": 190800,
        "expected_factor_base_digest": reference["factor_base_digest"],
        "root_table_slot_step_limit": 536870912,
        "timeout_seconds": 360,
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
