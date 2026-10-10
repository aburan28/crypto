"""Freeze the bounded-scatter C2800 memory, rank, and scalar control."""

import hashlib
import json
from pathlib import Path
import subprocess


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE = ROOT / "examples/s3_partial_prepared_fullrank_c2800stream_targetspan.rs"
BINARY = ROOT / "target/release/examples/s3_partial_prepared_fullrank_c2800stream_targetspan"
RUNNER = HERE / "run_feasibility.py"
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
REFERENCE_RECEIPT = ROOT / "experiments/koblitz-s3-c2800-prefix-20261009/local/feasibility_receipt.json"


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def main():
    parent = json.loads(PARENT.read_bytes())
    target = next(item for item in parent["targets"]
                  if item["label"] == "T025" and item["n"] == 53)
    source_commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT,
                                            text=True).strip()
    reference = json.loads(REFERENCE_RECEIPT.read_bytes())
    assert reference["scalar_verified"] and reference["full_rank_before_target"]
    plan = {
        "kind": "s3_n53_c2800stream_fullrank_feasibility_freeze_v1",
        "curve_id": "EC1N53Ckb1h738a572b91f2",
        "field_degree_n": 53,
        "orbit_columns": 2800,
        "prepare_rank": 2800,
        "expected_regular_states_if_all_regular": 2800 * 2800 * 53,
        "expected_input_roots_if_all_regular": 2 * 2800 * 2800 * 53,
        "root_table_slot_step_limit": 1073741824,
        "expected_root_table_slots": 1073741824,
        "root_table_capacity_policy": "uniform_initial_next_power_of_two_with_exact_growth_per_shard",
        "root_scatter_policy": "ordered_state_chunks_parallel_scatter_and_shard_insert",
        "root_scatter_chunk_state_limit": 4194304,
        "memory_envelope_bytes": 34359738368,
        "relation_seed": 20260928,
        "relation_workers": 14,
        "target_workers": 1,
        "timeout_seconds": 900,
        "resource_probe": "/usr/bin/time -l whole-process maximum resident set size",
        "target_label": target["label"],
        "point_file_relative": target["point_file_relative"],
        "public_point": target["public_point"],
        "fixture_scalar": target["fixture_scalar"],
        "point_sha256": target["point_sha256"],
        "parent_freeze_sha256": sha(PARENT.read_bytes()),
        "reference_receipt_sha256": sha(REFERENCE_RECEIPT.read_bytes()),
        "reference_raw_report_sha256": reference["raw_report_sha256"],
        "source_commit": source_commit,
        "source_sha256": sha(SOURCE.read_bytes()),
        "binary_sha256": sha(BINARY.read_bytes()),
        "runner_sha256": sha(RUNNER.read_bytes()),
        "actual_usable_points_B": None,
        "candidate_id": None,
        "host_exclusive_partition_receipt": None,
        "controlled_speedup": None,
    }
    (HERE / "FEASIBILITY_FREEZE.json").write_text(
        json.dumps(plan, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"source_commit": source_commit,
                      "expected_regular_states": plan["expected_regular_states_if_all_regular"],
                      "expected_input_roots": plan["expected_input_roots_if_all_regular"]},
                     sort_keys=True))


if __name__ == "__main__":
    main()
