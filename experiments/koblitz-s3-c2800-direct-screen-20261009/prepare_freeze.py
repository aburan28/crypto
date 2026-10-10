"""Freeze fresh C1200/C2800 direct-fill identities and paired screen."""

import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE = ROOT / "examples/s3_partial_prepared_fullrank_c2800direct_targetspan.rs"
BINARY = ROOT / "target/release/examples/s3_partial_prepared_fullrank_c2800direct_targetspan"
RUNNER = HERE / "run_paired.py"
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
OLD_BASELINE = ROOT / "experiments/koblitz-s3-c2800-prefix-screen-20261009/identities/candidates/IC1N53Ckb1fb127200PDP4rootRCsampleLAgaussTDdirectISO0haeee7ec59767.json"
EXPANDED = ROOT / "experiments/koblitz-s3-c2800-direct-20261009/identities/candidates/IC1N53Ckb1fb296800PDP4rootRCsampleLAgaussTDdirectISO0hf4206807ca48.json"


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()


def replace_exact(value, old_source, new_source, old_commit, new_commit):
    if isinstance(value, dict):
        return {key: replace_exact(item, old_source, new_source, old_commit, new_commit)
                for key, item in value.items()}
    if isinstance(value, list):
        return [replace_exact(item, old_source, new_source, old_commit, new_commit)
                for item in value]
    if value == old_source:
        return new_source
    if value == old_commit:
        return new_commit
    return value


def arm_record(candidate, manifest_path):
    return {
        "candidate_id": candidate["candidate_id"],
        "candidate_manifest": str(manifest_path.relative_to(ROOT)),
        "columns": candidate["factor_base"]["effective_relation_columns"],
        "factor_base_digest": candidate["factor_base"]["enumerated_set_digest"],
        "factor_base_points": candidate["factor_base"]["actual_usable_points_B"],
        "prepare_rank": candidate["factor_base"]["effective_relation_columns"],
        "root_table_slots": candidate["implementation"]["algorithm_flags"]["root_table_slots"],
        "root_table_capacity_policy": candidate["point_decomposition"]["root_table_capacity_policy"],
        "root_scatter_policy": candidate["point_decomposition"]["root_scatter_policy"],
        "root_scatter_chunk_state_limit": candidate["point_decomposition"]["root_scatter_chunk_state_limit"],
        "state_generation_policy": candidate["point_decomposition"]["state_generation_policy"],
        "target_search_schedule": candidate["point_decomposition"]["target_search_schedule"],
        "target_search_prefix_columns": 1200,
    }


def main():
    source_sha = sha(SOURCE.read_bytes())
    expanded = json.loads(EXPANDED.read_bytes())
    expanded_id = expanded.pop("candidate_id")
    assert expanded_id.endswith("h" + sha(canonical(expanded))[:12])
    assert expanded["implementation"]["example_source_sha256"] == source_sha
    source_commit = expanded["implementation"]["frozen_source_commit"]
    expanded["candidate_id"] = expanded_id

    old = json.loads(OLD_BASELINE.read_bytes())
    old_id = old.pop("candidate_id")
    assert old_id.endswith("h" + sha(canonical(old))[:12])
    baseline = replace_exact(
        old, old["implementation"]["example_source_sha256"], source_sha,
        old["implementation"]["frozen_source_commit"], source_commit)
    for section in (baseline["point_decomposition"],
                    baseline["implementation"]["algorithm_flags"]):
        section["root_table_capacity_policy"] = expanded["point_decomposition"]["root_table_capacity_policy"]
        section["root_scatter_policy"] = expanded["point_decomposition"]["root_scatter_policy"]
        section["root_scatter_chunk_state_limit"] = expanded["point_decomposition"]["root_scatter_chunk_state_limit"]
        section["state_generation_policy"] = "variable_yield_worker_chunks_serial_append"
    baseline_id = ("IC1N53Ckb1fb127200PDP4rootRCsampleLAgaussTDdirectISO0h"
                   + sha(canonical(baseline))[:12])
    baseline["candidate_id"] = baseline_id
    baseline_path = HERE / "identities/candidates" / (baseline_id + ".json")
    baseline_path.parent.mkdir(parents=True, exist_ok=True)
    baseline_path.write_bytes(canonical(baseline) + b"\n")
    arms = {"baseline": arm_record(baseline, baseline_path),
            "expanded": arm_record(expanded, EXPANDED)}
    assert arms["baseline"]["columns"] == 1200
    assert arms["expanded"]["columns"] == 2800
    plan = {
        "kind": "s3_n53_c1200_c2800direct_paired_call_screen_freeze_v1",
        "parent_freeze_sha256": sha(PARENT.read_bytes()),
        "source_commit": source_commit,
        "source_sha256": source_sha,
        "binary_sha256": sha(BINARY.read_bytes()),
        "runner_sha256": sha(RUNNER.read_bytes()),
        "arms": arms,
        "target_labels": [f"T{i:03}" for i in range(157, 165)],
        "repetitions_per_arm_per_target": 1,
        "run_order": "AB_even_target_BA_odd_target",
        "relation_seed": 20260928,
        "relation_workers": 14,
        "target_workers": 1,
        "timeout_seconds_per_run": 600,
        "memory_envelope_bytes": 34359738368,
        "resource_probe": "/usr/bin/time -l per child; full-process RSS, instructions, cycles, and context switches",
        "screen_advance_root_call_geomean_ratio_min": 3,
        "stop_after_first_invalid_or_over_budget_run": True,
        "host_exclusive_partition_receipt": None,
        "controlled_speedup": None,
    }
    (HERE / "SCREEN_FREEZE.json").write_text(json.dumps(plan, indent=2, sort_keys=True) + "\n")
    print(json.dumps({arm: value["candidate_id"] for arm, value in arms.items()}, sort_keys=True))


if __name__ == "__main__":
    main()
