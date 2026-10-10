"""Freeze fresh C1200/C2800 candidate identities and one-repetition screen."""

import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE = ROOT / "examples/s3_partial_prepared_fullrank_c2800prefix_targetspan.rs"
BINARY = ROOT / "target/release/examples/s3_partial_prepared_fullrank_c2800prefix_targetspan"
RUNNER = HERE / "run_paired.py"
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
OLD_BASELINE = ROOT / "experiments/koblitz-s3-c1800-prefix-20261009/identities/candidates/IC1N53Ckb1fb127200PDP4rootRCsampleLAgaussTDdirectISO0hee573707e568.json"
EXPANDED = ROOT / "experiments/koblitz-s3-c2800-prefix-20261009/identities/candidates/IC1N53Ckb1fb296800PDP4rootRCsampleLAgaussTDdirectISO0h1f41cb1aeaa9.json"


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
        "kind": "s3_n53_c1200_c2800prefix_paired_call_screen_freeze_v1",
        "parent_freeze_sha256": sha(PARENT.read_bytes()),
        "source_commit": source_commit,
        "source_sha256": source_sha,
        "binary_sha256": sha(BINARY.read_bytes()),
        "runner_sha256": sha(RUNNER.read_bytes()),
        "arms": arms,
        "target_labels": [f"T{i:03}" for i in range(149, 157)],
        "repetitions_per_arm_per_target": 1,
        "run_order": "AB_even_target_BA_odd_target",
        "relation_seed": 20260928,
        "relation_workers": 14,
        "target_workers": 1,
        "timeout_seconds_per_run": 600,
        "memory_envelope_bytes": 34359738368,
        "resource_probe": "/usr/bin/time -l per child; full-process RSS, instructions, cycles, and context switches",
        "screen_advance_root_call_geomean_ratio_min": 3,
        "host_exclusive_partition_receipt": None,
        "controlled_speedup": None,
    }
    (HERE / "SCREEN_FREEZE.json").write_text(json.dumps(plan, indent=2, sort_keys=True) + "\n")
    print(json.dumps({arm: value["candidate_id"] for arm, value in arms.items()}, sort_keys=True))


if __name__ == "__main__":
    main()
