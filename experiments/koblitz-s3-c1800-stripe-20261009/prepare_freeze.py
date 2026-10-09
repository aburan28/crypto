"""Freeze C1200/C1800 striped identities and a used-point control."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE = ROOT / "examples/s3_partial_prepared_fullrank_c1800stripe_targetspan.rs"
RUNNER = HERE / "run_diagnostic.py"
BINARY = ROOT / "target/release/examples/s3_partial_prepared_fullrank_c1800stripe_targetspan"
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
OLD_MANIFESTS = {
    "baseline": ROOT / "experiments/koblitz-s3-c1800-prefix-20261009/identities/candidates/IC1N53Ckb1fb127200PDP4rootRCsampleLAgaussTDdirectISO0hee573707e568.json",
    "expanded": ROOT / "experiments/koblitz-s3-c1800-prefix-20261009/identities/candidates/IC1N53Ckb1fb190800PDP4rootRCsampleLAgaussTDdirectISO0h98f1ea052683.json",
}


def digest(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def canonical(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()


def replace_sources(value: object, old_sha: str, new_sha: str, old_commit: str,
                    new_commit: str) -> object:
    if isinstance(value, dict):
        return {key: replace_sources(item, old_sha, new_sha, old_commit, new_commit)
                for key, item in value.items()}
    if isinstance(value, list):
        return [replace_sources(item, old_sha, new_sha, old_commit, new_commit)
                for item in value]
    if value == old_sha:
        return new_sha
    if value == old_commit:
        return new_commit
    return value


def main():
    source_commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT,
                                            text=True).strip()
    source_sha = digest(SOURCE.read_bytes())
    binary_sha = digest(BINARY.read_bytes())
    arms = {}
    for arm, path in OLD_MANIFESTS.items():
        old = json.loads(path.read_bytes())
        old_id = old.pop("candidate_id")
        assert old_id.endswith("h" + digest(canonical(old))[:12])
        old_sha = old["implementation"]["example_source_sha256"]
        old_commit = old["implementation"]["frozen_source_commit"]
        manifest = replace_sources(old, old_sha, source_sha, old_commit, source_commit)
        assert isinstance(manifest, dict)
        prefix = arm == "expanded"
        schedule = ("c1200_prefix256_then_striped8_c1200_then_complete_c1800" if prefix
                    else "native_rotating_state_order")
        prefix_columns = 1200
        manifest["implementation"]["algorithm_flags"]["target_search_schedule"] = schedule
        manifest["implementation"]["algorithm_flags"]["target_search_prefix_columns"] = prefix_columns
        manifest["implementation"]["algorithm_flags"]["target_search_lead_states"] = 256 if prefix else 0
        manifest["implementation"]["algorithm_flags"]["target_search_stripes"] = 8 if prefix else 1
        manifest["point_decomposition"]["target_search_schedule"] = schedule
        manifest["point_decomposition"]["target_search_prefix_columns"] = prefix_columns
        manifest["point_decomposition"]["target_search_lead_states"] = 256 if prefix else 0
        manifest["point_decomposition"]["target_search_stripes"] = 8 if prefix else 1
        manifest["point_decomposition"]["target_query_policy"] = (
            "point_seeded_c1200_prefix256_then_striped8_c1200_then_complete_c1800_"
            "with_first_verified_relation_and_worker_drain" if prefix else
            "point_seeded_native_rotating_state_order_"
            "with_first_verified_relation_and_worker_drain")
        candidate_id = old_id.rsplit("h", 1)[0] + "h" + digest(canonical(manifest))[:12]
        manifest["candidate_id"] = candidate_id
        destination = HERE / "identities/candidates" / (candidate_id + ".json")
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(canonical(manifest) + b"\n")
        arms[arm] = {
            "candidate_id": candidate_id,
            "candidate_manifest": str(destination.relative_to(ROOT)),
            "columns": manifest["factor_base"]["effective_relation_columns"],
            "factor_base_digest": manifest["factor_base"]["enumerated_set_digest"],
            "factor_base_points": manifest["factor_base"]["actual_usable_points_B"],
            "prepare_rank": manifest["factor_base"]["effective_relation_columns"],
            "root_table_slots": manifest["implementation"]["algorithm_flags"]["root_table_slots"],
            "target_search_schedule": schedule,
            "target_search_prefix_columns": prefix_columns,
            "target_search_lead_states": 256 if prefix else 0,
            "target_search_stripes": 8 if prefix else 1,
        }
    assert arms["baseline"]["columns"] == 1200
    assert arms["expanded"]["columns"] == 1800
    plan = {
        "kind": "s3_n53_c1800stripe_T110_diagnostic_freeze_v1",
        "parent_freeze_sha256": digest(PARENT.read_bytes()),
        "source_commit": source_commit,
        "source_sha256": source_sha,
        "binary_sha256": binary_sha,
        "runner_sha256": digest(RUNNER.read_bytes()),
        "arms": arms,
        "target_labels": ["T110"],
        "repetitions_per_arm_per_target": 1,
        "run_order": "baseline_then_expanded",
        "relation_seed": 20260928,
        "relation_workers": 14,
        "target_workers": 1,
        "timeout_seconds_per_run": 180,
        "memory_envelope_bytes": 34359738368,
        "prior_used_target_control": True,
        "host_exclusive_partition_receipt": None,
        "controlled_speedup": None,
    }
    (HERE / "DIAGNOSTIC_FREEZE.json").write_text(json.dumps(plan, indent=2, sort_keys=True) + "\n")
    print(json.dumps({arm: item["candidate_id"] for arm, item in arms.items()}, sort_keys=True))


if __name__ == "__main__":
    main()
