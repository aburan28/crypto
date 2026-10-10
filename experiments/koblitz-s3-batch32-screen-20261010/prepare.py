"""Freeze a source-matched batch-eight versus batch-32 S3 target screen."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BASE_NAME = "s3_partial_prepared_fullrank_prefetch_projective_c1800prefix_targetspan"
NEW_NAME = BASE_NAME + "_batch32"
PARENTS = {
    41: (ROOT / "experiments/koblitz-s3-targetkernel-n41-20261010", "improved"),
    53: (ROOT / "experiments/koblitz-s3-c1800-matchedsource-n53-20261010", "expanded"),
}


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def canonical(value: dict) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode()


def replace_source(value, old: str, new: str):
    if isinstance(value, dict):
        return {key: replace_source(item, old, new) for key, item in value.items()}
    if isinstance(value, list):
        return [replace_source(item, old, new) for item in value]
    return new if value == old else value


def main(binaries: Path) -> None:
    if (HERE / "FREEZE.json").exists():
        raise FileExistsError(HERE / "FREEZE.json")
    source_hashes = {name: sha((ROOT / "examples" / f"{name}.rs").read_bytes())
                     for name in (BASE_NAME, NEW_NAME)}
    binary_hashes = {name: sha((binaries / name).read_bytes())
                     for name in (BASE_NAME, NEW_NAME)}
    assert source_hashes[BASE_NAME] == "79f77fef1857c9aed10b8d8913ea796e501f9b0afd615345746c453f8d4907b2"
    assert source_hashes[NEW_NAME] == "aa575e1b9540e1ea8cb77a49aeaa1b4ec0ef8a36bc5ea7e8a48b5d64a446bd18"
    assert binary_hashes[BASE_NAME] == "d10625fe3aba0e7b61c44f001f119c1bb3b7f4fdde4ee558df1501f1f1ab1164"
    assert binary_hashes[NEW_NAME] == "180493c3aaa371784890b3ec1bff8b9365feeba9772e277e466e1474a8100175"
    source_commit = subprocess.check_output(
        ["git", "-C", str(ROOT), "log", "-1", "--format=%H", "--",
         f"examples/{NEW_NAME}.rs"], text=True).strip()
    assert source_commit == "15192dfebca31a0f6ec1dc980deffbc511f3280b"
    tree_oid = subprocess.check_output(
        ["git", "-C", str(ROOT), "rev-parse", "HEAD:src"], text=True).strip()
    assert tree_oid == "0a574800502e67882d8f8f88a846c1e80ee795fd"
    freeze = {
        "kind": "s3_batch32_source_matched_screen_v1",
        "comparison_variable": "target batch inversion and lookup prefetch width, eight versus 32 states",
        "measured_interval": "one public target, target-dependent query through independently checked scalar",
        "controlled_speedup": None,
        "host_exclusive_partition_receipt": None,
        "development_gate_online_ratio_at_least": "1.10",
        "per_run_timeout_seconds": 180,
        "observed_peak_rss_envelope_bytes": 34359738368,
        "relation_seed": 20260928,
        "relation_workers": 14,
        "target_workers": 1,
        "source_tree_git_oid": tree_oid,
        "cargo_lock_sha256": sha((ROOT / "Cargo.lock").read_bytes()),
        "preparer_sha256": sha(Path(__file__).read_bytes()),
        "runner_sha256": sha((HERE / "run.py").read_bytes()),
        "analyzer_sha256": sha((HERE / "analyze.py").read_bytes()),
        "replay_sha256": sha((HERE / "replay_sage.py").read_bytes()),
        "source_sha256": source_hashes,
        "binary_sha256": binary_hashes,
        "parent_freeze_sha256": {},
        "curves": {},
        "selected_targets": [],
        "run_order": [],
    }
    identities = HERE / "identities/candidates"
    identities.mkdir(parents=True, exist_ok=True)
    for n, (parent, parent_arm) in PARENTS.items():
        parent_freeze_path = parent / "FREEZE.json"
        parent_freeze = json.loads(parent_freeze_path.read_bytes())
        parent_arm_record = parent_freeze["arms"][parent_arm]
        assert parent_arm_record["binary_name"] == BASE_NAME
        assert parent_arm_record["binary_sha256"] == binary_hashes[BASE_NAME]
        assert parent_arm_record["source_sha256"] == source_hashes[BASE_NAME]
        freeze["parent_freeze_sha256"][str(n)] = sha(parent_freeze_path.read_bytes())
        baseline_path = parent / parent_arm_record["candidate_manifest_relative"]
        baseline = json.loads(baseline_path.read_bytes())
        assert baseline["candidate_id"] == parent_arm_record["candidate_id"]
        candidate = replace_source(baseline, source_hashes[BASE_NAME], source_hashes[NEW_NAME])
        del candidate["candidate_id"]
        candidate["implementation"]["frozen_source_commit"] = source_commit
        flags = candidate["implementation"]["algorithm_flags"]
        flags["target_query_batch_width"] = 32
        flags["target_query_inversion_policy"] = "one_normal_basis_inverse_per_32_regular_state_pair_products"
        flags["target_lookup_prefetch_policy"] = "prefetch_all_candidate_home_slots_after_32_state_root_batch_then_original_lookup_order"
        pdp = candidate["point_decomposition"]
        pdp["implementation"] = pdp["implementation"].replace("batch8", "batch32")
        pdp["target_query_batch_width"] = 32
        pdp["target_query_inversion_policy"] = flags["target_query_inversion_policy"]
        pdp["target_lookup_prefetch_policy"] = flags["target_lookup_prefetch_policy"]
        prefix = parent_arm_record["candidate_id"].rsplit("h", 1)[0]
        candidate_id = prefix + "h" + sha(canonical(candidate))[:12]
        candidate["candidate_id"] = candidate_id
        manifest_path = identities / f"{candidate_id}.json"
        manifest_path.write_bytes(canonical(candidate) + b"\n")
        freeze["curves"][str(n)] = {
            "columns": parent_arm_record["columns"],
            "factor_base_points": parent_arm_record["factor_base_points"],
            "factor_base_digest": parent_arm_record["factor_base_digest"],
            "prepare_rank": parent_arm_record["prepare_rank"],
            "root_table_slots": parent_arm_record["root_table_slots"],
            "target_search_schedule": parent_arm_record.get(
                "target_search_schedule", "native_rotating_state_order"),
            "target_search_prefix_columns": parent_arm_record.get(
                "target_search_prefix_columns", parent_arm_record["columns"]),
            "arms": {
                "batch8": {"binary_name": BASE_NAME,
                           "candidate_id": parent_arm_record["candidate_id"],
                           "candidate_manifest_sha256": sha(baseline_path.read_bytes())},
                "batch32": {"binary_name": NEW_NAME, "candidate_id": candidate_id,
                            "candidate_manifest_relative": str(manifest_path.relative_to(HERE)),
                            "candidate_manifest_sha256": sha(manifest_path.read_bytes())},
            },
        }
        for target in parent_freeze["selected_targets"][:8]:
            point_path = (parent.parent / "koblitz-s3-crosscurve-20261009" /
                          target["point_file_relative"] if n == 41 else
                          ROOT / target["point_file_relative"])
            assert sha(point_path.read_bytes()) == target["point_sha256"]
            freeze["selected_targets"].append({
                "n": n, "label": target["label"],
                "point_file_relative": str(point_path.relative_to(ROOT)),
                "point_sha256": target["point_sha256"],
                "public_point": target["public_point"],
                "fixture_scalar": target["fixture_scalar"],
                "workload_id": target["workload_id"],
                "workload_identity_sha256": target["workload_identity_sha256"],
            })
    by_curve = {n: [t for t in freeze["selected_targets"] if t["n"] == n]
                for n in (41, 53)}
    for index in range(8):
        for n in (41, 53):
            target = by_curve[n][index]
            order = ("batch8", "batch32") if (index + n) % 2 else ("batch32", "batch8")
            for arm in order:
                freeze["run_order"].append({"n": n, "label": target["label"],
                                            "arm": arm, "repetition": 1})
    freeze["planned_run_count"] = len(freeze["run_order"])
    assert freeze["planned_run_count"] == 32
    (HERE / "FREEZE.json").write_text(json.dumps(freeze, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"freeze_sha256": sha((HERE / "FREEZE.json").read_bytes()),
                      "planned_runs": 32,
                      "candidate_ids": {n: freeze["curves"][str(n)]["arms"]["batch32"]["candidate_id"]
                                        for n in (41, 53)}}, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binaries", type=Path, required=True)
    args = parser.parse_args()
    main(args.binaries.resolve())
