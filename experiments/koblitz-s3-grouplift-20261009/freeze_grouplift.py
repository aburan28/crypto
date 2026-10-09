"""Freeze IC1 identities for the prelifted group-sum S3 index."""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import subprocess


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
SOURCE = ROOT / "examples/s3_partial_prepared_grouplift_targetspan.rs"
BINARY = Path("/private/tmp/s3-prepared-target-20261009/release/examples/s3_partial_prepared_grouplift_targetspan")
PRIOR = ROOT / "experiments/koblitz-s3-partial-prepared-20261009"
PRIOR_MANIFEST_DIR = {
    41: PRIOR / "identities/candidates",
    53: PRIOR / "rank243/identities/candidates",
}


def canonical(value: dict) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode("utf-8")


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def write_exact(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        assert path.read_bytes() == data, f"frozen artifact changed: {path}"
    else:
        path.write_bytes(data)


def main() -> None:
    parent_bytes = PARENT.read_bytes()
    parent = json.loads(parent_bytes)
    source_sha = sha(SOURCE.read_bytes())
    source_tree = subprocess.check_output(
        ["git", "-C", str(ROOT), "rev-parse", "HEAD:src"], text=True).strip()
    candidates = {}
    for n in (41, 53):
        base_id = parent["curves"][str(n)]["rank243_candidate_id"]
        base_path = PRIOR_MANIFEST_DIR[n] / f"{base_id}.json"
        base = json.loads(base_path.read_bytes())
        assert base.pop("candidate_id") == base_id
        assert sha(canonical(base))[:12] == base_id[-12:]
        assert base["implementation"]["example_source_sha256"] == parent["partial_source_sha256"]
        assert base["implementation"]["rust_library_source_tree_git_oid"] == source_tree
        assert base["relation_collection"]["pre_target_rank_threshold"] == 243
        record = copy.deepcopy(base)
        for section in ("factor_base", "relation_collection", "relation_linear_algebra",
                        "target_descent"):
            record[section]["source_sha256"] = source_sha
        record["point_decomposition"]["implementation_source_sha256"] = source_sha
        record["point_decomposition"]["implementation"] = (
            "parallel_sharded_s3_index_prelifted_group_sum_roots")
        record["point_decomposition"]["index_state_y_policy"] = (
            "prelift_both_s3_x_roots_to_curve_points_before_target_input")
        record["point_decomposition"]["root_query_order"] = (
            "indexed_point_plus_minus_target_then_original_order_witness_checks")
        record["implementation"]["example_source_sha256"] = source_sha
        record["implementation"]["frozen_source_commit"] = None
        record["implementation"]["algorithm_flags"]["prelifted_group_s3_roots"] = True
        record["implementation"]["algorithm_flags"]["shared_inversion_across_two_index_points"] = True
        ident = (f"IC1N{n}Ckb1fb{parent['curves'][str(n)]['factor_base_points']}"
                 f"PDP4rootRCsampleLAgaussTDdirectISO0h{sha(canonical(record))[:12]}")
        write_exact(HERE / "identities/candidates" / f"{ident}.json",
                    canonical({**record, "candidate_id": ident}) + b"\n")
        candidates[str(n)] = {"grouplift_candidate_id": ident,
                              "prepared_candidate_id": base_id,
                              "curve_id": parent["curves"][str(n)]["curve_id"]}
    selected = [{"n": target["n"], "label": target["label"],
                 "workload_id": target["workload_id"],
                 "point_sha256": target["point_sha256"]}
                for target in parent["targets"] if int(target["label"][1:]) <= 24]
    assert len(selected) == 48
    freeze = {"kind": "s3_grouplift_matched_curve_pilot_v1",
              "parent_plan_sha256": sha(parent_bytes),
              "source_sha256": source_sha,
              "local_binary_sha256": sha(BINARY.read_bytes()),
              "source_tree_git_oid": source_tree,
              "relation_seed": parent["relation_seed"],
              "workers": parent["workers"],
              "prepare_rank": 243,
              "timeout_seconds_per_run": parent["timeout_seconds_per_run"],
              "target_count_per_curve": 24,
              "repetitions_per_arm_per_target": 2,
              "order": "ABBA_or_BAAB_rotated_by_index_and_curve",
              "curves": candidates,
              "targets": selected}
    write_exact(HERE / "FREEZE.json", canonical(freeze) + b"\n")
    print(json.dumps({"source_sha256": source_sha,
                      "binary_sha256": freeze["local_binary_sha256"],
                      "candidates": candidates,
                      "plan_sha256": sha((HERE / "FREEZE.json").read_bytes())}, sort_keys=True))


if __name__ == "__main__":
    main()
