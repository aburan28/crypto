"""Freeze batched S3 group-inversion candidates on both binary curves."""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import subprocess


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
PARALLEL = ROOT / "experiments/koblitz-s3-parallel-target-20261009/FREEZE.json"
SOURCE = ROOT / "examples/s3_partial_prepared_batched_group_targetspan.rs"
BINARY = Path("/private/tmp/s3-prepared-target-20261009/release/examples/s3_partial_prepared_batched_group_targetspan")


def canonical(value: dict) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode()


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
    parallel_bytes = PARALLEL.read_bytes()
    parallel = json.loads(parallel_bytes)
    assert parallel["parent_plan_sha256"] == sha(parent_bytes)
    source_tree = subprocess.check_output(
        ["git", "-C", str(ROOT), "rev-parse", "HEAD:src"], text=True).strip()
    assert parallel["source_tree_git_oid"] == source_tree
    source_sha = sha(SOURCE.read_bytes())
    curves = {}
    for n in (41, 53):
        old_id = parallel["curves"][str(n)]["parallel_candidate_id"]
        old_path = PARALLEL.parent / "identities/candidates" / f"{old_id}.json"
        record = json.loads(old_path.read_bytes())
        assert record.pop("candidate_id") == old_id
        assert sha(canonical(record))[:12] == old_id[-12:]
        record = copy.deepcopy(record)
        for section in ("factor_base", "relation_collection", "relation_linear_algebra",
                        "target_descent"):
            record[section]["source_sha256"] = source_sha
        point = record["point_decomposition"]
        point["implementation_source_sha256"] = source_sha
        point["implementation"] = (
            "parallel_sharded_s3_index_prelifted_group_sum_roots_batched_inversion_target")
        point["target_query_policy"] = (
            "independent_rotating_starts_first_verified_relation_batched_group_inversion_cooperative_cancel_and_worker_drain")
        point["target_group_inversion_batch_states"] = 16
        record["implementation"]["example_source_sha256"] = source_sha
        record["implementation"]["frozen_source_commit"] = None
        record["implementation"]["algorithm_flags"]["target_group_inversion_batch_states"] = 16
        ident = (f"IC1N{n}Ckb1fb{parent['curves'][str(n)]['factor_base_points']}"
                 f"PDP4rootRCsampleLAgaussTDdirectISO0h{sha(canonical(record))[:12]}")
        write_exact(HERE / "identities/candidates" / f"{ident}.json",
                    canonical({**record, "candidate_id": ident}) + b"\n")
        curves[str(n)] = {
            "curve_id": parent["curves"][str(n)]["curve_id"],
            "parallel_candidate_id": old_id,
            "batched_candidate_id": ident,
            "target_query_workers": parallel["target_workers_by_curve"][str(n)]}
    selected = [{"n": target["n"], "label": target["label"],
                 "workload_id": target["workload_id"],
                 "point_sha256": target["point_sha256"]}
                for target in parent["targets"]
                if 25 <= int(target["label"][1:]) <= 32]
    assert len(selected) == 16
    freeze = {"kind": "s3_batched_group_inversion_paired_v1",
              "parent_plan_sha256": sha(parent_bytes),
              "parallel_plan_sha256": sha(parallel_bytes),
              "source_sha256": source_sha,
              "local_binary_sha256": sha(BINARY.read_bytes()),
              "source_tree_git_oid": source_tree,
              "relation_seed": parent["relation_seed"],
              "workers": parent["workers"],
              "prepare_rank": 243,
              "target_group_inversion_batch_states": 16,
              "timeout_seconds_per_run": parent["timeout_seconds_per_run"],
              "repetitions_per_arm_per_target": 2,
              "curves": curves,
              "targets": selected,
              "controlled_speedup": None}
    write_exact(HERE / "FREEZE.json", canonical(freeze) + b"\n")
    print(json.dumps({"source_sha256": source_sha,
                      "binary_sha256": freeze["local_binary_sha256"],
                      "curves": curves,
                      "plan_sha256": sha((HERE / "FREEZE.json").read_bytes())}, sort_keys=True))


if __name__ == "__main__":
    main()
