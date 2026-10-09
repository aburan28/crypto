"""Freeze the persistent 2,048-state target tile and paired reference IDs."""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
TILED = ROOT / "experiments/koblitz-s3-tiled-target-20261009/FREEZE.json"
PARALLEL = ROOT / "experiments/koblitz-s3-parallel-target-20261009/FREEZE.json"
SOURCE = ROOT / "examples/s3_partial_prepared_persistent_tiled_targetspan.rs"
BINARY = Path("/private/tmp/s3-prepared-target-20261009/release/examples/s3_partial_prepared_persistent_tiled_targetspan")


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
    tiled_bytes = TILED.read_bytes()
    tiled = json.loads(tiled_bytes)
    parallel_bytes = PARALLEL.read_bytes()
    parallel = json.loads(parallel_bytes)
    assert tiled["parent_plan_sha256"] == parallel["parent_plan_sha256"] == sha(parent_bytes)
    tiled_id = tiled["tiled_candidate_ids_by_tile_states"]["2048"]
    tiled_path = TILED.parent / "identities/candidates" / f"{tiled_id}.json"
    base_record = json.loads(tiled_path.read_bytes())
    assert base_record.pop("candidate_id") == tiled_id
    assert sha(canonical(base_record))[:12] == tiled_id[-12:]
    source_sha = sha(SOURCE.read_bytes())
    record = copy.deepcopy(base_record)
    for section in ("factor_base", "relation_collection", "relation_linear_algebra",
                    "target_descent"):
        record[section]["source_sha256"] = source_sha
    point = record["point_decomposition"]
    point["implementation_source_sha256"] = source_sha
    point["implementation"] = (
        "parallel_sharded_s3_index_prelifted_group_sum_roots_persistent_tiled_target")
    point["target_query_policy"] = (
        "fixed_serial_order_tiles_persistent_workers_earliest_verified_relation_and_worker_drain")
    point["persistent_target_workers"] = True
    record["target_descent"]["persistent_target_workers"] = True
    record["implementation"]["example_source_sha256"] = source_sha
    record["implementation"]["frozen_source_commit"] = None
    record["implementation"]["algorithm_flags"]["persistent_target_workers"] = True
    persistent_id = f"IC1N53Ckb1fb25864PDP4rootRCsampleLAgaussTDdirectISO0h{sha(canonical(record))[:12]}"
    write_exact(HERE / "identities/candidates" / f"{persistent_id}.json",
                canonical({**record, "candidate_id": persistent_id}) + b"\n")
    selected = [{"label": target["label"], "workload_id": target["workload_id"],
                 "point_sha256": target["point_sha256"]}
                for target in parent["targets"]
                if target["n"] == 53 and 25 <= int(target["label"][1:]) <= 32]
    assert len(selected) == 8
    freeze = {"kind": "s3_persistent_tiled_paired_v1",
              "parent_plan_sha256": sha(parent_bytes),
              "tiled_plan_sha256": sha(tiled_bytes),
              "parallel_plan_sha256": sha(parallel_bytes),
              "source_sha256": source_sha,
              "local_binary_sha256": sha(BINARY.read_bytes()),
              "curve_id": parent["curves"]["53"]["curve_id"],
              "candidate_ids_by_arm": {
                  "parallel": parallel["curves"]["53"]["parallel_candidate_id"],
                  "tiled": tiled_id, "persistent": persistent_id},
              "target_tile_states": 2048,
              "target_query_workers": 14,
              "relation_workers": 14,
              "prepare_rank": 243,
              "relation_seed": parent["relation_seed"],
              "repetitions_per_arm_per_target": 2,
              "targets": selected,
              "controlled_speedup": None}
    write_exact(HERE / "FREEZE.json", canonical(freeze) + b"\n")
    print(json.dumps({"source_sha256": source_sha,
                      "binary_sha256": freeze["local_binary_sha256"],
                      "candidate_id": persistent_id,
                      "plan_sha256": sha((HERE / "FREEZE.json").read_bytes())}, sort_keys=True))


if __name__ == "__main__":
    main()
