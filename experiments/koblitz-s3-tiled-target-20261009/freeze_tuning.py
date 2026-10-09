"""Freeze IC1 identities for three deterministic S3 target tile sizes."""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PARALLEL = ROOT / "experiments/koblitz-s3-parallel-target-20261009/FREEZE.json"
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
SOURCE = ROOT / "examples/s3_partial_prepared_tiled_targetspan.rs"
BINARY = Path("/private/tmp/s3-prepared-target-20261009/release/examples/s3_partial_prepared_tiled_targetspan")
TILES = (2048, 4096, 8192)


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
    old_id = parallel["curves"]["53"]["parallel_candidate_id"]
    old_path = PARALLEL.parent / "identities/candidates" / f"{old_id}.json"
    old = json.loads(old_path.read_bytes())
    assert old.pop("candidate_id") == old_id
    assert sha(canonical(old))[:12] == old_id[-12:]
    source_sha = sha(SOURCE.read_bytes())
    identities = {}
    policy = ("fixed_serial_order_tiles_parallel_subranges_earliest_verified_"
              "relation_and_worker_drain")
    for tile_states in TILES:
        record = copy.deepcopy(old)
        for section in ("factor_base", "relation_collection", "relation_linear_algebra",
                        "target_descent"):
            record[section]["source_sha256"] = source_sha
        point = record["point_decomposition"]
        point["implementation_source_sha256"] = source_sha
        point["implementation"] = (
            "parallel_sharded_s3_index_prelifted_group_sum_roots_tiled_target")
        point["target_query_policy"] = policy
        point["target_tile_states"] = tile_states
        record["target_descent"]["target_tile_states"] = tile_states
        record["implementation"]["example_source_sha256"] = source_sha
        record["implementation"]["frozen_source_commit"] = None
        record["implementation"]["algorithm_flags"]["cooperative_target_cancellation"] = False
        record["implementation"]["algorithm_flags"]["target_tile_states"] = tile_states
        candidate = f"IC1N53Ckb1fb25864PDP4rootRCsampleLAgaussTDdirectISO0h{sha(canonical(record))[:12]}"
        write_exact(HERE / "identities/candidates" / f"{candidate}.json",
                    canonical({**record, "candidate_id": candidate}) + b"\n")
        identities[str(tile_states)] = candidate
    selected = [{"n": 53, "label": target["label"],
                 "workload_id": target["workload_id"],
                 "point_sha256": target["point_sha256"]}
                for target in parent["targets"]
                if target["n"] == 53 and 25 <= int(target["label"][1:]) <= 32]
    assert len(selected) == 8
    freeze = {"kind": "s3_tiled_target_size_tuning_v1",
              "parent_plan_sha256": sha(parent_bytes),
              "parallel_plan_sha256": sha(parallel_bytes),
              "source_sha256": source_sha,
              "local_binary_sha256": sha(BINARY.read_bytes()),
              "curve_id": parent["curves"]["53"]["curve_id"],
              "reference_candidate_id": old_id,
              "tiled_candidate_ids_by_tile_states": identities,
              "target_query_workers": 14,
              "relation_workers": 14,
              "prepare_rank": 243,
              "relation_seed": parent["relation_seed"],
              "target_count": 8,
              "repetitions_per_arm_per_target": 2,
              "targets": selected,
              "controlled_speedup": None}
    write_exact(HERE / "FREEZE.json", canonical(freeze) + b"\n")
    print(json.dumps({"source_sha256": source_sha,
                      "binary_sha256": freeze["local_binary_sha256"],
                      "plan_sha256": sha((HERE / "FREEZE.json").read_bytes()),
                      "candidates": identities}, sort_keys=True))


if __name__ == "__main__":
    main()
