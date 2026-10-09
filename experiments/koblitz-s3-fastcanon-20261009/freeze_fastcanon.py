"""Freeze IC1 identities for bitset Frobenius canonicalization."""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import subprocess


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
GROUPLIFT = ROOT / "experiments/koblitz-s3-grouplift-20261009/FREEZE.json"
GROUPLIFT_DIR = GROUPLIFT.parent
SOURCE = ROOT / "examples/s3_partial_prepared_fastcanon_targetspan.rs"
BINARY = Path("/private/tmp/s3-prepared-target-20261009/release/examples/s3_partial_prepared_fastcanon_targetspan")


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
    grouplift_bytes = GROUPLIFT.read_bytes()
    grouplift = json.loads(grouplift_bytes)
    source_sha = sha(SOURCE.read_bytes())
    source_tree = subprocess.check_output(
        ["git", "-C", str(ROOT), "rev-parse", "HEAD:src"], text=True).strip()
    assert grouplift["parent_plan_sha256"] == sha(parent_bytes)
    assert grouplift["source_tree_git_oid"] == source_tree
    candidates = {}
    for n in (41, 53):
        old_id = grouplift["curves"][str(n)]["grouplift_candidate_id"]
        old_path = GROUPLIFT_DIR / "identities/candidates" / f"{old_id}.json"
        old = json.loads(old_path.read_bytes())
        assert old.pop("candidate_id") == old_id
        assert sha(canonical(old))[:12] == old_id[-12:]
        assert old["implementation"]["example_source_sha256"] == grouplift["source_sha256"]
        record = copy.deepcopy(old)
        for section in ("factor_base", "relation_collection", "relation_linear_algebra",
                        "target_descent"):
            record[section]["source_sha256"] = source_sha
        point = record["point_decomposition"]
        point["implementation_source_sha256"] = source_sha
        point["implementation"] = (
            "parallel_sharded_s3_index_prelifted_group_sum_roots_fastcanon")
        point["frobenius_canonicalization"] = (
            "bitset_longest_zero_run_then_exact_min_rotation_with_earliest_shift_tie")
        record["implementation"]["example_source_sha256"] = source_sha
        record["implementation"]["frozen_source_commit"] = None
        record["implementation"]["algorithm_flags"]["bitset_frobenius_canonicalization"] = True
        ident = (f"IC1N{n}Ckb1fb{parent['curves'][str(n)]['factor_base_points']}"
                 f"PDP4rootRCsampleLAgaussTDdirectISO0h{sha(canonical(record))[:12]}")
        write_exact(HERE / "identities/candidates" / f"{ident}.json",
                    canonical({**record, "candidate_id": ident}) + b"\n")
        candidates[str(n)] = {"fastcanon_candidate_id": ident,
                              "grouplift_candidate_id": old_id,
                              "curve_id": parent["curves"][str(n)]["curve_id"]}
    selected = [{"n": target["n"], "label": target["label"],
                 "workload_id": target["workload_id"],
                 "point_sha256": target["point_sha256"]}
                for target in parent["targets"] if int(target["label"][1:]) <= 24]
    assert len(selected) == 48
    freeze = {"kind": "s3_fastcanon_matched_curve_pilot_v1",
              "parent_plan_sha256": sha(parent_bytes),
              "grouplift_plan_sha256": sha(grouplift_bytes),
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
              "targets": selected,
              "controlled_speedup": None}
    write_exact(HERE / "FREEZE.json", canonical(freeze) + b"\n")
    print(json.dumps({"source_sha256": source_sha,
                      "binary_sha256": freeze["local_binary_sha256"],
                      "candidates": candidates,
                      "plan_sha256": sha((HERE / "FREEZE.json").read_bytes())}, sort_keys=True))


if __name__ == "__main__":
    main()
