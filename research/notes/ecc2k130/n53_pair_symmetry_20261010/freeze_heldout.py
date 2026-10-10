#!/usr/bin/env python3
"""Bind the once-generated Q to six one-target paired workloads."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from freeze import canonical, save_new_or_equal, sha


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
PARENT = ROOT / "research/notes/ecc2k130/n53_fixed_base_rank_restarts_20261010"


def main() -> None:
    prep = json.loads((HERE / "FROZEN_PREP.json").read_text())
    generation = json.loads((HERE / "HELDOUT_GENERATION.json").read_text())
    receipt = json.loads((HERE / "inputs/heldout_generation_receipt.json").read_text())
    point_path = HERE / "inputs/heldout_q.jsonl"
    fixture_path = HERE / "inputs/heldout_fixture_verifier_only.jsonl"
    point = json.loads(point_path.read_text())
    fixture = json.loads(fixture_path.read_text())
    assert receipt["status"] == "ELIGIBLE" and not receipt["prior_committed_matches"]
    assert receipt["manifest_sha256"] == sha(HERE / "HELDOUT_GENERATION.json")
    assert receipt["frozen_prep_sha256"] == sha(HERE / "FROZEN_PREP.json")
    assert receipt["point"] == point == fixture["published_q"]
    assert receipt["point_sha256"] == sha(point_path)
    assert receipt["fixture_sha256"] == sha(fixture_path)
    assert generation["rank_seeds"] == prep["rank_seeds"]
    assert generation["rho_seeds"] == prep["rho_seeds"]
    workload_ids = {}
    for rank_seed, rho_seed in zip(prep["rank_seeds"], prep["rho_seeds"], strict=True):
        record = {
            "curve_id": prep["curve_id"], "subgroup_order": "21044858204113",
            "target": point, "target_kind": "public_hash_to_curve_cofactor",
            "public_hash_seed": prep["public_hash_seed"],
            "public_hash_counter": fixture["public_hash_counter"],
            "target_point_sha256": sha(point_path),
            "input_law": "one unseen public point Q, scalar withheld from both timed arms",
            "rank_seed": rank_seed, "rho_seed": rho_seed, "target_count": 1,
            "cache_state": "cold",
        }
        workload_id = "W" + hashlib.sha256(canonical(record)).hexdigest()[:12]
        workload_ids[str(rank_seed)] = workload_id
        save_new_or_equal(HERE / f"heldout_workload_{rank_seed}.json",
                          {"workload_id": workload_id, "record": record})
    frozen = {
        "schema": "n53-pair-symmetry-heldout-freeze-v1",
        "generation_rule_sha256": sha(HERE / "HELDOUT_GENERATION.json"),
        "generation_receipt_sha256": sha(HERE / "inputs/heldout_generation_receipt.json"),
        "frozen_prep_sha256": sha(HERE / "FROZEN_PREP.json"),
        "public_point": point, "public_point_sha256": sha(point_path),
        "fixture_verifier_only_sha256": sha(fixture_path),
        "public_hash_counter": fixture["public_hash_counter"],
        "ic_candidate_ids": prep["candidate_ids"],
        "ic_source_sha256": prep["source_sha256"],
        "ic_binary_sha256": prep["ic_binary_sha256"],
        "rho_source_sha256": prep["rho_source_sha256"],
        "rho_binary_sha256": prep["rho_binary_sha256"],
        "runner_sha256": sha(HERE / "run_heldout.py"),
        "rho_verifier_sha256": sha(PARENT / "verify_rho_public.py"),
        "ic_verifier_sha256": sha(PARENT / "verify_public_target.py"),
        "rank_seeds": prep["rank_seeds"], "rho_seeds": prep["rho_seeds"],
        "workload_ids": workload_ids, "resource_envelope": prep["resource_envelope"],
    }
    save_new_or_equal(HERE / "HELDOUT_FROZEN.json", frozen)
    print(json.dumps({"point": point, "workload_ids": workload_ids}, sort_keys=True))


if __name__ == "__main__":
    main()
