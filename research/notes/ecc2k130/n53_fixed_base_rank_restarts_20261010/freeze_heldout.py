#!/usr/bin/env python3
"""Bind the generated held-out Q and six paired seed workloads before timing."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from freeze import digest, write_frozen


HERE = Path(__file__).resolve().parent
RANK_SEEDS = tuple(range(531053, 531059))
RHO_SEEDS = tuple(range(530153, 530159))


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    rule_path = HERE / "HELDOUT_GENERATION.json"
    rule = json.loads(rule_path.read_text())
    amendment_path = HERE / "HELDOUT_GENERATION_AMENDMENT.json"
    amendment = json.loads(amendment_path.read_text())
    assert amendment["original_generation_rule_sha256"] == sha(rule_path)
    rule.update(amendment["overrides"])
    receipt_path = HERE / "inputs/heldout_generation_receipt.json"
    receipt = json.loads(receipt_path.read_text())
    assert receipt["status"] == "ELIGIBLE" and not receipt["prior_committed_matches"]
    assert receipt["manifest_sha256"] == sha(rule_path)
    assert receipt["amendment_sha256"] == sha(amendment_path)
    assert rule["rank_seeds"] == list(RANK_SEEDS)
    assert rule["rho_seeds"] == list(RHO_SEEDS)
    assert rule["selected_rank_probe_cap"] == 400000
    q_path = HERE / "inputs/heldout_q.jsonl"
    fixture_path = HERE / "inputs/heldout_fixture_verifier_only.jsonl"
    point = json.loads(q_path.read_text())
    fixture = json.loads(fixture_path.read_text())
    assert point == receipt["point"] == fixture["published_q"]
    assert sha(q_path) == receipt["point_sha256"]
    assert sha(fixture_path) == receipt["fixture_sha256"]
    pilot = json.loads((HERE / "PILOT_ANALYSIS.json").read_text())
    assert pilot["selected_cap_for_heldout"] == "cap400000"
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    assert rule["ic_candidate_ids"] == {
        label: frozen["candidate_ids"][label] for label in ("control", "cap400000")
    }
    workload_ids = {}
    for rank_seed, rho_seed in zip(RANK_SEEDS, RHO_SEEDS, strict=True):
        record = {
            "curve_id": rule["curve_id"],
            "subgroup_order": "21044858204113",
            "target": point,
            "target_kind": "public_hash_to_curve_cofactor",
            "public_hash_seed": rule["public_hash_seed"],
            "public_hash_counter": fixture["public_hash_counter"],
            "target_point_sha256": sha(q_path),
            "input_law": "one unseen public point Q, scalar withheld from both timed arms",
            "rank_seed": rank_seed,
            "rho_seed": rho_seed,
            "target_count": 1,
            "cache_state": "cold",
        }
        workload_id = "W" + digest(record)[:12]
        workload_ids[str(rank_seed)] = workload_id
        write_frozen(HERE / f"heldout_workload_{rank_seed}.json",
                     {"workload_id": workload_id, "record": record})
    write_frozen(HERE / "HELDOUT_FROZEN.json", {
        "schema": "n53-rank-restart-heldout-freeze-v1",
        "generation_rule_sha256": sha(rule_path),
        "generation_amendment_sha256": sha(amendment_path),
        "generation_receipt_sha256": sha(receipt_path),
        "fixture_verifier_only_sha256": sha(fixture_path),
        "public_point_sha256": sha(q_path),
        "public_point": point,
        "public_hash_seed": rule["public_hash_seed"],
        "public_hash_counter": fixture["public_hash_counter"],
        "selected_cap": 400000,
        "ic_candidate_ids": rule["ic_candidate_ids"],
        "ic_source_sha256": frozen["source_sha256"],
        "ic_binary_sha256": frozen["binary_sha256"],
        "rho_source_sha256": rule["rho_online_source_sha256"],
        "rho_binary_sha256": rule["rho_online_binary_sha256"],
        "rank_seeds": list(RANK_SEEDS),
        "rho_seeds": list(RHO_SEEDS),
        "workload_ids": workload_ids,
        "resource_envelope": rule["per_arm_resource_envelope"],
    })
    print(json.dumps({"point": point, "workload_ids": workload_ids}, sort_keys=True))


if __name__ == "__main__":
    main()
