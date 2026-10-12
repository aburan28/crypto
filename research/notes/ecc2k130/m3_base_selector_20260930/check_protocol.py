#!/usr/bin/env python3
"""Read-only pre-outcome checks; does not generate candidate bases or Q."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
OLD = ROOT / "research/notes/ecc2k130/m3_four_policy_20260930"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check() -> dict:
    config = json.loads((HERE / "CONFIG.json").read_text())
    old = json.loads((OLD / "CONFIG.json").read_text())
    prior_result = OLD / "evidence_run_36722040881/result.json"
    assert config["schema"] == "ecc2k130-degree7-m3-base-selector-protocol-v1"
    assert config["pre_outcome"] is True
    assert config["source_lock_required_before_fixture_generation"] is True
    assert sha(prior_result) == config["reference_result_sha256"]
    assert (config["field_degree"], config["field_modulus_hex"],
            config["subgroup_order"], config["cofactor"],
            config["isogeny_degree"], config["arity"],
            config["base_size"]) == (
                old["field_degree"], old["field_modulus_hex"],
                old["subgroup_order"], old["cofactor"],
                old["isogeny_degree"], old["arity"], old["base_size"])
    assert len(config["base_seeds"]) == len(set(config["base_seeds"])) == 4
    assert set(config["base_seeds"]).isdisjoint(old["base_seeds"])
    assert config["candidates_per_role_seed"] == 2
    assert config["source_signed_orbit_quota"] == 2
    assert config["holdouts"] == ["A", "B"]
    assert config["accepted_relation_targets_per_holdout"] == 512
    assert config["max_target_draws_per_holdout"] == 8192
    assert config["max_base_x_trials_per_candidate"] == 4096
    assert set(config["labels"]) == {"secret", "base", "target"}
    assert all(config["labels"][key] != old["labels"][key]
               for key in ("secret", "base", "target"))
    assert (config["score"]["pair_additions_per_candidate"],
            config["score"]["third_additions_per_candidate"],
            config["score"]["eligible_base_row_rank"]) == (36, 288, 8)
    assert config["producer_wall_cap_seconds"] == 240
    assert config["verifier_wall_cap_seconds"] == 360
    assert config["process_rss_cap_bytes"] == 1 << 30
    return {"status": "PASS_PRE_OUTCOME", "fixtures_generated": False,
            "config_sha256": sha(HERE / "CONFIG.json"),
            "protocol_sha256": sha(HERE / "PROTOCOL.md"),
            "reference_result_sha256": sha(prior_result),
            "reference_config_sha256": sha(OLD / "CONFIG.json")}


if __name__ == "__main__":
    print(json.dumps(check(), indent=2, sort_keys=True))
