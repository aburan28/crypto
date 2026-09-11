#!/usr/bin/env python3
"""Verify the committed independent replay of the Stage-21 pair oracle."""

from __future__ import annotations

import json
import math
from pathlib import Path

import run_koblitz_blind_pdp_phase_b as phase_b


REPO = Path(__file__).resolve().parents[1]
STAGE = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
ROOT = STAGE / "stage-28-independent-relation-yield-replay-20260911"


class VerificationError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def load(path: Path) -> dict:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"JSON object required: {path}")
    return value


def verify() -> dict:
    seal = load(ROOT / "replay-seal.json")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    require(
        seal.get("schema") == "koblitz_stage28_independent_relation_yield_replay_seal.v1"
        and seal.get("status") == "replay_frozen"
        and claimed == phase_b.canonical_sha256(payload),
        "Stage-28 replay seal is invalid",
    )
    inventory = phase_b.all_regular_inventory(ROOT, {"replay-seal.json"})
    require(inventory == seal.get("inventory") and phase_b.canonical_sha256(inventory) == seal.get("inventory_sha256"), "Stage-28 replay inventory changed")
    stage21_seal = load(STAGE / "stage-21-relation-yield-result-seal-20260910.json")
    input_identity = stage21_seal.get("production_artifacts", {}).get("yield_result", {})
    input_path = ROOT / "input-yield-result.json"
    require(
        input_path.stat().st_size == input_identity.get("bytes")
        and phase_b.sha256_file(input_path, "Stage-28 input") == input_identity.get("sha256")
        and seal.get("source_stage21_yield_result_sha256") == input_identity.get("sha256")
        and seal.get("source_stage21_run_inventory_sha256") == stage21_seal.get("cross_bindings", {}).get("run_inventory_sha256"),
        "Stage-28 input differs from the sealed Stage-21 producer result",
    )
    source = load(input_path)
    replay = load(ROOT / "replay.json")
    require(source.get("schema") == "koblitz_relation_yield_bridge.v1" and source.get("status") == "complete", "Stage-28 source schema changed")
    require(replay.get("schema") == "koblitz_stage28_independent_relation_yield_replay.v1" and replay.get("status") == "PASS", "Stage-28 replay status changed")
    require(
        replay.get("factor_base_materialization_replayed") is True
        and replay.get("canonical_pair_oracle_replayed") is True
        and replay.get("witness_curve_readdition_replayed") is True,
        "Stage-28 did not complete all three payload replays",
    )
    require(replay.get("target_subgroup_enumerated") is False and replay.get("discrete_log_labels_used") is False, "Stage-28 widened the information boundary")
    factor = replay.get("factor_base", {})
    source_factor = source.get("factor_base", {})
    require(
        factor == {
            "abscissae": 4096,
            "dimension": 12,
            "factor_base_blake3": source_factor.get("factor_base_blake3"),
            "rational_points": 4281,
        },
        "Stage-28 factor-base replay changed",
    )
    pair = replay.get("pair_oracle", {})
    source_pair = source.get("exact_pair_table", {})
    for replay_field, source_field in (
        ("canonical_pairs", "canonical_pairs"),
        ("enumerated_pairs", "enumerated_pairs"),
        ("unique_target_entries", "unique_target_entries"),
        ("duplicate_pair_sums", "duplicate_pair_sums"),
        ("canonical_pair_transcript_blake3", "canonical_pair_transcript_blake3"),
    ):
        require(pair.get(replay_field) == source_pair.get(source_field), f"Stage-28 pair replay changed at {replay_field}")
    require(replay.get("rows") == 384 and replay.get("checked_witnesses") == 227, "Stage-28 row or witness count changed")
    require(
        replay.get("arms")
        == {
            "natural": {"hits": 163, "targets": 256, "verified_witnesses": 163},
            "planted_sat": {"hits": 64, "targets": 64, "verified_witnesses": 64},
            "proven_unsat": {"hits": 0, "targets": 64, "verified_witnesses": 0},
        },
        "Stage-28 arm replay changed",
    )
    metrics = load(ROOT / "replay.metrics.json")
    resources = metrics.get("metrics", {})
    require(
        metrics.get("returncode") == 0
        and metrics.get("timed_out") is False
        and metrics.get("orphan_group_terminated") is False
        and isinstance(resources.get("total_core_seconds"), (int, float))
        and math.isfinite(resources["total_core_seconds"])
        and resources["total_core_seconds"] > 0
        and isinstance(resources.get("wall_seconds"), (int, float))
        and resources["wall_seconds"] > 0
        and isinstance(resources.get("peak_rss_bytes"), int)
        and resources["peak_rss_bytes"] > 0,
        "Stage-28 process receipt is incomplete",
    )
    require(replay.get("independent_external_reproduction_satisfied") is False and replay.get("full_cost_gate_passed") is False and replay.get("koblitz_index_calculus_sota") is False, "Stage-28 replay widened its conclusion")
    return {
        "schema": "koblitz_stage28_relation_yield_replay_verification.v1",
        "status": "independent_internal_factor_base_pair_oracle_and_witness_replay_verified",
        "input_sha256": input_identity["sha256"],
        "replay_sha256": phase_b.sha256_file(ROOT / "replay.json", "Stage-28 replay"),
        "inventory_sha256": seal["inventory_sha256"],
        "canonical_pairs": pair["canonical_pairs"],
        "checked_witnesses": replay["checked_witnesses"],
        "total_core_seconds": resources["total_core_seconds"],
        "wall_seconds": resources["wall_seconds"],
        "peak_rss_bytes": resources["peak_rss_bytes"],
        "scientific_measurement_admission": "finite_public_synthetic_internal_replay_complete",
        "independent_external_reproduction_satisfied": False,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
    }


if __name__ == "__main__":
    try:
        print(json.dumps(verify(), indent=2, sort_keys=True))
    except (OSError, ValueError, KeyError, VerificationError, phase_b.PhaseBError) as error:
        raise SystemExit(f"stage28-replay: {error}")
