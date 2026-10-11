#!/usr/bin/env python3
"""Audit 18 held-out raw cells and evaluate the frozen paired rank gate."""

from __future__ import annotations

import hashlib
import itertools
import json
from pathlib import Path
import statistics


HERE = Path(__file__).resolve().parent
LABELS = ("control", "cap400000", "rho")


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path: Path) -> dict:
    return json.loads(path.read_text())


def one(path: Path) -> dict:
    rows = [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
    assert len(rows) == 1, path
    return rows[0]


def bootstrap_median_interval(values: list[float]) -> list[float]:
    """Exact six-draw with-replacement descriptive percentile interval."""
    assert len(values) == 6
    medians = sorted(statistics.median(draw) for draw in itertools.product(values, repeat=6))
    return [medians[round((len(medians) - 1) * tail)] for tail in (0.025, 0.975)]


def main() -> None:
    frozen = load(HERE / "FROZEN.json")
    heldout = load(HERE / "HELDOUT_FROZEN.json")
    q_path = HERE / "inputs/heldout_q.jsonl"
    q = json.loads(q_path.read_text())
    fixture = one(HERE / "inputs/heldout_fixture_verifier_only.jsonl")
    assert heldout["public_point"] == fixture["published_q"] == q
    assert sha(q_path) == heldout["public_point_sha256"]
    assert heldout["selected_cap"] == 400000
    assert len(heldout["rank_seeds"]) == len(heldout["rho_seeds"]) == 6
    rows: dict[tuple[int, str], dict] = {}
    scalar = fixture["recovered_fixture_scalar"]
    for rank_seed, rho_seed in zip(heldout["rank_seeds"], heldout["rho_seeds"], strict=True):
        workload = load(HERE / f"heldout_workload_{rank_seed}.json")
        assert workload["workload_id"] == heldout["workload_ids"][str(rank_seed)]
        assert workload["record"]["target"] == q
        assert workload["record"]["rank_seed"] == rank_seed
        assert workload["record"]["rho_seed"] == rho_seed
        for label in LABELS:
            directory = HERE / "runs/heldout" / str(rank_seed) / label
            status = load(directory / "status.json")
            intent = load(directory / "intent.json")
            assert status["status"] == "VERIFIED"
            assert status["exit_code"] == 0 and status["timed_out"] is False
            assert status["peak_rss_bytes"] <= heldout["resource_envelope"]["observed_rss_cap_bytes"]
            assert status["workload_id"] == workload["workload_id"]
            assert status["rank_seed"] == rank_seed and status["rho_seed"] == rho_seed
            assert status["runner_sha256"] == heldout["runner_sha256"]
            assert status["public_point_sha256"] == heldout["public_point_sha256"]
            assert all(status[key] == value for key, value in intent.items())
            if label == "rho":
                assert status["source_sha256"] == heldout["rho_source_sha256"]
                assert status["binary_sha256"] == heldout["rho_binary_sha256"]
                rho_rows = [json.loads(line) for line in (directory / "rho.jsonl").read_text().splitlines()]
                assert len(rho_rows) == 2
                rho, summary = rho_rows
                replay = load(directory / "rho_replay.json")
                assert replay["status"] == "PASS" and replay["point"] == q
                assert replay["recovered_scalar"] == status["recovered_scalar"] == scalar
                assert rho["batch_seed"] == summary["batch_seed"] == rho_seed
                assert rho["published_q"] == q and rho["target_source"] == "public_point"
                assert (summary["rung"], summary["lanes"], summary["dp_bits"]) == (3, 32, 4)
                assert rho["verified"] is True and summary["all_verified"] is True
                assert abs(rho["online_ms"] - sum(rho[key] for key in
                           ("walk_ms", "collision_ms", "recovery_check_ms"))) < 1e-6
                assert status["rho_online_ms"] == rho["online_ms"]
                assert sha(directory / "rho.jsonl") == status["rho_sha256"]
                assert sha(directory / "rho_replay.json") == status["rho_replay_sha256"]
                rows[(rank_seed, label)] = {
                    "run_id": status["run_id"], "status": status["status"],
                    "rho_online_ms": rho["online_ms"], "rho_walk_steps": rho["walk_steps"],
                    "rho_online_phases_ms": {key: rho[key] for key in
                                             ("walk_ms", "collision_ms", "recovery_check_ms")},
                    "rho_charges": summary["charges"],
                    "external_wall_ms": status["external_wall_ms"],
                    "peak_rss_bytes": status["peak_rss_bytes"],
                    "rho_sha256": status["rho_sha256"],
                }
            else:
                assert status["source_sha256"] == heldout["ic_source_sha256"]
                assert status["binary_sha256"] == heldout["ic_binary_sha256"]
                assert status["run_id"].startswith(heldout["ic_candidate_ids"][label])
                base = one(directory / "base.jsonl")
                summary = one(directory / "summary.jsonl")
                target = one(directory / "targets.jsonl")
                trace = [json.loads(line) for line in (directory / "rank.jsonl").read_text().splitlines()]
                attempts = trace[1:-1]
                cap = None if label == "control" else 400000
                assert base["base_hash"] == frozen["base_hash_blake3"]
                assert sha(directory / "base.jsonl") == frozen["base_file_sha256"]
                assert base["factor_base_points"] == 23320 and base["orbit_columns"] == 220
                assert summary["rank"] == trace[-1]["rank"] == 220
                assert summary["rank_probe_cap"] == trace[0]["rank_probe_cap"] == cap
                assert summary["rank_attempts"] == len(attempts)
                assert summary["rank_relations"] == sum(row["found"] for row in attempts) == 220
                assert summary["rank_probes_total"] == sum(row["probes"] for row in attempts)
                assert summary["rank_capped_attempts"] == sum(row["capped"] for row in attempts)
                assert all(row["rank_after"] >= row["rank_before"] for row in attempts)
                assert all(not row["capped"] or (row["probes"] == cap and not row["found"])
                           for row in attempts)
                assert target["target"] == q and target["published_fixture_scalar"] is None
                assert target["group_verified"] is True
                assert target["recovered_scalar"] == status["recovered_scalar"] == scalar
                assert abs(sum(summary["cold_phase_ms"].values()) - summary["cold_in_process_ms"]) < 1e-6
                assert abs(target["target_phase_sum_ms"] - target["online_ms"]) < 1e-6
                assert load(directory / "rank_replay.json")["status"] == "PASS"
                assert load(directory / "target_replay.json")["status"] == "PASS"
                for name, field in (("rank.jsonl", "rank_trace_sha256"),
                                    ("targets.jsonl", "target_sha256"),
                                    ("rank_replay.json", "rank_replay_sha256"),
                                    ("target_replay.json", "target_replay_sha256")):
                    assert sha(directory / name) == status[field]
                rows[(rank_seed, label)] = {
                    "run_id": status["run_id"], "status": status["status"],
                    "rank_attempts": summary["rank_attempts"],
                    "rank_capped_attempts": summary["rank_capped_attempts"],
                    "rank_probes_total": summary["rank_probes_total"],
                    "rank_relations": summary["rank_relations"],
                    "rank": summary["rank"],
                    "cold_in_process_ms": summary["cold_in_process_ms"],
                    "cold_phase_ms": summary["cold_phase_ms"],
                    "target_online_ms": target["online_ms"],
                    "target_online_phases_ms": {key: target[key] for key in
                                                ("target_query_ms", "target_pdp_ms",
                                                 "target_relation_check_ms", "target_descent_ms",
                                                 "target_recovery_check_ms")},
                    "external_wall_ms": status["external_wall_ms"],
                    "peak_rss_bytes": status["peak_rss_bytes"],
                    "rank_trace_sha256": status["rank_trace_sha256"],
                }

    paired = []
    for seed in heldout["rank_seeds"]:
        control = rows[(seed, "control")]
        capped = rows[(seed, "cap400000")]
        rho = rows[(seed, "rho")]
        paired.append({
            "rank_seed": seed,
            "rho_seed": heldout["rho_seeds"][seed - heldout["rank_seeds"][0]],
            "probe_reduction_fraction": 1 - capped["rank_probes_total"] / control["rank_probes_total"],
            "cold_ratio_capped_to_control": capped["cold_in_process_ms"] / control["cold_in_process_ms"],
            "observed_rho_to_control_online_ratio": rho["rho_online_ms"] / control["target_online_ms"],
            "observed_rho_to_capped_online_ratio": rho["rho_online_ms"] / capped["target_online_ms"],
        })
    probe_reductions = [row["probe_reduction_fraction"] for row in paired]
    cold_ratios = [row["cold_ratio_capped_to_control"] for row in paired]
    gate = {
        "all_18_cells_verified": len(rows) == 18,
        "median_paired_probe_reduction_fraction": statistics.median(probe_reductions),
        "median_paired_cold_ratio": statistics.median(cold_ratios),
        "probe_threshold_fraction": 0.15,
        "requires_no_paired_cold_increase": True,
        "passed": (len(rows) == 18 and statistics.median(probe_reductions) >= 0.15
                   and statistics.median(cold_ratios) <= 1),
    }
    result = {
        "schema": "n53-fixed-base-rank-restart-heldout-analysis-v1",
        "status": "PASS",
        "timing_class": "exploratory_shared_host",
        "controlled_online_speedup": None,
        "public_target": q,
        "independently_replayed_scalar": scalar,
        "base_hash_blake3": frozen["base_hash_blake3"],
        "candidate_ids": heldout["ic_candidate_ids"],
        "rank_seeds": heldout["rank_seeds"],
        "rho_seeds": heldout["rho_seeds"],
        "cells": [{"seed": seed, "label": label, **rows[(seed, label)]}
                  for seed in heldout["rank_seeds"] for label in LABELS],
        "paired": paired,
        "probe_reduction_range": [min(probe_reductions), max(probe_reductions)],
        "cold_ratio_range": [min(cold_ratios), max(cold_ratios)],
        "descriptive_bootstrap_95_median_probe_reduction": bootstrap_median_interval(probe_reductions),
        "descriptive_bootstrap_95_median_cold_ratio": bootstrap_median_interval(cold_ratios),
        "gate": gate,
    }
    output = HERE / "HELDOUT_ANALYSIS.json"
    encoded = json.dumps(result, sort_keys=True, indent=2) + "\n"
    if output.exists():
        assert output.read_text() == encoded, "held-out analysis changed"
    else:
        output.write_text(encoded)
    print(json.dumps({"status": result["status"], "gate": gate}, sort_keys=True))


if __name__ == "__main__":
    main()
