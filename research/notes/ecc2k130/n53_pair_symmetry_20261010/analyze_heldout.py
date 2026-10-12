#!/usr/bin/env python3
"""Independently audit raw paired IC/rho cells and the frozen operation gate."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path
import statistics


HERE = Path(__file__).resolve().parent
LABELS = ("control", "symmetry", "rho")
ROW_FIELDS = ("scalar", "target", "pivotless_column", "point_indices", "x_codes",
              "pinned_intermediates", "row", "rank_before", "rank_after")


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path: Path) -> dict:
    return json.loads(path.read_text())


def one(path: Path) -> dict:
    rows = [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
    assert len(rows) == 1, path
    return rows[0]


def interval(values: list[float]) -> list[float]:
    """Exact six-draw bootstrap, descriptive because host isolation is absent."""
    assert len(values) == 6
    medians = sorted(statistics.median(draw) for draw in itertools.product(values, repeat=6))
    return [medians[round((len(medians) - 1) * q)] for q in (0.025, 0.975)]


def semantic_rows(trace: list[dict]) -> list[dict]:
    return [{key: row[key] for key in ROW_FIELDS} for row in trace[1:-1]
            if row["found"] and row["gained"]]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--runs-root", type=Path, default=HERE / "runs/heldout")
    args = parser.parse_args()
    prep = load(HERE / "FROZEN_PREP.json")
    frozen = load(HERE / "HELDOUT_FROZEN.json")
    q_path = HERE / "inputs/heldout_q.jsonl"
    q = json.loads(q_path.read_text())
    fixture = one(HERE / "inputs/heldout_fixture_verifier_only.jsonl")
    assert frozen["public_point"] == fixture["published_q"] == q
    assert sha(q_path) == frozen["public_point_sha256"]
    assert sha(HERE / "inputs/heldout_fixture_verifier_only.jsonl") == frozen["fixture_verifier_only_sha256"]
    assert frozen["ic_candidate_ids"] == prep["candidate_ids"]
    assert len(frozen["rank_seeds"]) == len(frozen["rho_seeds"]) == 6
    scalar = int(fixture["recovered_fixture_scalar"])
    cells = []
    rows = {}
    all_verified = True
    for seed, rho_seed in zip(frozen["rank_seeds"], frozen["rho_seeds"], strict=True):
        workload = load(HERE / f"heldout_workload_{seed}.json")
        assert workload["workload_id"] == frozen["workload_ids"][str(seed)]
        assert (workload["record"]["target"], workload["record"]["rank_seed"],
                workload["record"]["rho_seed"]) == (q, seed, rho_seed)
        for label in LABELS:
            directory = args.runs_root / str(seed) / label
            status = load(directory / "status.json")
            intent = load(directory / "intent.json")
            assert all(status[key] == value for key, value in intent.items())
            assert status["rank_seed"] == seed and status["rho_seed"] == rho_seed
            assert status["workload_id"] == workload["workload_id"]
            assert status["public_point_sha256"] == frozen["public_point_sha256"]
            assert status["runner_sha256"] == frozen["runner_sha256"]
            assert status["resource_envelope"] == frozen["resource_envelope"]
            cell = {"rank_seed": seed, "rho_seed": rho_seed, "label": label,
                    "run_id": status["run_id"], "status": status["status"],
                    "exit_code": status["exit_code"], "timed_out": status["timed_out"],
                    "external_wall_ms": status["external_wall_ms"],
                    "peak_rss_bytes": status["peak_rss_bytes"],
                    "status_sha256": sha(directory / "status.json")}
            if status["status"] != "VERIFIED":
                all_verified = False
                cell["error"] = status.get("error")
                cells.append(cell)
                continue
            assert status["exit_code"] == 0 and status["timed_out"] is False
            assert status["peak_rss_bytes"] <= frozen["resource_envelope"]["observed_rss_cap_bytes"]
            assert status["recovered_scalar"] == scalar
            if label == "rho":
                assert status["source_sha256"] == frozen["rho_source_sha256"]
                assert status["binary_sha256"] == frozen["rho_binary_sha256"]
                rho, summary = [json.loads(line) for line in (directory / "rho.jsonl").read_text().splitlines()]
                replay = load(directory / "rho_replay.json")
                assert replay["status"] == "PASS" and replay["point"] == q
                assert replay["recovered_scalar"] == scalar
                assert rho["batch_seed"] == summary["batch_seed"] == rho_seed
                assert rho["published_q"] == q and rho["target_source"] == "public_point"
                assert (summary["rung"], summary["lanes"], summary["dp_bits"]) == (3, 32, 4)
                assert rho["verified"] is True and summary["all_verified"] is True
                assert abs(rho["online_ms"] - sum(rho[key] for key in
                           ("walk_ms", "collision_ms", "recovery_check_ms"))) < 1e-6
                assert sha(directory / "rho.jsonl") == status["rho_sha256"]
                assert sha(directory / "rho_replay.json") == status["rho_replay_sha256"]
                cell.update({"rho_online_ms": rho["online_ms"],
                             "rho_walk_steps": rho["walk_steps"],
                             "rho_phases_ms": {key: rho[key] for key in
                                               ("walk_ms", "collision_ms", "recovery_check_ms")},
                             "rho_charges": summary["charges"]})
            else:
                assert status["source_sha256"] == frozen["ic_source_sha256"]
                assert status["binary_sha256"] == frozen["ic_binary_sha256"]
                assert status["run_id"].startswith(frozen["ic_candidate_ids"][label])
                base = one(directory / "base.jsonl")
                summary = one(directory / "summary.jsonl")
                target = one(directory / "targets.jsonl")
                trace = [json.loads(line) for line in (directory / "rank.jsonl").read_text().splitlines()]
                attempts = trace[1:-1]
                assert base["base_hash"] == prep["base_hash_blake3"]
                assert sha(directory / "base.jsonl") == prep["base_file_sha256"]
                assert base["factor_base_points"] == 23_320 and base["orbit_columns"] == 220
                assert summary["rank"] == trace[-1]["rank"] == 220
                assert summary["rank_probe_cap"] is None and trace[0]["rank_probe_cap"] is None
                assert summary["pair_symmetry"] is (label == "symmetry")
                assert trace[0]["pair_symmetry"] is (label == "symmetry")
                assert summary["regular_states"] == 2_565_200
                assert summary["root_table_entries"] == 2_564_528
                assert summary["index_pair_candidates"] == (1_282_710 if label == "symmetry" else 2_565_200)
                assert summary["rank_attempts"] == len(attempts)
                assert summary["rank_probes_total"] == sum(row["probes"] for row in attempts)
                assert summary["rank_root_calls_total"] == sum(row["root_calls"] for row in attempts)
                assert summary["rank_skipped_symmetric_states"] == sum(
                    row["skipped_symmetric_states"] for row in attempts)
                assert len(semantic_rows(trace)) == 220
                assert target["target"] == q and target["published_fixture_scalar"] is None
                assert target["group_verified"] is True and target["recovered_scalar"] == scalar
                assert abs(sum(summary["cold_phase_ms"].values()) - summary["cold_in_process_ms"]) < 1e-6
                assert abs(target["target_phase_sum_ms"] - target["online_ms"]) < 1e-6
                assert load(directory / "rank_replay.json")["status"] == "PASS"
                assert load(directory / "target_replay.json")["status"] == "PASS"
                for name, field in (("rank.jsonl", "rank_trace_sha256"),
                                    ("targets.jsonl", "target_sha256"),
                                    ("rank_replay.json", "rank_replay_sha256"),
                                    ("target_replay.json", "target_replay_sha256")):
                    assert sha(directory / name) == status[field]
                cell.update({"index_pair_candidates": summary["index_pair_candidates"],
                             "regular_states": summary["regular_states"],
                             "root_table_entries": summary["root_table_entries"],
                             "rank_attempts": summary["rank_attempts"],
                             "rank_probes_total": summary["rank_probes_total"],
                             "rank_root_calls_total": summary["rank_root_calls_total"],
                             "rank_skipped_symmetric_states": summary["rank_skipped_symmetric_states"],
                             "rank": summary["rank"], "cold_in_process_ms": summary["cold_in_process_ms"],
                             "cold_phase_ms": summary["cold_phase_ms"],
                             "target_online_ms": target["online_ms"],
                             "target_online_phases_ms": {key: target[key] for key in
                                                         ("target_query_ms", "target_pdp_ms",
                                                          "target_relation_check_ms", "target_descent_ms",
                                                          "target_recovery_check_ms")},
                             "rank_rows_sha256": hashlib.sha256(json.dumps(
                                 semantic_rows(trace), sort_keys=True,
                                 separators=(",", ":")).encode()).hexdigest()})
            rows[(seed, label)] = cell
            cells.append(cell)
    paired = []
    if all_verified and len(rows) == 18:
        for seed, rho_seed in zip(frozen["rank_seeds"], frozen["rho_seeds"], strict=True):
            control, symmetry, rho = (rows[(seed, label)] for label in LABELS)
            paired.append({
                "rank_seed": seed, "rho_seed": rho_seed,
                "same_220_rank_rows": control["rank_rows_sha256"] == symmetry["rank_rows_sha256"],
                "rank_root_calls_nonincrease": symmetry["rank_root_calls_total"] <= control["rank_root_calls_total"],
                "rank_root_call_reduction": control["rank_root_calls_total"] - symmetry["rank_root_calls_total"],
                "rank_probe_reduction": control["rank_probes_total"] - symmetry["rank_probes_total"],
                "index_ratio_symmetry_to_control": (
                    symmetry["cold_phase_ms"]["precompute_index"] /
                    control["cold_phase_ms"]["precompute_index"]),
                "rank_pdp_ratio_symmetry_to_control": (
                    symmetry["cold_phase_ms"]["rank_pdp"] /
                    control["cold_phase_ms"]["rank_pdp"]),
                "cold_ratio_symmetry_to_control": symmetry["cold_in_process_ms"] / control["cold_in_process_ms"],
                "target_online_ratio_symmetry_to_control": symmetry["target_online_ms"] / control["target_online_ms"],
                "rho_to_control_online_ratio": rho["rho_online_ms"] / control["target_online_ms"],
                "rho_to_symmetry_online_ratio": rho["rho_online_ms"] / symmetry["target_online_ms"],
            })
    gate = {"all_18_cells_verified": all_verified and len(rows) == 18,
            "same_rank_rows_each_pair": bool(paired) and all(row["same_220_rank_rows"] for row in paired),
            "rank_root_calls_nonincrease_each_pair": bool(paired) and all(
                row["rank_root_calls_nonincrease"] for row in paired),
            "index_pair_calls_control": 2_565_200,
            "index_pair_calls_symmetry": 1_282_710}
    gate["passed"] = all(gate[key] for key in ("all_18_cells_verified",
                                                  "same_rank_rows_each_pair",
                                                  "rank_root_calls_nonincrease_each_pair"))
    result = {"schema": "n53-pair-symmetry-heldout-analysis-v1",
              "status": "COMPLETE" if gate["all_18_cells_verified"] else "INCOMPLETE",
              "timing_class": "exploratory_shared_host", "controlled_online_speedup": None,
              "public_target": q, "independently_replayed_scalar": scalar,
              "base_hash_blake3": prep["base_hash_blake3"],
              "candidate_ids": frozen["ic_candidate_ids"],
              "rank_seeds": frozen["rank_seeds"], "rho_seeds": frozen["rho_seeds"],
              "cells": cells, "paired": paired, "gate": gate}
    if paired:
        for metric in ("index_ratio_symmetry_to_control", "rank_pdp_ratio_symmetry_to_control",
                       "cold_ratio_symmetry_to_control", "target_online_ratio_symmetry_to_control",
                       "rho_to_control_online_ratio", "rho_to_symmetry_online_ratio"):
            values = [row[metric] for row in paired]
            result[f"median_{metric}"] = statistics.median(values)
            result[f"descriptive_bootstrap_95_median_{metric}"] = interval(values)
    output = HERE / "HELDOUT_ANALYSIS.json"
    encoded = json.dumps(result, sort_keys=True, indent=2) + "\n"
    if args.check:
        assert output.read_text() == encoded, "held-out analysis changed"
    elif output.exists():
        assert output.read_text() == encoded, "refusing to replace held-out analysis"
    else:
        output.write_text(encoded)
    print(json.dumps({"status": result["status"], "gate": gate,
                      "paired_seeds": len(paired)}, sort_keys=True))


if __name__ == "__main__":
    main()
