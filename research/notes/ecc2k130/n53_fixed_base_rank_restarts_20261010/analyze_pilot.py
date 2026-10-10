#!/usr/bin/env python3
"""Audit every frozen pilot cell and select the preregistered capped arm."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import statistics


HERE = Path(__file__).resolve().parent
SEEDS = (530053, 530054, 530055)
LABELS = ("control", "cap200000", "cap400000")


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path: Path) -> dict:
    return json.loads(path.read_text())


def one_jsonl(path: Path) -> dict:
    rows = [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
    assert len(rows) == 1, path
    return rows[0]


def main() -> None:
    frozen = load(HERE / "FROZEN.json")
    point = json.loads((HERE / "inputs/development_q.jsonl").read_text())
    rows = {}
    scalar = None
    for seed in SEEDS:
        for label in LABELS:
            directory = HERE / "runs/pilot" / str(seed) / label
            status = load(directory / "status.json")
            summary = one_jsonl(directory / "summary.jsonl")
            target = one_jsonl(directory / "targets.jsonl")
            base = one_jsonl(directory / "base.jsonl")
            rank_replay = load(directory / "rank_replay.json")
            target_replay = load(directory / "target_replay.json")
            cap = None if label == "control" else int(label.removeprefix("cap"))
            assert status["status"] == "VERIFIED"
            assert status["exit_code"] == 0 and not status["timed_out"]
            assert status["peak_rss_bytes"] <= frozen["resource_envelope"]["observed_rss_cap_bytes"]
            assert status["source_sha256"] == frozen["source_sha256"]
            assert status["binary_sha256"] == frozen["binary_sha256"]
            assert status["candidate_id"] == frozen["candidate_ids"][label]
            assert status["workload_id"] == frozen["workload_ids"][str(seed)]
            assert status["rank_probe_cap"] == summary["rank_probe_cap"] == cap
            assert summary["rank"] == 220 and summary["targets_solved"] == 1
            assert summary["targets_failed"] == 0
            assert base["base_hash"] == frozen["base_hash_blake3"]
            assert sha(directory / "base.jsonl") == frozen["base_file_sha256"]
            assert target["target"] == point and target["group_verified"] is True
            assert target["published_fixture_scalar"] is None
            assert rank_replay["status"] == target_replay["status"] == "PASS"
            assert rank_replay["rank"] == 220
            assert target_replay["target"] == point
            assert abs(sum(summary["cold_phase_ms"].values()) - summary["cold_in_process_ms"]) < 1e-6
            assert abs(target["target_phase_sum_ms"] - target["online_ms"]) < 1e-6
            for name, field in (
                ("rank.jsonl", "rank_trace_sha256"),
                ("summary.jsonl", "summary_sha256"),
                ("targets.jsonl", "target_sha256"),
                ("rank_replay.json", "rank_replay_sha256"),
                ("target_replay.json", "target_replay_sha256"),
            ):
                assert sha(directory / name) == status[field]
            assert target["recovered_scalar"] == target_replay["recovered_scalar"]
            if scalar is None:
                scalar = target["recovered_scalar"]
            assert scalar == target["recovered_scalar"]
            rows[(seed, label)] = {
                "run_id": status["run_id"],
                "status": status["status"],
                "cold_in_process_ms": summary["cold_in_process_ms"],
                "external_wall_ms": status["external_wall_ms"],
                "peak_rss_bytes": status["peak_rss_bytes"],
                "rank_probes_total": summary["rank_probes_total"],
                "rank_attempts": summary["rank_attempts"],
                "rank_capped_attempts": summary["rank_capped_attempts"],
                "rank_relations": summary["rank_relations"],
                "rank_failures": summary["rank_failures"],
                "rank": summary["rank"],
                "target_online_ms": target["online_ms"],
                "cold_phase_ms": summary["cold_phase_ms"],
                "rank_trace_sha256": status["rank_trace_sha256"],
            }

    by_candidate = {}
    for label in LABELS:
        run_rows = [rows[(seed, label)] for seed in SEEDS]
        paired = [
            {
                "seed": seed,
                "probe_reduction_fraction": 1 - rows[(seed, label)]["rank_probes_total"]
                / rows[(seed, "control")]["rank_probes_total"],
                "cold_ratio_to_control": rows[(seed, label)]["cold_in_process_ms"]
                / rows[(seed, "control")]["cold_in_process_ms"],
            }
            for seed in SEEDS
        ]
        by_candidate[label] = {
            "candidate_id": frozen["candidate_ids"][label],
            "all_verified": all(row["status"] == "VERIFIED" for row in run_rows),
            "median_cold_ms": statistics.median(row["cold_in_process_ms"] for row in run_rows),
            "median_total_rank_probes": statistics.median(row["rank_probes_total"] for row in run_rows),
            "cold_range_ms": [min(row["cold_in_process_ms"] for row in run_rows),
                              max(row["cold_in_process_ms"] for row in run_rows)],
            "paired": paired,
            "median_paired_probe_reduction_fraction": statistics.median(
                item["probe_reduction_fraction"] for item in paired
            ),
            "median_paired_cold_ratio": statistics.median(
                item["cold_ratio_to_control"] for item in paired
            ),
        }
    caps = ("cap200000", "cap400000")
    eligible = [label for label in caps if by_candidate[label]["all_verified"]]
    selected = min(eligible, key=lambda label: (by_candidate[label]["median_cold_ms"], label))
    result = {
        "schema": "n53-fixed-base-rank-restart-pilot-analysis-v1",
        "status": "PASS",
        "timing_class": "exploratory_shared_host",
        "source_sha256": frozen["source_sha256"],
        "binary_sha256": frozen["binary_sha256"],
        "base_hash_blake3": frozen["base_hash_blake3"],
        "public_target": point,
        "independently_replayed_scalar": scalar,
        "rank_seeds": list(SEEDS),
        "cells": [
            {"seed": seed, "label": label, **rows[(seed, label)]}
            for seed in SEEDS for label in LABELS
        ],
        "by_candidate": by_candidate,
        "selected_cap_for_heldout": selected,
        "selection_rule": "lowest median complete cold ms among capped variants verified on all three pilot seeds",
    }
    output = HERE / "PILOT_ANALYSIS.json"
    encoded = json.dumps(result, sort_keys=True, indent=2) + "\n"
    if output.exists():
        assert output.read_text() == encoded, "pilot analysis changed"
    else:
        output.write_text(encoded)
    print(json.dumps({
        "status": result["status"],
        "selected": selected,
        "median_cold_ms": {label: by_candidate[label]["median_cold_ms"] for label in LABELS},
        "median_paired_probe_reduction": {
            label: by_candidate[label]["median_paired_probe_reduction_fraction"] for label in LABELS
        },
    }, sort_keys=True))


if __name__ == "__main__":
    main()
