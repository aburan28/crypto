#!/usr/bin/env python3
"""Reduce the two independently verified, five-block panel receipts."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import statistics

from verify_panel import ARMS, HERE, sha


def analyze(runs: Path) -> dict:
    cells = {}
    for n in (41, 53):
        path = runs / f"n{n}_L1024" / "receipt.json"
        receipt = json.loads(path.read_text())
        assert receipt["status"] == "PASS" and receipt["n"] == n
        assert receipt["targets_verified_per_arm"] == 5 * 1024
        assert receipt["rank_runs_verified"] == 5 * 3
        measurements = receipt["measurements"]
        assert set(measurements) == {str(i) for i in range(5)}
        checks = receipt["checks"]
        assert len(checks) == 5 * len(ARMS)
        table = {}
        for arm in ARMS:
            rows = [measurements[str(block)][arm] for block in range(5)]
            table[arm] = {
                "cpu_median_seconds": statistics.median(row["cpu_seconds"] for row in rows),
                "wall_median_seconds": statistics.median(row["wall_seconds"] for row in rows),
                "rss_median_mib": statistics.median(row["rss_kib"] for row in rows) / 1024,
                "logs_independently_verified": 5 * 1024,
                "S": None,
            }
            if arm != "rho_normal":
                arm_checks = [item for item in checks if item["arm"] == arm]
                assert len(arm_checks) == 5
                table[arm]["prefilter_bytes"] = statistics.median(
                    item["root_prefilter_bytes"] for item in arm_checks)
                table[arm]["phase_median_ms"] = {
                    phase: statistics.median(item["timing_ms"][phase] for item in arm_checks)
                    for phase in ("index_build", "rank_stage", "linear_algebra", "targets_total")}
                table[arm]["root_keys_considered_median"] = statistics.median(
                    item["rank_s3_counts"]["root_keys_considered"]
                    + item["target_s3_counts"]["root_keys_considered"]
                    for item in arm_checks)
                table[arm]["filter_definite_misses_median"] = statistics.median(
                    item["rank_s3_counts"]["prefilter_definite_misses"]
                    + item["target_s3_counts"]["prefilter_definite_misses"]
                    for item in arm_checks)
                table[arm]["filter_false_positives_median"] = statistics.median(
                    item["rank_s3_counts"]["prefilter_false_positives"]
                    + item["target_s3_counts"]["prefilter_false_positives"]
                    for item in arm_checks)
        cells[str(n)] = {
            "receipt_sha256": sha(path), "k": receipt["k"],
            "table": table, "paired": receipt["paired"],
            "aa_timing_valid": receipt["aa_timing_valid"],
            "uncontended": receipt["uncontended"],
            "timing_eligible": receipt["timing_eligible"],
            "engineering_win": receipt["engineering_win"],
            "timing_crossover_cell": receipt["timing_crossover_cell"],
        }
    return {
        "schema": "compact-s3-prefilter-analysis-v1",
        "status": "PASS", "frozen_sha256": sha(HERE / "FROZEN.json"),
        "cells": cells,
        "engineering_win_both_sizes": all(cell["engineering_win"] for cell in cells.values()),
        "timing_crossover_both_sizes": all(
            cell["timing_crossover_cell"] for cell in cells.values()),
        "S_off": None, "S_filter": None, "S_rho_normal": None,
        "method_level_speedup": None,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--runs", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    assert not args.out.exists(), "never overwrite analysis"
    result = analyze(args.runs.resolve())
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({name: value for name, value in result.items()
                      if name != "cells"}, sort_keys=True))


if __name__ == "__main__":
    main()
