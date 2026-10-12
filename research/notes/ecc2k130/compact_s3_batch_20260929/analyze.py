#!/usr/bin/env python3
"""Derive one-unit full-process CPU table and all pre-registered ratios."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import statistics


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def analyze(run_dir: Path) -> dict:
    receipt_path, host_path = run_dir / "receipt.json", run_dir / "host.json"
    receipt, host = json.loads(receipt_path.read_text()), json.loads(host_path.read_text())
    assert receipt["status"] == "PASS"
    assert run_dir.name == f"n{receipt['n']}_L1024" and receipt["L"] == 1024
    assert receipt["targets_verified_per_arm"] == 9 * 1024
    assert receipt["rank_runs_verified"] == 9 * 8
    assert len(receipt["checks"]) == 9 * 9
    assert sorted(int(block) for block in receipt["measurements"]) == list(range(9))
    grid = receipt["k_grid"]
    assert len(grid) == 2
    arms = [f"{policy}_k{k}" for k in grid
            for policy in ("scalar_a", "batch16", "batch64", "scalar_b")]
    arms.append("rho_normal")
    table = []
    for arm in arms:
        measures = [receipt["measurements"][str(block)][arm] for block in range(9)]
        checks = [check for check in receipt["checks"] if check["arm"] == arm]
        assert len(checks) == 9
        row = {
            "arm": arm,
            "cpu_median_seconds": statistics.median(m["cpu_seconds"] for m in measures),
            "wall_median_seconds": statistics.median(m["wall_seconds"] for m in measures),
            "rss_median_mib": statistics.median(m["rss_kib"] for m in measures) / 1024,
            "targets_verified": sum(c["targets_verified"] for c in checks),
        }
        if arm == "rho_normal":
            row["walk_steps_median"] = statistics.median(c["walk_steps"] for c in checks)
            row["batch_inversion_calls_median"] = statistics.median(
                c["charges"]["batch_inversion_calls"] for c in checks)
        else:
            phases = ("index_build", "rank_stage", "linear_algebra", "targets_total")
            row.update({
                "regular_states": checks[0]["regular_states"],
                "root_table_entries": checks[0]["root_table_entries"],
                "root_table_slots": checks[0]["root_table_slots"],
                "rank_attempts_median": statistics.median(c["rank_attempts"] for c in checks),
                "rank_probes_mean": statistics.mean(c["rank_probes_mean"] for c in checks),
                "target_probes_mean": statistics.mean(c["target_probes_mean"] for c in checks),
                "phase_median_seconds": {
                    phase: statistics.median(c["timing_ms"][phase] for c in checks) / 1000
                    for phase in phases},
                "s3_counts_median": {
                    phase: {key: statistics.median(c[f"{phase}_s3_counts"][key] for c in checks)
                            for key in checks[0][f"{phase}_s3_counts"]}
                    for phase in ("index", "rank", "target")},
            })
            k = int(arm.rsplit("_k", 1)[1])
            if arm.startswith("batch"):
                row["paired_to_scalar_a"] = receipt["paired"][f"{arm}_to_scalar_a"]
                row["paired_to_rho"] = receipt["paired"][f"{arm}_to_rho"]
            elif arm.startswith("scalar_b"):
                row["paired_to_scalar_a"] = receipt["paired"][f"scalar_b_k{k}_to_scalar_a"]
            else:
                row["paired_to_rho"] = receipt["paired"][f"scalar_a_k{k}_to_rho"]
        table.append(row)
    candidates = [row for row in table if row["arm"].startswith("batch")]
    best = min(candidates, key=lambda row: (row["cpu_median_seconds"], row["arm"]))
    return {
        "schema": "compact-s3-batch-derived-v1",
        "n": receipt["n"], "L": receipt["L"],
        "host_cpu_model": host["cpu_model"], "host_git_checkout": host["git_head"],
        "receipt_sha256": sha(receipt_path), "host_sha256": sha(host_path),
        "table": table,
        "best_exploratory_batch_arm": best["arm"],
        "aa_timing_valid": receipt["aa_timing_valid"],
        "all_batch_cpu_intervals_above_rho": all(
            row["paired_to_rho"]["cpu_ratio_95pct_log_t_interval"][0] > 1
            for row in candidates),
        "S_scalar": receipt["S_scalar"], "S_batch16": receipt["S_batch16"],
        "S_batch64": receipt["S_batch64"], "S_rho_normal": receipt["S_rho_normal"],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", type=Path, action="append", required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    assert not args.out.exists(), "never overwrite derived evidence"
    results = [analyze(path) for path in args.run_dir]
    assert sorted(row["n"] for row in results) == [41, 53]
    args.out.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n")
    print(json.dumps({row["n"]: {
        "best_batch_arm": row["best_exploratory_batch_arm"],
        "all_batch_cpu_intervals_above_rho": row["all_batch_cpu_intervals_above_rho"],
        "aa_timing_valid": row["aa_timing_valid"],
    } for row in results}, sort_keys=True))


if __name__ == "__main__":
    main()
