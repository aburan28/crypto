#!/usr/bin/env python3
"""Derive the full-cost K/policy table from independently checked receipts."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import statistics


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def median(values) -> float:
    return statistics.median(values)


def analyze(run_dir: Path) -> dict:
    receipt_path = run_dir / "receipt.json"
    host_path = run_dir / "host.json"
    receipt = json.loads(receipt_path.read_text())
    host = json.loads(host_path.read_text())
    assert receipt["status"] == "PASS"
    assert run_dir.name == f"n{receipt['n']}_L1024"
    assert receipt["L"] == 1024
    assert receipt["targets_verified_per_arm"] == 9 * 1024
    assert receipt["rank_runs_verified"] == 9 * 8
    assert len(receipt["checks"]) == 9 * 9
    assert sorted(int(block) for block in receipt["measurements"]) == list(range(9))
    grid = receipt["k_grid"]
    assert len(grid) == 4
    arms = [f"{policy}_k{k}" for k in grid for policy in ("base", "swap")]
    arms.append("rho_normal")
    table = []
    for arm in arms:
        measures = [receipt["measurements"][str(block)][arm] for block in range(9)]
        checks = [check for check in receipt["checks"] if check["arm"] == arm]
        assert len(checks) == 9
        row = {
            "arm": arm,
            "cpu_median_seconds": median(m["cpu_seconds"] for m in measures),
            "wall_median_seconds": median(m["wall_seconds"] for m in measures),
            "rss_median_mib": median(m["rss_kib"] for m in measures) / 1024,
            "targets_verified": sum(c["targets_verified"] for c in checks),
        }
        if arm == "rho_normal":
            row["walk_steps_median"] = median(c["walk_steps"] for c in checks)
            row["batch_inversion_calls_median"] = median(
                c["charges"]["batch_inversion_calls"] for c in checks)
        else:
            phases = ("index_build", "rank_stage", "linear_algebra", "targets_total")
            row.update({
                "regular_states": checks[0]["regular_states"],
                "root_table_entries": checks[0]["root_table_entries"],
                "root_table_slots": checks[0]["root_table_slots"],
                "rank_attempts_median": median(c["rank_attempts"] for c in checks),
                "rank_probes_mean": statistics.mean(c["rank_probes_mean"] for c in checks),
                "target_probes_mean": statistics.mean(c["target_probes_mean"] for c in checks),
                "phase_median_seconds": {
                    phase: median(c["timing_ms"][phase] for c in checks) / 1000
                    for phase in phases
                },
            })
            k = int(arm.split("_k")[1])
            ratio = (f"{arm}_to_rho" if arm.startswith("base_")
                     else f"swap_k{k}_to_rho")
            row["paired_to_rho"] = receipt["paired"][ratio]
            if arm.startswith("swap_"):
                row["paired_to_base"] = receipt["paired"][f"swap_k{k}_to_base"]
        table.append(row)
    swap = [r for r in table if r["arm"].startswith("swap_")]
    best = min(swap, key=lambda row: (row["cpu_median_seconds"], row["arm"]))
    return {
        "schema": "compact-swap-quotient-derived-v1",
        "n": receipt["n"],
        "L": receipt["L"],
        "host_cpu_model": host["cpu_model"],
        "host_git_checkout": host["git_head"],
        "receipt_sha256": sha(receipt_path),
        "host_sha256": sha(host_path),
        "table": table,
        "exploratory_grid_minimum_swap_arm": best["arm"],
        "all_swap_cpu_intervals_above_rho": all(
            row["paired_to_rho"]["cpu_ratio_95pct_log_t_interval"][0] > 1
            for row in swap),
        "S_base": receipt["S_base"],
        "S_swap": receipt["S_swap"],
        "S_rho_normal": receipt["S_rho_normal"],
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
        "best_swap_arm": row["exploratory_grid_minimum_swap_arm"],
        "all_swap_cpu_intervals_above_rho": row["all_swap_cpu_intervals_above_rho"],
    } for row in results}, sort_keys=True))


if __name__ == "__main__":
    main()
