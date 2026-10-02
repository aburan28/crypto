#!/usr/bin/env python3
"""Apply the preregistered four-window stop rule to verified charged CPU."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import statistics

from run_screen import config, schedule
from verify_screen import paired


def cpu(child: dict) -> float:
    assert child["exit_code"] == 0 and child["stopped_for"] is None
    value = child["child_user_cpu_seconds"] + child["child_system_cpu_seconds"]
    assert math.isfinite(value) and value > 0
    return value


def same_numbers(a: object, b: object) -> bool:
    if isinstance(a, dict) and isinstance(b, dict):
        return set(a) == set(b) and all(same_numbers(a[key], b[key]) for key in a)
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and all(same_numbers(x, y) for x, y in zip(a, b))
    if isinstance(a, (float, int)) and isinstance(b, (float, int)):
        return math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-12)
    return a == b


def analyze(report: dict, replay: dict) -> dict:
    cfg = config()
    assert report["schema"] == "ecc2k130-base-window-screen-run-v1"
    assert report["config"] == cfg
    assert replay["schema"] == "ecc2k130-base-window-screen-replay-v1"
    assert report["plan"] == [
        {"block": b, "arm": arm} for b, arm in schedule(cfg)
    ]
    base = {
        "schema": "ecc2k130-base-window-screen-analysis-v1",
        "classification": cfg["classification"],
        "curve_slug": cfg["curve_slug"],
        "n": cfg["n"],
        "L": cfg["public_targets_per_block"],
        "K": cfg["useful_orbit_columns"],
        "blocks": cfg["blocks"],
        "children_planned": 45,
        "timing_eligible": False,
        "method_speedup_claim": None,
        "n131_transfer_claim": None,
    }
    if replay["status"] != "PASS" or report["status"] != "PASS":
        return {
            **base,
            "decision": "CENSORED",
            "reason": "producer_or_independent_replay_failure",
            "producer_status": report["status"],
            "replay_status": replay["status"],
            "failure": replay.get("failure") or report.get("failure"),
            "children_recorded": len(report["runs"]),
        }
    assert len(report["runs"]) == 45
    by_block = {block: {} for block in range(cfg["blocks"])}
    for item in report["runs"]:
        block, arm = item["block"], item["arm"]
        children = item["children"]
        if arm == "rho":
            assert set(children) == {"rho"}
            total = cpu(children["rho"])
            phases = {"rho_cpu_seconds": total}
        else:
            assert set(children) == {"generator", "compact"}
            generation, compact = cpu(children["generator"]), cpu(children["compact"])
            total = generation + compact
            phases = {
                "generator_cpu_seconds": generation,
                "compact_cpu_seconds": compact,
            }
        by_block[block][arm] = {"charged_cpu_seconds": total, **phases}
    assert all(set(arms) == set(cfg["arm_order_before_rotation"])
               for arms in by_block.values())
    aa = {}
    fixed = {}
    paired_cpu = {str(window): [] for window in range(4)}
    rho_cpu = []
    ratios_by_block = {block: {} for block in range(cfg["blocks"])}
    for block in range(cfg["blocks"]):
        arms = by_block[block]
        reference = arms["rho"]["charged_cpu_seconds"]
        rho_cpu.append(reference)
        for window in range(4):
            a = arms[f"w{window}_a"]["charged_cpu_seconds"]
            b = arms[f"w{window}_b"]["charged_cpu_seconds"]
            candidate = math.sqrt(a * b)
            paired_cpu[str(window)].append(candidate)
            ratios_by_block[block][str(window)] = candidate / reference
    for window in range(4):
        key = str(window)
        aa_values = [
            by_block[block][f"w{window}_b"]["charged_cpu_seconds"] /
            by_block[block][f"w{window}_a"]["charged_cpu_seconds"]
            for block in range(cfg["blocks"])
        ]
        fixed_values = [ratios_by_block[block][key] for block in range(cfg["blocks"])]
        aa[key] = paired(aa_values, cfg["log_t_critical_df4"])
        fixed[key] = paired(fixed_values, cfg["log_t_critical_df4"])
        assert same_numbers(aa[key], replay["aa"][key])
        assert same_numbers(fixed[key], replay["fixed_window_over_rho"][key])
    assert replay["independent_rank_replays"] == 40
    assert replay["verified_compact_target_logs"] == 40960
    assert replay["verified_rho_target_logs"] == 5120
    aa_valid = {
        key: cfg["aa_ratio_min"] <= row["median"] <= cfg["aa_ratio_max"]
        and row["interval_95pct"][0] <= 1 <= row["interval_95pct"][1]
        for key, row in aa.items()
    }
    assert aa_valid == replay["aa_valid"]
    eligible = replay["uncontended"] and all(aa_valid.values())
    assert replay["timing_eligible"] == eligible
    if not eligible:
        return {
            **base,
            "decision": "CENSORED",
            "reason": "isolation_or_aa_gate",
            "children_recorded": 45,
            "uncontended": replay["uncontended"],
            "aa_valid": aa_valid,
            "aa": aa,
            "fixed_window_over_rho": fixed,
        }
    hindsight = [
        min(ratios_by_block[block].values()) for block in range(cfg["blocks"])
    ]
    h_stats = paired(hindsight, cfg["log_t_critical_df4"])
    crossing_fixed = [
        window for window in range(4)
        if fixed[str(window)]["interval_95pct"][1] < cfg["parity_ratio"]
    ]
    if all(value > 1 for value in hindsight) and h_stats["interval_95pct"][0] > 1:
        decision = "NO_WINDOW_OPPORTUNITY"
    elif crossing_fixed:
        decision = "FIXED_WINDOW_CROSSOVER_CANDIDATE"
    elif statistics.median(hindsight) < 1:
        decision = "HINDSIGHT_VARIATION_ONLY"
    else:
        decision = "INCONCLUSIVE"
    return {
        **base,
        "timing_eligible": True,
        "decision": decision,
        "children_recorded": 45,
        "independent_rank_replays": 40,
        "verified_compact_target_logs": 40960,
        "verified_rho_target_logs": 5120,
        "charged_cpu_seconds_by_block_and_arm": {
            str(block): by_block[block] for block in range(cfg["blocks"])
        },
        "rho_cpu_seconds_by_block": rho_cpu,
        "paired_window_cpu_seconds_by_block": paired_cpu,
        "fixed_window_over_rho": fixed,
        "aa": aa,
        "hindsight_lower_envelope": {
            **h_stats,
            "role": "post_hoc_optimistic_lower_bound_not_a_policy",
        },
        "crossing_fixed_windows_requiring_confirmation": crossing_fixed,
        "selection_cpu_charged": None,
        "common_operation_equivalent_S": None,
        "second_host_and_new_Q_confirmation": None,
        "next_direction_if_no_opportunity": (
            "natural high-arity descendant-native PDP yield and non-eager index/query policy"
            if decision == "NO_WINDOW_OPPORTUNITY" else None
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--replay", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    assert not args.out.exists(), "never overwrite an analysis"
    report = json.loads((args.run_dir / "screen_run.json").read_text())
    replay = json.loads(args.replay.read_text())
    result = analyze(report, replay)
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "decision": result["decision"],
        "timing_eligible": result["timing_eligible"],
        "children_recorded": result["children_recorded"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
