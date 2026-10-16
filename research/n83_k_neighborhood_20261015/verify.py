#!/usr/bin/env python3
"""Independently replay the bounded N83 root11 K-neighborhood."""

import argparse
import importlib.util
import json
import random
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "research/n83_indexed5_20261013/analyze.py"
RUNNER = ROOT / "research/n83_k_neighborhood_20261015/run.py"
HOLDOUT = ROOT / "research/n83_factor_base_confirm_20261015/targets.jsonl"
ANSWERS = ROOT / "research/n83_factor_base_confirm_20261015/answers.jsonl"
ANSWERS_SHA256 = "8e18c9c97ba9f57e5eb43f1b39b9e9575b5cfab87840c25c428dac150c868101"

spec = importlib.util.spec_from_file_location("n83_panel_group", SOURCE)
group = importlib.util.module_from_spec(spec)
spec.loader.exec_module(group)
runner_spec = importlib.util.spec_from_file_location("n83_panel_runner", RUNNER)
runner = importlib.util.module_from_spec(runner_spec)
runner_spec.loader.exec_module(runner)


def one_json(path):
    return json.loads(path.read_text())


def json_lines(path):
    return [json.loads(line) for line in path.read_text().splitlines() if line]


def frozen_targets():
    runner.frozen_inputs()
    if runner.parent.digest(ANSWERS) != ANSWERS_SHA256:
        raise ValueError("4000-point answer sidecar changed")
    points = [group.point_from_json(row) for row in json_lines(HOLDOUT)]
    answers = json_lines(ANSWERS)
    if len(points) != len(answers) or len(points) != 4000 or len(set(points)) != len(points):
        raise ValueError("4000-point target panel is incomplete or has duplicate points")
    scalars = []
    for index, (point, answer) in enumerate(zip(points, answers)):
        scalar = int(answer["scalar"])
        if answer["fixture_index"] != index or group.point_from_json(answer["point"]) != point:
            raise ValueError("target and answer sidecar pairing differs")
        if (not group.on_curve(point) or group.multiply(point, group.R) is not None
                or group.multiply(group.G, scalar) != point):
            raise ValueError("archived target scalar or subgroup verification failed")
        scalars.append(scalar)
    return points, scalars


def target_block_mean_interval(values, seed):
    if len(values) != 4000:
        raise ValueError("neighborhood comparison requires 4000 paired targets")
    block_means = [statistics.mean(values[i:i + 50])
                   for i in range(0, len(values), 50)]
    rng = random.Random(seed)
    samples = [statistics.mean(block_means[rng.randrange(len(block_means))]
                               for _ in block_means)
               for _ in range(10_000)]
    samples.sort()
    return [samples[249], samples[9749]]


def verify_arm(output, k, variant, points, scalars):
    arm = f"k{k}_{variant}"
    folder = output / "runs" / arm
    if not (folder / "usage.json").exists():
        return {"arm": arm, "status": "NOT_RUN"}
    usage = one_json(folder / "usage.json")
    if usage["status"] != "PASS_PRODUCER":
        return {"arm": arm, "status": usage["status"], "usage": usage}
    summary = one_json(folder / "stdout.jsonl")
    target_rows = json_lines(folder / "targets_result.jsonl")
    dump = one_json(folder / "base_dump.jsonl")
    layout = {"root4": "root_only83_11byte_states_reconstructed_second_root",
              "root5": "root_only83_11byte_states_reconstructed_second_root",
              "root11": "root_only83_11byte_states_reconstructed_second_root",
              "root21": "root_only83_21byte_states",
              "root22": "root_only83_22byte_states"}[variant]
    state_bytes = {"root4": 11, "root5": 11, "root11": 11,
                   "root21": 21, "root22": 22}[variant]
    four = variant == "root4"
    fixed5 = variant in ("root5", "root11", "root21")
    slot_bytes = 4 if four else 5
    state_bound = 83 * k * (k - 1) // 2 + 42 * k
    if (summary["n"] != 83 or summary["a"] != 1
            or summary["targets"] != len(points)
            or summary["targets_solved"] != len(points)
            or summary["targets_failed"] != 0
            or summary["orbit_columns"] != k or summary["rank"] != k
            or summary["root_table_layout"] !=
            ("indexed83_4byte_exact" if four else "indexed83_5byte_exact")
            or summary["root_table_slot_bytes"] != slot_bytes
            or summary["root_table_indexed4_enabled"] is not four
            or summary["root_table_indexed5_enabled"] is not (not four)
            or summary["root_state_layout"] != layout
            or summary["root_state_slot_bytes"] != state_bytes
            or summary["regular_states"] != state_bound
            or summary["symmetric_s3_states_enabled"] is not True
            or summary["indexed5_load75_enabled"] is not (variant == "root22")
            or summary["indexed5_fixed30_enabled"] is not fixed5
            or summary["indexed4_fixed31_enabled"] is not four
            or summary["root_only11_enabled"] is not
            (variant in ("root4", "root5", "root11"))
            or len(target_rows) != len(points)):
        raise ValueError(f"{arm}: source layout, base, rank or panel size mismatch")
    if (summary["root_table_entries"] >=
            summary["root_table_allocated_bytes"] // slot_bytes):
        raise ValueError(f"{arm}: exact table occupancy invalid")
    if four and (summary["root_table_request_entries"] != 1 << 30
                 or summary["root_table_allocated_bytes"] != 4 * (1 << 31)):
        raise ValueError(f"{arm}: fixed four-byte table allocation mismatch")
    if fixed5 and (summary["root_table_request_entries"] != 1 << 29
                   or summary["root_table_allocated_bytes"] != 5 * (1 << 30)):
        raise ValueError(f"{arm}: fixed table allocation mismatch")
    base = group.replay_base(dump, summary)
    if len(base) != 166 * k:
        raise ValueError(f"{arm}: actual factor-base size mismatch")
    targets = []
    for index, (point, answer) in enumerate(zip(points, target_rows)):
        scalar = int(answer["recovered_scalar"])
        if (answer["fixture_index"] != index
                or group.point_from_json(answer["target"]) != point
                or scalar != scalars[index]
                or answer["group_verified"] is not True
                or group.multiply(group.G, scalar) != point):
            raise ValueError(f"{arm}: target {index} scalar or point replay failed")
        indices = answer["point_indices"]
        if (len(indices) != 4 or any(not isinstance(i, int) or i < 0 or i >= len(base)
                                     for i in indices)):
            raise ValueError(f"{arm}: target {index} invalid four-point indices")
        relation = None
        for base_index in indices:
            relation = group.add(relation, base[base_index])
        if relation != point:
            raise ValueError(f"{arm}: target {index} four-point relation failed")
        phase_ms = sum(float(answer[key]) for key in
                       ("target_query_ms", "target_pdp_ms",
                        "target_relation_check_ms", "target_descent_ms",
                        "target_recovery_check_ms"))
        online_ms = float(answer["target_ms"])
        if abs(phase_ms - online_ms) > 1.0:
            raise ValueError(f"{arm}: target {index} exclusive online phases differ")
        targets.append({"index": index, "online_ms_exploratory": online_ms,
                        "probes": answer["probes"], "scalar_replay": "PASS",
                        "relation_replay": "PASS"})
    times = [row["online_ms_exploratory"] for row in targets]
    return {"arm": arm, "status": "PASS_INDEPENDENT_REPLAY",
            "actual_factor_base_points": len(base), "effective_columns": k,
            "base_hash": summary["base_hash"],
            "root_table_entries": summary["root_table_entries"],
            "root_table_allocated_bytes": summary["root_table_allocated_bytes"],
            "root_state_allocated_bytes": summary["root_state_allocated_bytes"],
            "peak_rss_bytes": usage["peak_rss_bytes"],
            "resource_envelope_status": (
                "WITHIN_14_GIB" if usage["peak_rss_bytes"] <= runner.parent.RSS_CAP_BYTES
                else "ABOVE_14_GIB"),
            "setup_ms": summary["timing_ms"]["setup_total"],
            "mean_online_ms_exploratory": statistics.mean(times),
            "mean_online_target_block_bootstrap_95pct_ms":
                target_block_mean_interval(times, 2026101400 + k),
            "median_online_ms_exploratory": statistics.median(times),
            "max_online_ms_exploratory": max(times),
            "targets": targets}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    output = args.output.resolve()
    points, scalars = frozen_targets()
    arms = []
    for k, variant, _ in runner.ARMS:
        try:
            arms.append(verify_arm(output, k, variant, points, scalars))
        except Exception as error:
            arms.append({"arm": f"k{k}_{variant}",
                         "status": "VERIFICATION_FAILED", "error": str(error)})
    controls = {row["arm"]: row for row in arms}
    baseline = controls["k3300_root11"]
    baseline_hash_ok = (baseline["base_hash"] ==
                        "a4ba33098d2776520b98a326585a2301f8ccd5bc005df7c1a70fbe8496c61d8b"
                        if baseline["status"] == "PASS_INDEPENDENT_REPLAY" else None)
    if baseline_hash_ok is False:
        baseline["status"] = "FROZEN_BASE_MISMATCH"
    eligible = [row for row in arms if row["status"] == "PASS_INDEPENDENT_REPLAY"
                and row["resource_envelope_status"] == "WITHIN_14_GIB"]
    eligible_by_arm = {row["arm"]: row for row in eligible}
    comparisons = {}
    if "k3300_root11" in eligible_by_arm:
        base_targets = baseline["targets"]
        for k in (3200, 3250, 3350):
            arm = f"k{k}_root11"
            if arm not in eligible_by_arm:
                comparisons[arm] = None
                continue
            deltas = [other["online_ms_exploratory"] - base["online_ms_exploratory"]
                      for base, other in zip(base_targets,
                                             eligible_by_arm[arm]["targets"])]
            interval = target_block_mean_interval(deltas, 2026101599 + k)
            comparisons[arm] = {
                "contrast": f"{arm}_minus_k3300_root11",
                "mean_ms": statistics.mean(deltas),
                "contiguous_block_size": 50,
                "contiguous_block_bootstrap_95pct_ms": interval,
                "confirmed_faster_than_k3300": interval[1] < 0,
            }
    confirmed = [eligible_by_arm[arm] for arm, result in comparisons.items()
                 if result is not None and result["confirmed_faster_than_k3300"]]
    selected = (min(confirmed, key=lambda row: row["mean_online_ms_exploratory"])
                if confirmed else eligible_by_arm.get("k3300_root11"))
    complete = len(eligible) == len(runner.ARMS) and baseline_hash_ok is True
    source = one_json(output / "source.json")
    if (source["schema"] != "n83.k-neighborhood-source/v1"
            or source["hashes"] != runner.frozen_inputs()
            or source["panel_runner_sha256"] != runner.parent.digest(RUNNER)):
        raise ValueError("executed neighborhood source or fixture changed")
    analysis = {"schema": "n83.k-neighborhood-analysis/v1",
                "hardware_tier": "cloud_instance_exploratory",
                "workload_role": "secondary_k_neighborhood_reusing_confirmation_points",
                "target_count": len(points), "arms": arms,
                "baseline_base_hash_matches": baseline_hash_ok,
                "selection_rule": "K3200/K3250/K3350 minus K3300 paired online means; 50-target contiguous-block bootstrap, 10000 resamples, seed 2026101599+K; promote fastest arm whose upper 95% endpoint is below zero",
                "selected_arm_confirmed": selected["arm"] if selected else None,
                "paired_comparisons": comparisons,
                "status": "COMPLETE" if complete else "INCOMPLETE",
                "bootstrap_resamples": 10_000,
                "controlled_speedup": None}
    (output / "analysis.json").write_text(json.dumps(analysis, indent=2,
                                                    sort_keys=True) + "\n")
    print(json.dumps(analysis, sort_keys=True))
    return 0 if complete else 1


if __name__ == "__main__":
    raise SystemExit(main())
