#!/usr/bin/env python3
"""Independently replay every completed 4000-target N83 K3300 layout arm."""

import argparse
import importlib.util
import json
import random
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "research/n83_indexed5_20261013/analyze.py"
RUNNER = ROOT / "research/n83_root4_vs_root11_20261015/run.py"
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
        raise ValueError("confirmation requires 4000 paired targets")
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
    expected_bases = {
        "k3300_root4": "a4ba33098d2776520b98a326585a2301f8ccd5bc005df7c1a70fbe8496c61d8b",
        "k3300_root11": "a4ba33098d2776520b98a326585a2301f8ccd5bc005df7c1a70fbe8496c61d8b",
    }
    base_matches = {}
    for arm, expected in expected_bases.items():
        row = controls[arm]
        base_matches[arm] = (row["base_hash"] == expected
                             if row["status"] == "PASS_INDEPENDENT_REPLAY" else None)
        if base_matches[arm] is False:
            row["status"] = "FROZEN_BASE_MISMATCH"
    same_base_and_key_count = None
    if all(row["status"] == "PASS_INDEPENDENT_REPLAY" for row in arms):
        root4, root11 = controls["k3300_root4"], controls["k3300_root11"]
        same_base_and_key_count = (
            root4["base_hash"] == root11["base_hash"]
            and root4["root_table_entries"] == root11["root_table_entries"])
        if not same_base_and_key_count:
            root4["status"] = "CONTROL_PARITY_FAILED"
    eligible = [row for row in arms if row["status"] == "PASS_INDEPENDENT_REPLAY"
                and row["resource_envelope_status"] == "WITHIN_14_GIB"]
    selected = (min(eligible, key=lambda row: row["mean_online_ms_exploratory"])
                if eligible else None)
    paired_difference = None
    if len(eligible) == 2:
        a, b = controls["k3300_root4"], controls["k3300_root11"]
        deltas = [b_row["online_ms_exploratory"] - a_row["online_ms_exploratory"]
                  for a_row, b_row in zip(a["targets"], b["targets"])]
        paired_difference = {
            "contrast": "k3300_root11_minus_k3300_root4",
            "mean_ms": statistics.mean(deltas),
            "contiguous_block_size": 50,
            "contiguous_block_bootstrap_95pct_ms":
                target_block_mean_interval(deltas, 2026101599),
        }
    interval = (paired_difference["contiguous_block_bootstrap_95pct_ms"]
                if paired_difference else None)
    confirmation = ("PASS_ROOT4_FASTER" if interval and interval[0] > 0 else
                    "PASS_ROOT11_FASTER" if interval and interval[1] < 0 else
                    "INCONCLUSIVE" if interval else "INCOMPLETE")
    analysis = {"schema": "n83.root4-vs-root11-analysis/v1",
                "hardware_tier": "cloud_instance_exploratory",
                "workload_role": "secondary_layout_confirmation",
                "target_count": len(points), "arms": arms,
                "frozen_base_matches": base_matches,
                "same_base_and_key_count": same_base_and_key_count,
                "selection_rule": "predeclared root11-minus-root4 paired online mean with 50-target contiguous-block bootstrap; lower 95% endpoint above zero confirms root4",
                "selected_arm_descriptive": selected["arm"] if selected else None,
                "paired_difference": paired_difference,
                "confirmation_status": confirmation,
                "bootstrap_resamples": 10_000,
                "controlled_speedup": None}
    (output / "analysis.json").write_text(json.dumps(analysis, indent=2,
                                                    sort_keys=True) + "\n")
    print(json.dumps(analysis, sort_keys=True))
    return 0 if (selected is not None and same_base_and_key_count is True
                 and all(value is True for value in base_matches.values()) and all(
        row["status"] == "PASS_INDEPENDENT_REPLAY" for row in arms)) else 1


if __name__ == "__main__":
    raise SystemExit(main())
