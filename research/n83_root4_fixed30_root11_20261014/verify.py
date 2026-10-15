#!/usr/bin/env python3
"""Independently replay the N83 four-byte/11-byte capacity panel."""

import argparse
import importlib.util
import json
import random
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "research/n83_indexed5_20261013/analyze.py"
RUNNER = ROOT / "research/n83_root4_fixed30_root11_20261014/run.py"
HOLDOUT = ROOT / "research/n83_root4_fixed30_20261014/targets.jsonl"
ANSWERS = ROOT / "research/n83_root4_fixed30_20261014/answers.jsonl"
ANSWERS_SHA256 = "f8cbc90f7f592aae5014ca4f2693a3327a72e935970f22e8f02868e4aeed7f71"

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
    if runner.base.parent.digest(ANSWERS) != ANSWERS_SHA256:
        raise ValueError("1000-point answer sidecar changed")
    points = [group.point_from_json(row) for row in json_lines(HOLDOUT)]
    answers = json_lines(ANSWERS)
    if len(points) != len(answers) or len(points) != 1000 or len(set(points)) != len(points):
        raise ValueError("1000-point target panel is incomplete or has duplicate points")
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
    rng = random.Random(seed)
    samples = [statistics.mean(values[rng.randrange(len(values))] for _ in values)
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
    layout = {"root4f30r11": "root_only83_11byte_states_reconstructed_second_root",
              "root4f30": "root_only83_21byte_states",
              "root11": "root_only83_11byte_states_reconstructed_second_root",
              "root21": "root_only83_21byte_states",
              "root22": "root_only83_22byte_states"}[variant]
    state_bytes = {"root4f30r11": 11, "root4f30": 21, "root11": 11,
                   "root21": 21, "root22": 22}[variant]
    four = variant.startswith("root4f30")
    fixed5 = variant == "root11"
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
            or summary["indexed4_fixed30_enabled"] is not four
            or summary["indexed4_fixed31_enabled"] is not False
            or summary["root_only11_enabled"] is not
            (variant in ("root11", "root4f30r11"))
            or len(target_rows) != len(points)):
        raise ValueError(f"{arm}: source layout, base, rank or panel size mismatch")
    if (summary["root_table_entries"] >=
            summary["root_table_allocated_bytes"] // slot_bytes):
        raise ValueError(f"{arm}: exact table occupancy invalid")
    if four and (summary["root_table_request_entries"] != 1 << 29
                 or summary["root_table_allocated_bytes"] != 4 * (1 << 30)):
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
                "WITHIN_14_GIB" if usage["peak_rss_bytes"] <= runner.base.parent.RSS_CAP_BYTES
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
    parity = {}
    for baseline in ("root11", "root4f30"):
        a, b = controls["k3300_root4f30r11"], controls[f"k3300_{baseline}"]
        equal = None
        if a["status"] == b["status"] == "PASS_INDEPENDENT_REPLAY":
            equal = (a["base_hash"] == b["base_hash"] and
                     a["root_table_entries"] == b["root_table_entries"] and
                     [t["probes"] for t in a["targets"]] ==
                     [t["probes"] for t in b["targets"]])
            if not equal:
                a["status"] = "CONTROL_PARITY_FAILED"
        parity[f"k3300_root4f30r11_vs_{baseline}"] = equal
    eligible = [row for row in arms if row["status"] == "PASS_INDEPENDENT_REPLAY"
                and row["resource_envelope_status"] == "WITHIN_14_GIB"]
    selected = (min(eligible, key=lambda row: row["mean_online_ms_exploratory"])
                if eligible else None)
    runner_up = (sorted(eligible, key=lambda row: row["mean_online_ms_exploratory"])[1]
                 if len(eligible) >= 2 else None)
    paired_difference = None
    if selected and runner_up:
        deltas = [b["online_ms_exploratory"] - a["online_ms_exploratory"]
                  for a, b in zip(selected["targets"], runner_up["targets"])]
        paired_difference = {
            "runner_up_arm": runner_up["arm"],
            "mean_runner_up_minus_selected_ms": statistics.mean(deltas),
            "target_block_bootstrap_95pct_ms":
                target_block_mean_interval(deltas, 2026101499),
        }
    analysis = {"schema": "n83.root4-fixed30-root11-analysis/v1",
                "hardware_tier": "cloud_instance_exploratory",
                "workload_role": "secondary_factor_base_selection",
                "target_count": len(points), "arms": arms,
                "same_k_base_key_probe_parity": parity,
                "selection_rule": "lowest arithmetic mean of 1000 verified online intervals among arms with observed peak <=14 GiB",
                "selected_arm_descriptive": selected["arm"] if selected else None,
                "selected_vs_runner_up_descriptive": paired_difference,
                "bootstrap_resamples": 10_000,
                "controlled_speedup": None}
    (output / "analysis.json").write_text(json.dumps(analysis, indent=2,
                                                    sort_keys=True) + "\n")
    print(json.dumps(analysis, sort_keys=True))
    return 0 if (selected is not None and all(value is True for value in parity.values()) and all(
        row["status"] == "PASS_INDEPENDENT_REPLAY" for row in arms)) else 1


if __name__ == "__main__":
    raise SystemExit(main())
