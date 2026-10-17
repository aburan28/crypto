#!/usr/bin/env python3
"""Independently replay the N83 same-base Bloom prefilter screen."""

import argparse
import importlib.util
import json
import random
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
GROUP_SOURCE = ROOT / "research/n83_indexed5_20261013/analyze.py"
RUNNER_SOURCE = ROOT / "research/n83_root11_bloom_20261015/run.py"
ANSWERS = ROOT / "research/n83_factor_base_refine_20261014/answers.jsonl"
ANSWERS_SHA256 = "f8cbc90f7f592aae5014ca4f2693a3327a72e935970f22e8f02868e4aeed7f71"
BASE_HASH = "a4ba33098d2776520b98a326585a2301f8ccd5bc005df7c1a70fbe8496c61d8b"
ROOT_KEYS = 903_852_075
BLOOM_BYTES = 1 << 29

group_spec = importlib.util.spec_from_file_location("n83_bloom_group", GROUP_SOURCE)
group = importlib.util.module_from_spec(group_spec)
group_spec.loader.exec_module(group)
runner_spec = importlib.util.spec_from_file_location("n83_bloom_runner", RUNNER_SOURCE)
runner = importlib.util.module_from_spec(runner_spec)
runner_spec.loader.exec_module(runner)


def read(path):
    return json.loads(path.read_text())


def read_lines(path):
    return [json.loads(line) for line in path.read_text().splitlines() if line]


def frozen_targets():
    runner.frozen_inputs()
    if runner.parent.digest(ANSWERS) != ANSWERS_SHA256:
        raise ValueError("frozen scalar sidecar changed")
    points = [group.point_from_json(row) for row in read_lines(runner.POINTS)]
    answers = read_lines(ANSWERS)
    if len(points) != 1000 or len(answers) != 1000 or len(set(points)) != 1000:
        raise ValueError("frozen target count or uniqueness differs")
    scalars = []
    for index, (point, answer) in enumerate(zip(points, answers)):
        scalar = int(answer["scalar"])
        if (answer["fixture_index"] != index
                or group.point_from_json(answer["point"]) != point
                or not group.on_curve(point)
                or group.multiply(point, group.R) is not None
                or group.multiply(group.G, scalar) != point):
            raise ValueError(f"frozen target {index} failed scalar replay")
        scalars.append(scalar)
    return points, scalars


def block_interval(values):
    if len(values) != 1000:
        raise ValueError("paired Bloom contrast requires 1000 targets")
    blocks = [statistics.mean(values[i:i + 50]) for i in range(0, 1000, 50)]
    randomizer = random.Random(2026101507)
    samples = [statistics.mean(blocks[randomizer.randrange(len(blocks))]
                               for _ in blocks) for _ in range(10_000)]
    samples.sort()
    return [samples[249], samples[9749]]


def verify_arm(output, arm, bloom, points, scalars):
    folder = output / "runs" / arm
    if not (folder / "usage.json").exists():
        return {"arm": arm, "status": "NOT_RUN"}
    usage = read(folder / "usage.json")
    if usage["status"] != "PASS_PRODUCER":
        return {"arm": arm, "status": usage["status"], "usage": usage}
    summary = read(folder / "stdout.jsonl")
    answers = read_lines(folder / "targets_result.jsonl")
    dump = read(folder / "base_dump.jsonl")
    expected_states = 83 * 3300 * 3299 // 2 + 42 * 3300
    if (summary["n"] != 83 or summary["a"] != 1
            or summary["targets"] != 1000
            or summary["targets_solved"] != 1000
            or summary["targets_failed"] != 0
            or summary["orbit_columns"] != 3300 or summary["rank"] != 3300
            or summary["root_table_layout"] != "indexed83_5byte_exact"
            or summary["root_table_slot_bytes"] != 5
            or summary["root_table_request_entries"] != 1 << 29
            or summary["root_table_allocated_bytes"] != 5 * (1 << 30)
            or summary["root_table_entries"] != ROOT_KEYS
            or summary["root_state_layout"] !=
            "root_only83_11byte_states_reconstructed_second_root"
            or summary["root_state_slot_bytes"] != 11
            or summary["regular_states"] != expected_states
            or summary["symmetric_s3_states_enabled"] is not True
            or summary["indexed5_fixed30_enabled"] is not True
            or summary["indexed5_load75_enabled"] is not False
            or summary["root_only11_enabled"] is not True
            or summary["root_bloom_enabled"] is not bloom
            or summary["root_bloom_bytes"] != (BLOOM_BYTES if bloom else 0)
            or summary["base_hash"] != BASE_HASH
            or len(answers) != 1000):
        raise ValueError(f"{arm}: index, Bloom, base or panel differs")
    base = group.replay_base(dump, summary)
    if len(base) != 547_800:
        raise ValueError(f"{arm}: factor-base size differs")
    times = []
    probes = []
    for index, (point, answer) in enumerate(zip(points, answers)):
        scalar = int(answer["recovered_scalar"])
        indices = answer["point_indices"]
        if (answer["fixture_index"] != index
                or group.point_from_json(answer["target"]) != point
                or scalar != scalars[index]
                or answer["group_verified"] is not True
                or group.multiply(group.G, scalar) != point
                or len(indices) != 4
                or any(not isinstance(i, int) or i < 0 or i >= len(base)
                       for i in indices)):
            raise ValueError(f"{arm}: target {index} scalar or indices failed")
        relation = None
        for base_index in indices:
            relation = group.add(relation, base[base_index])
        if relation != point:
            raise ValueError(f"{arm}: target {index} four-point relation failed")
        phase_total = sum(float(answer[name]) for name in
                          ("target_query_ms", "target_pdp_ms",
                           "target_relation_check_ms", "target_descent_ms",
                           "target_recovery_check_ms"))
        online = float(answer["target_ms"])
        if abs(phase_total - online) > 1:
            raise ValueError(f"{arm}: target {index} online phases differ")
        times.append(online)
        probes.append(answer["probes"])
    return {"arm": arm, "status": "PASS_INDEPENDENT_REPLAY",
            "base_hash": summary["base_hash"], "actual_factor_base_points": len(base),
            "effective_columns": 3300, "root_table_entries": summary["root_table_entries"],
            "root_bloom_bytes": summary["root_bloom_bytes"],
            "root_table_allocated_bytes": summary["root_table_allocated_bytes"],
            "root_state_allocated_bytes": summary["root_state_allocated_bytes"],
            "peak_rss_bytes": usage["peak_rss_bytes"],
            "resource_envelope_status": (
                "WITHIN_14_GIB" if usage["peak_rss_bytes"] <= runner.parent.RSS_CAP_BYTES
                else "ABOVE_14_GIB"),
            "setup_ms": summary["timing_ms"]["setup_total"],
            "mean_online_ms_exploratory": statistics.mean(times),
            "median_online_ms_exploratory": statistics.median(times),
            "mean_probes": statistics.mean(probes),
            "online_ms_per_target_exploratory": times}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    output = args.output.resolve()
    points, scalars = frozen_targets()
    arms = []
    for arm, bloom, _ in runner.ARMS:
        try:
            arms.append(verify_arm(output, arm, bloom, points, scalars))
        except Exception as error:
            arms.append({"arm": arm, "status": "VERIFICATION_FAILED",
                         "error": str(error)})
    by_arm = {row["arm"]: row for row in arms}
    usable = all(row["status"] == "PASS_INDEPENDENT_REPLAY"
                 and row["resource_envelope_status"] == "WITHIN_14_GIB"
                 for row in arms)
    same_base = (by_arm["bloom_off"].get("base_hash") ==
                 by_arm["bloom_on"].get("base_hash") == BASE_HASH
                 and by_arm["bloom_off"].get("root_table_entries") ==
                 by_arm["bloom_on"].get("root_table_entries") == ROOT_KEYS)
    contrast = None
    if usable and same_base:
        deltas = [enabled - control for enabled, control in zip(
            by_arm["bloom_on"]["online_ms_per_target_exploratory"],
            by_arm["bloom_off"]["online_ms_per_target_exploratory"])]
        interval = block_interval(deltas)
        contrast = {"definition": "bloom_on_minus_bloom_off_online_ms",
                    "mean_ms": statistics.mean(deltas),
                    "contiguous_block_size": 50, "resamples": 10_000,
                    "seed": 2026101507,
                    "bootstrap_95pct_ms": interval,
                    "confirmed_faster": interval[1] < 0}
    source = read(output / "source.json")
    if (source["schema"] != "n83.root11-bloom-source/v1"
            or source["hashes"] != runner.frozen_inputs()
            or source["screen_runner_sha256"] != runner.parent.digest(RUNNER_SOURCE)):
        raise ValueError("screen source or inputs differ")
    analysis = {"schema": "n83.root11-bloom-analysis/v1",
                "hardware_tier": "cloud_instance_exploratory",
                "workload_role": "secondary_reused_1000_point_layout_screen",
                "target_count": len(points), "arms": arms,
                "same_base_and_key_count": same_base,
                "paired_difference": contrast,
                "selected_arm_exploratory": (
                    "bloom_on" if contrast and contrast["confirmed_faster"]
                    else "bloom_off" if usable and same_base else None),
                "status": "COMPLETE" if usable and same_base else "INCOMPLETE",
                "controlled_speedup": None}
    (output / "analysis.json").write_text(json.dumps(analysis, indent=2,
                                                    sort_keys=True) + "\n")
    print(json.dumps(analysis, sort_keys=True))
    return 0 if analysis["status"] == "COMPLETE" else 1


if __name__ == "__main__":
    raise SystemExit(main())
