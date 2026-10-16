#!/usr/bin/env python3
"""Independently replay six balanced-order N83 root Bloom profile arms."""

import argparse
import importlib.util
import json
import random
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
SOURCE = ROOT / "research/n83_indexed5_20261013/analyze.py"
RUNNER = HERE / "run.py"
HOLDOUT = ROOT / "research/n83_k_abba_20261015/targets.jsonl"
ANSWERS = ROOT / "research/n83_k_abba_20261015/answers.jsonl"
ANSWERS_SHA256 = "6e1181ba958708803316903a645463425bc79130074041d39c8094778e2f4379"
BASE_HASH = "a4ba33098d2776520b98a326585a2301f8ccd5bc005df7c1a70fbe8496c61d8b"
ROOT_KEYS = 903852075
PROFILE = {"standard": (3, 8, 536870912),
           "wide3": (3, 16, 1073741824),
           "wide5": (5, 16, 1073741824)}

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


def cgroup_counter(receipt, section, key):
    lines = (receipt.get(section) or "").splitlines()
    for line in lines:
        fields = line.split()
        if len(fields) == 2 and fields[0] == key:
            return int(fields[1])
    return None


def frozen_targets():
    runner.frozen_inputs()
    if runner.parent.digest(ANSWERS) != ANSWERS_SHA256:
        raise ValueError("4000-point answer sidecar changed")
    points = [group.point_from_json(row) for row in json_lines(HOLDOUT)]
    answers = json_lines(ANSWERS)
    if len(points) != len(answers) or len(points) != 4000 or len(set(points)) != len(points):
        raise ValueError("4000-point target panel is incomplete or has duplicate points")
    if not group.on_curve(group.G) or group.multiply(group.G, group.R) is not None:
        raise ValueError("frozen subgroup generator or order changed")
    scalars = []
    for index, (point, answer) in enumerate(zip(points, answers)):
        scalar = int(answer["scalar"])
        if answer["fixture_index"] != index or group.point_from_json(answer["point"]) != point:
            raise ValueError("target and answer sidecar pairing differs")
        if not group.on_curve(point) or group.multiply(group.G, scalar) != point:
            raise ValueError("archived target scalar or subgroup verification failed")
        scalars.append(scalar)
    return points, scalars


def target_block_mean_interval(values, seed):
    if len(values) != 4000:
        raise ValueError("balanced Bloom comparison requires 4000 paired targets")
    block_means = [statistics.mean(values[i:i + 50])
                   for i in range(0, len(values), 50)]
    rng = random.Random(seed)
    samples = [statistics.mean(block_means[rng.randrange(len(block_means))]
                               for _ in block_means)
               for _ in range(10_000)]
    samples.sort()
    return [samples[49], samples[9949]]


def verify_arm(output, arm, profile, points, scalars):
    k, variant = 3300, "root11"
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
            or summary["sparse_canonical_enabled"] is not True
            or summary["root_bloom_enabled"] is not True
            or (summary["root_bloom_hashes"],
                summary["root_bloom_bits_per_request"],
                summary["root_bloom_bytes"]) != PROFILE[profile]
            or summary["root_bloom_profile"] != profile
            or summary["base_hash"] != BASE_HASH
            or summary["root_table_entries"] != ROOT_KEYS
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
                or answer["group_verified"] is not True):
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
    before, after = usage["cgroup_before"], usage["cgroup_after"]
    throttled_before = cgroup_counter(before, "cpu.stat", "nr_throttled")
    throttled_after = cgroup_counter(after, "cpu.stat", "nr_throttled")
    oom_before = cgroup_counter(before, "memory.events", "oom_kill")
    oom_after = cgroup_counter(after, "memory.events", "oom_kill")
    throttled = (throttled_after - throttled_before
                 if throttled_before is not None and throttled_after is not None else None)
    oom_kills = (oom_after - oom_before
                 if oom_before is not None and oom_after is not None else None)
    return {"arm": arm, "profile": profile, "status": "PASS_INDEPENDENT_REPLAY",
            "actual_factor_base_points": len(base), "effective_columns": k,
            "base_hash": summary["base_hash"],
            "root_table_entries": summary["root_table_entries"],
            "root_table_allocated_bytes": summary["root_table_allocated_bytes"],
            "root_state_allocated_bytes": summary["root_state_allocated_bytes"],
            "peak_rss_bytes": usage["peak_rss_bytes"],
            "cpu_throttled_periods_delta": throttled,
            "oom_kills_delta": oom_kills,
            "resource_envelope_status": (
                "WITHIN_14_GIB" if usage["peak_rss_bytes"] <= runner.parent.RSS_CAP_BYTES
                else "ABOVE_14_GIB"),
            "setup_ms": summary["timing_ms"]["setup_total"],
            "mean_online_ms_exploratory": statistics.mean(times),
            "mean_online_target_block_bootstrap_99pct_ms":
                target_block_mean_interval(times, 2026101550 + runner.ARMS.index(next(a for a in runner.ARMS if a[0] == arm))),
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
    for arm, profile, _ in runner.ARMS:
        try:
            arms.append(verify_arm(output, arm, profile, points, scalars))
        except Exception as error:
            arms.append({"arm": arm,
                         "status": "VERIFICATION_FAILED", "error": str(error)})
    eligible = [row for row in arms if row["status"] == "PASS_INDEPENDENT_REPLAY"
                and row["resource_envelope_status"] == "WITHIN_14_GIB"
                and row["cpu_throttled_periods_delta"] == 0
                and row["oom_kills_delta"] == 0]
    eligible_by_arm = {row["arm"]: row for row in eligible}
    complete = len(eligible) == len(runner.ARMS)
    balanced = None
    selected_profile = None
    if complete:
        online = {arm: [target["online_ms_exploratory"] for target in row["targets"]]
                  for arm, row in eligible_by_arm.items()}
        profiles = ("standard", "wide3", "wide5")
        means = {profile: [(online[f"{profile}_a"][i] + online[f"{profile}_b"][i]) / 2
                           for i in range(len(points))] for profile in profiles}
        pairwise = {}
        for i, candidate in enumerate(profiles):
            for j, other in enumerate(profiles):
                if candidate == other:
                    continue
                delta = [left - right for left, right in zip(means[candidate], means[other])]
                pairwise[f"{candidate}_minus_{other}"] = {
                    "mean_ms": statistics.mean(delta),
                    "contiguous_block_bootstrap_99pct_ms":
                        target_block_mean_interval(delta, 2026101540 + i * 3 + j),
                }
        fastest = min(profiles, key=lambda profile: statistics.mean(means[profile]))
        if all(pairwise[f"{fastest}_minus_{other}"]["contiguous_block_bootstrap_99pct_ms"][1] < 0
               for other in profiles if other != fastest):
            selected_profile = fastest
        repeatability = {}
        for index, profile in enumerate(profiles):
            same_profile = [online[f"{profile}_a"][i] - online[f"{profile}_b"][i]
                            for i in range(len(points))]
            repeatability[f"{profile}_a_minus_b"] = {
                "mean_ms": statistics.mean(same_profile),
                "contiguous_block_bootstrap_99pct_ms":
                    target_block_mean_interval(same_profile, 2026101560 + index),
            }
        balanced = {
            "profile_mean_ms": {profile: statistics.mean(means[profile]) for profile in profiles},
            "contiguous_block_size": 50,
            "pairwise": pairwise,
            "repeatability": repeatability,
            "selected_profile": selected_profile,
        }
    source = one_json(output / "source.json")
    if (source["schema"] != "n83.bloom-profile-source/v1"
            or source["hashes"] != runner.frozen_inputs()
            or source["arm_order"] != [arm for arm, _, _ in runner.ARMS]
            or source["panel_runner_sha256"] != runner.parent.digest(RUNNER)):
        raise ValueError("executed Bloom-profile source or fixture changed")
    analysis = {"schema": "n83.bloom-profile-analysis/v1",
                "hardware_tier": "cloud_instance_exploratory",
                "workload_role": "secondary_reused_4000_point_balanced_order_bloom_profile_comparison",
                "target_count": len(points), "arms": arms,
                "arm_order": source["arm_order"],
                "selection_rule": "standard/wide3/wide5/wide5/wide3/standard with fresh empty caches; average two runs per profile on each public point; 80 contiguous 50-target blocks, 10000 resamples and fixed pair seeds 2026101540+i*3+j; select the lowest-mean profile only when its upper 99% paired endpoint is below zero versus both others",
                "balanced_comparison": balanced,
                "selected_profile": selected_profile,
                "status": "COMPLETE" if complete else "INCOMPLETE",
                "bootstrap_resamples": 10_000,
                "controlled_speedup": None}
    (output / "analysis.json").write_text(json.dumps(analysis, indent=2,
                                                    sort_keys=True) + "\n")
    print(json.dumps(analysis, sort_keys=True))
    return 0 if complete else 1


if __name__ == "__main__":
    raise SystemExit(main())
