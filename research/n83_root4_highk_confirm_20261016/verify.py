#!/usr/bin/env python3
"""Independently replay fresh N83 K3500/K3596 confirmation."""

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
HOLDOUT = HERE / "targets.jsonl"
ANSWERS = HERE / "answers.jsonl"
ANSWERS_SHA256 = "da9cd25bfe87eedd2a2c74d39f5c7c7ed5eaff296026a580f07aba76d2e444a8"
BASE_HASHES = {}
ROOT_KEYS = {}
PIPELINES = ("root4k3500std", "root4k3596std")

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
        raise ValueError("8000-point answer sidecar changed")
    points = [group.point_from_json(row) for row in json_lines(HOLDOUT)]
    answers = json_lines(ANSWERS)
    if len(points) != len(answers) or len(points) != 8000 or len(set(points)) != len(points):
        raise ValueError("8000-point target panel is incomplete or has duplicate points")
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
    if len(values) != 8000:
        raise ValueError("high-K confirmation requires 8000 paired targets")
    block_means = [statistics.mean(values[i:i + 25])
                   for i in range(0, len(values), 25)]
    rng = random.Random(seed)
    samples = [statistics.mean(block_means[rng.randrange(len(block_means))]
                               for _ in block_means)
               for _ in range(10_000)]
    samples.sort()
    return [samples[49], samples[9949]]


def verify_arm(output, arm, k, variant, profile, points, scalars):
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
    if state_bound + 1 >= 1 << 29:
        raise ValueError(f"{arm}: four-byte state ID exceeds 29 bits")
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
            or summary["root_bloom_enabled"] is not (profile != "off")
            or (summary["root_bloom_profile"], summary["root_bloom_hashes"],
                summary["root_bloom_bits_per_request"], summary["root_bloom_bytes"]) !=
               (("wide3", 3, 16, 1073741824) if profile == "wide3"
                else ("standard", 3, 8, 1073741824) if profile == "standard"
                else ("none", 0, 0, 0))
            or (k in BASE_HASHES and summary["base_hash"] != BASE_HASHES[k])
            or (k in ROOT_KEYS and summary["root_table_entries"] != ROOT_KEYS[k])
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
    if variant == "root22" and (
            summary["root_table_request_entries"] != (state_bound * 4 + 2) // 3
            or summary["root_table_allocated_bytes"] != 5 * (1 << 30)):
        raise ValueError(f"{arm}: load-75 root22 table allocation mismatch")
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
    return {"arm": arm, "pipeline": arm.rsplit("_", 1)[0],
            "k": k, "variant": variant, "bloom_profile": profile,
            "status": "PASS_INDEPENDENT_REPLAY",
            "actual_factor_base_points": len(base), "effective_columns": k,
            "base_hash": summary["base_hash"],
            "root_table_entries": summary["root_table_entries"],
            "root_table_allocated_bytes": summary["root_table_allocated_bytes"],
            "root_state_allocated_bytes": summary["root_state_allocated_bytes"],
            "peak_rss_bytes": usage["peak_rss_bytes"],
            "cpu_throttled_periods_delta": throttled,
            "oom_kills_delta": oom_kills,
            "resource_envelope_status": (
                "WITHIN_15_25_GIB" if usage["peak_rss_bytes"] <= runner.parent.RSS_CAP_BYTES
                else "ABOVE_15_25_GIB"),
            "setup_ms": summary["timing_ms"]["setup_total"],
            "mean_online_ms_exploratory": statistics.mean(times),
            "mean_online_target_block_bootstrap_99pct_ms":
                target_block_mean_interval(times, 2026101617 + k),
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
    for arm, k, variant, profile, _ in runner.ARMS:
        try:
            arms.append(verify_arm(output, arm, k, variant, profile, points, scalars))
        except Exception as error:
            arms.append({"arm": arm,
                         "status": "VERIFICATION_FAILED", "error": str(error)})
    for pipeline in PIPELINES:
        repetitions = [row for row in arms
                       if row["arm"] in (f"{pipeline}_a", f"{pipeline}_b")
                       and row["status"] == "PASS_INDEPENDENT_REPLAY"]
        if len(repetitions) == 2 and (
                repetitions[0]["base_hash"] != repetitions[1]["base_hash"]
                or repetitions[0]["root_table_entries"] != repetitions[1]["root_table_entries"]):
            for row in repetitions:
                row["status"] = "REPEATED_BASE_OR_INDEX_MISMATCH"
    eligible = [row for row in arms if row["status"] == "PASS_INDEPENDENT_REPLAY"
                and row["resource_envelope_status"] == "WITHIN_15_25_GIB"
                and row["cpu_throttled_periods_delta"] == 0
                and row["oom_kills_delta"] == 0]
    eligible_by_arm = {row["arm"]: row for row in eligible}
    complete = len(eligible) == len(runner.ARMS)
    balanced = None
    selected_pipeline = None
    if complete:
        online = {arm: [target["online_ms_exploratory"] for target in row["targets"]]
                  for arm, row in eligible_by_arm.items()}
        means = {pipeline: [(online[f"{pipeline}_a"][i] + online[f"{pipeline}_b"][i]) / 2
                            for i in range(len(points))] for pipeline in PIPELINES}
        pairwise = {}
        for i, candidate in enumerate(PIPELINES):
            for j, other in enumerate(PIPELINES):
                if candidate == other:
                    continue
                delta = [left - right for left, right in zip(means[candidate], means[other])]
                pairwise[f"{candidate}_minus_{other}"] = {
                    "mean_ms": statistics.mean(delta),
                    "contiguous_block_bootstrap_99pct_ms":
                        target_block_mean_interval(delta, 2026101651 + i * 2 + j),
                }
        fastest = min(PIPELINES, key=lambda pipeline: statistics.mean(means[pipeline]))
        if all(pairwise[f"{fastest}_minus_{other}"]["contiguous_block_bootstrap_99pct_ms"][1] < 0
               for other in PIPELINES if other != fastest):
            selected_pipeline = fastest
        repeatability = {}
        for index, pipeline in enumerate(PIPELINES):
            same_pipeline = [online[f"{pipeline}_a"][i] - online[f"{pipeline}_b"][i]
                             for i in range(len(points))]
            repeatability[f"{pipeline}_a_minus_b"] = {
                "mean_ms": statistics.mean(same_pipeline),
                "contiguous_block_bootstrap_99pct_ms":
                    target_block_mean_interval(same_pipeline, 2026101661 + index),
            }
        balanced = {
            "pipeline_mean_ms": {pipeline: statistics.mean(means[pipeline])
                                 for pipeline in PIPELINES},
            "contiguous_block_size": 25,
            "pairwise": pairwise,
            "repeatability": repeatability,
            "selected_pipeline": selected_pipeline,
        }
    source = one_json(output / "source.json")
    if (source["schema"] != "n83.root4-highk-confirm-source/v1"
            or source["hashes"] != runner.frozen_inputs()
            or source["arm_order"] != [arm for arm, _, _, _, _ in runner.ARMS]
            or source["panel_runner_sha256"] != runner.parent.digest(RUNNER)):
        raise ValueError("executed balanced-order source or fixture changed")
    analysis = {"schema": "n83.root4-highk-confirm-analysis/v1",
                "hardware_tier": "cloud_instance_exploratory",
                "workload_role": "secondary_fresh_8000_point_balanced_highk_root4_confirmation",
                "target_count": len(points), "arms": arms,
                "arm_order": source["arm_order"],
                "selection_rule": "root4k3500std/root4k3596std/root4k3596std/root4k3500std with fresh empty index and logs; 15.25 GiB observed RSS guard; average mirrored runs per pipeline on each of 8000 fresh public points; 320 contiguous 25-target blocks, 10000 resamples and fixed pair seeds 2026101651+i*2+j; select the lower-mean pipeline only when its upper 99% paired endpoint is below zero versus the other",
                "balanced_comparison": balanced,
                "selected_pipeline": selected_pipeline,
                "status": "COMPLETE" if complete else "INCOMPLETE",
                "bootstrap_resamples": 10_000,
                "controlled_speedup": None}
    (output / "analysis.json").write_text(json.dumps(analysis, indent=2,
                                                    sort_keys=True) + "\n")
    print(json.dumps(analysis, sort_keys=True))
    return 0 if complete else 1


if __name__ == "__main__":
    raise SystemExit(main())
