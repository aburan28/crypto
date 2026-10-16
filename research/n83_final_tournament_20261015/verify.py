#!/usr/bin/env python3
"""Replay finalist bases and 4,000 public targets, then apply the frozen rule."""

import argparse
import importlib.util
import json
import random
import statistics
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE_VERIFIER = HERE.parent / "n83_root4_bloom_20261015/verify.py"
RUNNER = HERE / "run.py"
TARGET_COUNT = 4000
BLOCK_SIZE = 50
RESAMPLES = 10_000

rspec = importlib.util.spec_from_file_location("n83_final_runner", RUNNER)
runner = importlib.util.module_from_spec(rspec)
rspec.loader.exec_module(runner)
vspec = importlib.util.spec_from_file_location("n83_base_verifier", BASE_VERIFIER)
base = importlib.util.module_from_spec(vspec)
vspec.loader.exec_module(base)
base.runner = runner


def block_interval(values, seed):
    if len(values) != TARGET_COUNT:
        raise ValueError("final comparison requires 4000 paired public points")
    blocks = [statistics.mean(values[i:i + BLOCK_SIZE])
              for i in range(0, TARGET_COUNT, BLOCK_SIZE)]
    rng = random.Random(seed)
    samples = [statistics.mean(blocks[rng.randrange(len(blocks))]
                               for _ in blocks)
               for _ in range(RESAMPLES)]
    samples.sort()
    return [samples[49], samples[9949]]


base.target_block_mean_interval = block_interval


def targets():
    admission, _ = runner.frozen()
    points = [base.group.point_from_json(row)
              for row in base.json_lines(runner.TARGETS)]
    answers = base.json_lines(runner.ANSWERS)
    if (len(points) != TARGET_COUNT or len(answers) != TARGET_COUNT
            or len(set(points)) != TARGET_COUNT
            or not base.group.on_curve(base.group.G)
            or base.group.multiply(base.group.G, base.group.R) is not None):
        raise ValueError("final target panel or subgroup is invalid")
    scalars = []
    for index, (point, answer) in enumerate(zip(points, answers)):
        scalar = int(answer["scalar"])
        if (answer["fixture_index"] != index
                or base.group.point_from_json(answer["point"]) != point
                or not base.group.on_curve(point)
                or base.group.multiply(base.group.G, scalar) != point):
            raise ValueError(f"final public target {index} fails scalar replay")
        scalars.append(scalar)
    return admission, points, scalars


def verify_arm(output, arm, candidate, points, scalars):
    row = base.verify_arm(output, arm, candidate["k"], candidate["variant"],
                          candidate["bloom_profile"], points, scalars)
    if row["status"] == "PASS_INDEPENDENT_REPLAY":
        if (row["base_hash"] != candidate["base_hash"]
                or row["root_table_entries"] != candidate["root_table_entries"]
                or row["actual_factor_base_points"] != candidate["actual_usable_points"]):
            raise ValueError(f"{arm}: finalist base or root index changed")
        row["candidate_id"] = candidate["candidate_id"]
    before = output / f"{arm}_system_before.json"
    after = output / f"{arm}_system_after.json"
    if before.exists() and after.exists():
        first = base.one_json(before)["proc_stat"]
        last = base.one_json(after)["proc_stat"]
        def steal_ticks(proc_stat):
            fields = proc_stat.splitlines()[0].split()
            if fields[0] != "cpu" or len(fields) < 9:
                raise ValueError("system steal counter unavailable")
            return int(fields[8])
        row["system_steal_ticks_delta"] = steal_ticks(last) - steal_ticks(first)
    else:
        row["system_steal_ticks_delta"] = None
        if row["status"] == "PASS_INDEPENDENT_REPLAY":
            row["status"] = "SYSTEM_COUNTER_MISSING"
    return row


def calibration_receipt(output):
    before = output / "calibration_before.json"
    after = output / "calibration_after.json"
    if not (before.exists() and after.exists()):
        return None
    first, last = base.one_json(before), base.one_json(after)
    if ([item["sha256"] for item in first["checksums_and_cpu"]] !=
            [item["sha256"] for item in last["checksums_and_cpu"]]):
        raise ValueError("calibration checksums changed")
    return {"before": first, "after": last,
            "median_ratio_after_over_before":
                last["median_ms"] / first["median_ms"]}


def comparison(arms, candidates):
    eligible = {row["arm"]: row for row in arms
                if (row["status"] == "PASS_INDEPENDENT_REPLAY"
                    and row["resource_envelope_status"] == "WITHIN_14_5_GIB"
                    and row["cpu_throttled_periods_delta"] == 0
                    and row["oom_kills_delta"] == 0)}
    expected = {arm for arm, _ in runner.arms(candidates)}
    if len(arms) != len(expected) or set(eligible) != expected:
        return None
    averages = {}
    repeatability = {}
    for i, candidate in enumerate(candidates):
        pipeline = candidate["pipeline"]
        first = eligible[f"{pipeline}_a"]["targets"]
        second = eligible[f"{pipeline}_b"]["targets"]
        a = [row["online_ms_exploratory"] for row in first]
        b = [row["online_ms_exploratory"] for row in second]
        averages[pipeline] = [(left + right) / 2 for left, right in zip(a, b)]
        repeatability[f"{pipeline}_a_minus_b"] = {
            "mean_ms": statistics.mean(x - y for x, y in zip(a, b)),
            "contiguous_block_bootstrap_99pct_ms":
                block_interval([x - y for x, y in zip(a, b)], 20261015950 + i),
        }
    pairwise = {}
    for i, left in enumerate(candidates):
        for j, right in enumerate(candidates):
            if i == j:
                continue
            a, b = averages[left["pipeline"]], averages[right["pipeline"]]
            delta = [x - y for x, y in zip(a, b)]
            pairwise[f"{left['pipeline']}_minus_{right['pipeline']}"] = {
                "mean_ms": statistics.mean(delta),
                "contiguous_block_bootstrap_99pct_ms":
                    block_interval(delta, 20261015960 + i * 16 + j),
            }
    fastest = min(candidates, key=lambda row: (
        statistics.mean(averages[row["pipeline"]]), row["pipeline"]))
    selected = fastest if all(
        pairwise[f"{fastest['pipeline']}_minus_{other['pipeline']}"][
            "contiguous_block_bootstrap_99pct_ms"][1] < 0
        for other in candidates if other != fastest) else None
    return {
        "pipeline_mean_ms": {key: statistics.mean(value)
                             for key, value in averages.items()},
        "contiguous_block_size": BLOCK_SIZE,
        "bootstrap_resamples": RESAMPLES,
        "pairwise": pairwise,
        "repeatability": repeatability,
        "descriptive_lowest_pipeline": fastest["pipeline"],
        "selected_pipeline": selected["pipeline"] if selected else None,
        "selected_candidate_id": selected["candidate_id"] if selected else None,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    output = args.output.resolve()
    admission, points, scalars = targets()
    candidates = admission["candidate_order"]
    arms = []
    for arm, candidate in runner.arms(candidates):
        try:
            arms.append(verify_arm(output, arm, candidate, points, scalars))
        except Exception as error:
            arms.append({"arm": arm, "candidate_id": candidate["candidate_id"],
                         "status": "VERIFICATION_FAILED", "error": str(error)})
    source = base.one_json(output / "source.json")
    _, hashes = runner.frozen()
    if (source["schema"] != "n83.final-tournament-source/v1"
            or source["hashes"] != hashes
            or source["admission_sha256"] != hashes[str(runner.ADMISSION.relative_to(runner.ROOT))]
            or source["arm_order"] != [arm for arm, _ in runner.arms(candidates)]
            or source["candidate_ids"] != [row["candidate_id"] for row in candidates]
            or source["panel_runner_sha256"] != runner.digest(RUNNER)):
        raise ValueError("executed source, admission or public fixture differs")
    calibration = calibration_receipt(output)
    balanced = comparison(arms, candidates) if calibration is not None else None
    selected = balanced["selected_pipeline"] if balanced else None
    analysis = {
        "schema": "n83.final-tournament-analysis/v1",
        "hardware_tier": "cloud_instance_exploratory",
        "workload_role": "secondary_fresh_4000_point_final_factor_base_selection",
        "target_count": TARGET_COUNT,
        "candidate_order": [row["pipeline"] for row in candidates],
        "arm_order": source["arm_order"],
        "arms": arms,
        "calibration": calibration,
        "balanced_comparison": balanced,
        "selected_pipeline": selected,
        "status": "COMPLETE" if balanced else "INCOMPLETE",
        "decision_status": "SELECTED" if selected else "BOUNDED_UNKNOWN",
        "controlled_speedup": None,
    }
    (output / "analysis.json").write_text(json.dumps(analysis, indent=2,
                                                    sort_keys=True) + "\n")
    print(json.dumps({"status": analysis["status"],
                      "decision_status": analysis["decision_status"],
                      "selected_pipeline": selected}, sort_keys=True))
    return 0 if balanced else 1


if __name__ == "__main__":
    raise SystemExit(main())
