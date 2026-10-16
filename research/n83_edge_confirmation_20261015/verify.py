#!/usr/bin/env python3
"""Independently replay the fresh confirmation and apply its frozen rule."""

import argparse
import importlib.util
import json
import random
import statistics
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUNNER = HERE / "run.py"
SOURCE = HERE.parent / "n83_final_tournament_20261015/verify.py"
TARGET_COUNT = 4000
BLOCK_SIZE = 50
RESAMPLES = 10_000

rspec = importlib.util.spec_from_file_location("n83_confirmation_runner", RUNNER)
runner = importlib.util.module_from_spec(rspec)
rspec.loader.exec_module(runner)
vspec = importlib.util.spec_from_file_location("n83_confirmation_verifier_source", SOURCE)
source = importlib.util.module_from_spec(vspec)
vspec.loader.exec_module(source)
source.runner = runner
source.base.runner = runner


def interval(values, seed):
    if len(values) != TARGET_COUNT:
        raise ValueError("confirmation requires 4000 paired points")
    blocks = [statistics.mean(values[i:i + BLOCK_SIZE])
              for i in range(0, TARGET_COUNT, BLOCK_SIZE)]
    rng = random.Random(seed)
    samples = [statistics.mean(blocks[rng.randrange(len(blocks))]
                               for _ in blocks)
               for _ in range(RESAMPLES)]
    samples.sort()
    return [samples[49], samples[9949]]


def calibration_status(calibration):
    if calibration is None:
        return "MISSING"
    before, after = calibration["before"], calibration["after"]
    for receipt in (before, after):
        samples = receipt["sample_ms"]
        if (len(samples) != 3 or min(samples) <= 0
                or max(samples) / min(samples) > 1.15):
            return "NOISY"
    ratio = calibration["median_ratio_after_over_before"]
    return "PASS" if 1 / 1.15 <= ratio <= 1.15 else "NOISY"


def comparison(arms, candidates):
    eligible = {row["arm"]: row for row in arms
                if (row["status"] == "PASS_INDEPENDENT_REPLAY"
                    and row["resource_envelope_status"] == "WITHIN_14_5_GIB"
                    and row["cpu_throttled_periods_delta"] == 0
                    and row["oom_kills_delta"] == 0)}
    if len(arms) != 2 * len(candidates) or set(eligible) != {
            arm for arm, _ in runner.arms(candidates)}:
        return None
    averages = {}
    repeatability = {}
    for i, candidate in enumerate(candidates):
        name = candidate["pipeline"]
        first = [row["online_ms_exploratory"] for row in eligible[f"{name}_a"]["targets"]]
        second = [row["online_ms_exploratory"] for row in eligible[f"{name}_b"]["targets"]]
        averages[name] = [(a + b) / 2 for a, b in zip(first, second)]
        repeatability[f"{name}_a_minus_b"] = {
            "mean_ms": statistics.mean(a - b for a, b in zip(first, second)),
            "contiguous_block_bootstrap_99pct_ms":
                interval([a - b for a, b in zip(first, second)], 20261015600 + i),
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
                    interval(delta, 20261015601 + 16 * i + j),
            }
    fastest = min(candidates, key=lambda row: (
        statistics.mean(averages[row["pipeline"]]), row["pipeline"]))
    selected = fastest if all(
        pairwise[f"{fastest['pipeline']}_minus_{other['pipeline']}"][
            "contiguous_block_bootstrap_99pct_ms"][1] < 0
        for other in candidates if other != fastest) else None
    return {
        "pipeline_mean_ms": {name: statistics.mean(values)
                             for name, values in averages.items()},
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
    admission, points, scalars = source.targets()
    candidates = admission["candidate_order"]
    arms = []
    for arm, candidate in runner.arms(candidates):
        try:
            arms.append(source.verify_arm(output, arm, candidate, points, scalars))
        except Exception as error:
            arms.append({"arm": arm, "candidate_id": candidate["candidate_id"],
                         "status": "VERIFICATION_FAILED", "error": str(error)})
    receipt = source.base.one_json(output / "source.json")
    _, hashes = runner.frozen()
    if (receipt["schema"] != "n83.edge-confirmation-source/v1"
            or receipt["hashes"] != hashes
            or receipt["admission_sha256"] != runner.digest(runner.ADMISSION)
            or receipt["arm_order"] != [arm for arm, _ in runner.arms(candidates)]
            or receipt["candidate_ids"] != [row["candidate_id"] for row in candidates]
            or receipt["panel_runner_sha256"] != runner.digest(RUNNER)):
        raise ValueError("executed source or public fixture differs")
    calibration = source.calibration_receipt(output)
    noise = calibration_status(calibration)
    balanced = comparison(arms, candidates) if noise == "PASS" else None
    selected = balanced["selected_pipeline"] if balanced else None
    analysis = {
        "schema": "n83.edge-confirmation-analysis/v1",
        "hardware_tier": "cloud_instance_exploratory",
        "workload_role": "secondary_fresh_4000_point_factor_base_confirmation",
        "target_count": TARGET_COUNT,
        "candidate_order": [row["pipeline"] for row in candidates],
        "arm_order": receipt["arm_order"],
        "arms": arms,
        "calibration": calibration,
        "calibration_status": noise,
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
