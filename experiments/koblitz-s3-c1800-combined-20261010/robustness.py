"""Reanalyze frozen N53 target and repeat variance without dropping runs."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import random
import statistics


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "analysis.json"
SEED = 20261010
BOOTSTRAP_DRAWS = 100_000


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def geometric(values: list[float]) -> float:
    return math.exp(statistics.mean(map(math.log, values)))


def percentile(ordered: list[float], fraction: float) -> float:
    position = (len(ordered) - 1) * fraction
    low = int(position)
    return ordered[low] + (ordered[min(low + 1, len(ordered) - 1)] -
                           ordered[low]) * (position - low)


def analyze(output: Path) -> dict:
    if output.exists():
        raise FileExistsError(output)
    source_raw = SOURCE.read_bytes()
    source = json.loads(source_raw)
    assert source["run_count"] == source["verified_run_count"] == 32
    assert source["complete_point_count"] == 8
    assert source["independent_sage_replay_verified"]
    points = source["points"]
    assert len(points) == 8
    assert all(len(point[f"{arm}_online_ms"]) == 2
               for point in points for arm in ("baseline", "combined"))

    point_logs = [math.log(point["online_ratio"]) for point in points]
    assert abs(math.exp(statistics.mean(point_logs)) -
               source["online_geomean_baseline_over_combined"]) < 1e-12
    leave_one_out = [
        {"omitted_label": point["label"],
         "remaining_target_geometric_ratio":
             math.exp(statistics.mean(point_logs[:index] + point_logs[index + 1:]))}
        for index, point in enumerate(points)
    ]

    # Resample target identities and the two fresh process timings in each
    # arm. All original runs stay eligible; this interval adds run noise to
    # the target-only bootstrap in analysis.json.
    rng = random.Random(SEED)
    bootstrap = []
    for _ in range(BOOTSTRAP_DRAWS):
        draw_logs = []
        for _ in points:
            point = points[rng.randrange(len(points))]
            baseline = point["baseline_online_ms"]
            combined = point["combined_online_ms"]
            a = math.sqrt(baseline[rng.randrange(2)] * baseline[rng.randrange(2)])
            b = math.sqrt(combined[rng.randrange(2)] * combined[rng.randrange(2)])
            draw_logs.append(math.log(a / b))
        bootstrap.append(math.exp(statistics.mean(draw_logs)))
    bootstrap.sort()

    # Exact variance of the empirical target/repeat resampling law. For a
    # target, each arm's two log timings are sampled twice with replacement.
    # A two-point mean contributes squared log-repeat difference / 8.
    conditional_means = []
    conditional_variances = []
    for point in points:
        a = list(map(math.log, point["baseline_online_ms"]))
        b = list(map(math.log, point["combined_online_ms"]))
        conditional_means.append(statistics.mean(a) - statistics.mean(b))
        conditional_variances.append(((a[0] - a[1]) ** 2 +
                                      (b[0] - b[1]) ** 2) / 8)
    log_mean = statistics.mean(conditional_means)
    between_target_variance = statistics.pvariance(conditional_means)
    within_target_repeat_variance = statistics.mean(conditional_variances)
    log_sd = math.sqrt(between_target_variance + within_target_repeat_variance)
    distance_to_twofold = log_mean - math.log(2)
    assert distance_to_twofold > 0
    expected_n = math.ceil((1.96 * log_sd / distance_to_twofold) ** 2)
    power80_n = math.ceil(((1.96 + .84) * log_sd / distance_to_twofold) ** 2)

    result = {
        "kind": "s3_n53_c1800_combined_robustness_v1",
        "source_analysis_sha256": sha(source_raw),
        "point_count": len(points),
        "run_count": source["run_count"],
        "all_native_and_sage_verified": True,
        "target_only_bootstrap95": source["online_bootstrap95_target_resample"],
        "hierarchical_bootstrap": {
            "seed": SEED, "draws": BOOTSTRAP_DRAWS,
            "unit": "resample target, then two process timings per arm",
            "geometric_ratio_95": [percentile(bootstrap, .025),
                                   percentile(bootstrap, .975)],
            "resampled_fraction_ratio_at_least_2":
                sum(value >= 2 for value in bootstrap) / len(bootstrap),
        },
        "leave_one_target_out": leave_one_out,
        "leave_one_target_out_ratio_range": [
            min(row["remaining_target_geometric_ratio"] for row in leave_one_out),
            max(row["remaining_target_geometric_ratio"] for row in leave_one_out)],
        "targets_with_ratio_at_least_2": sum(point["online_ratio"] >= 2 for point in points),
        "sample_size_planning": {
            "model": "normal_approximation_to_empirical_target_plus_repeat_log_ratio_law",
            "log_ratio_mean": log_mean,
            "between_target_log_variance": between_target_variance,
            "within_target_repeat_log_variance": within_target_repeat_variance,
            "combined_log_standard_deviation": log_sd,
            "threshold_log_ratio": math.log(2),
            "estimated_distinct_targets_for_expected_two_sided_95_lower_above_2": expected_n,
            "estimated_distinct_targets_for_80_percent_power": power80_n,
            "process_runs_per_target": 4,
            "caveat": "local host variance may not transfer to the confirmation host",
        },
        "host_exclusive_partition_receipt": None,
        "controlled_speedup": None,
    }
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "hierarchical_bootstrap95": result["hierarchical_bootstrap"]["geometric_ratio_95"],
        "planning_distinct_targets_80_percent_power": power80_n,
        "leave_one_out_range": result["leave_one_target_out_ratio_range"],
    }, sort_keys=True))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    analyze(args.output)
