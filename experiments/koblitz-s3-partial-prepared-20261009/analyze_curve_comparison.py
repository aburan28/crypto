"""Summarize 100 one-target pairs per curve and size a controlled replication."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import random
import statistics

from analyze_holdout import sha, wilson


def main(first: Path, extension: Path, output: Path) -> None:
    if output.exists():
        raise FileExistsError(output)
    inputs = [json.loads(path.read_text()) for path in (first, extension)]
    assert [item["target_count_per_curve"] for item in inputs] == [24, 76]
    assert [item["first_target_index"] for item in inputs] == [1, 25]
    assert all(item["failure_count"] == 0 for item in inputs)
    by_curve = {}
    for n in (41, 53):
        pairs = [pair for item in inputs for curve in item["curves"]
                 if curve["n"] == n for pair in curve["pairs"]]
        assert len(pairs) == 100 and len({pair["label"] for pair in pairs}) == 100
        assert all(pair["pair_verified"] for pair in pairs)
        ratios = [pair["exploratory_online_ratio"] for pair in pairs]
        logs = list(map(math.log, ratios))
        fallback_free = sum(pair["partial_fallback_attempts"] == 0 for pair in pairs)
        at_least_2 = sum(ratio >= 2 for ratio in ratios)
        by_curve[n] = {"target_count": 100,
                       "fallback_count": 100 - fallback_free,
                       "fallback_free_wilson_95": wilson(fallback_free, 100),
                       "at_least_2_count": at_least_2,
                       "at_least_2_wilson_95": wilson(at_least_2, 100),
                       "below_2_labels": [pair["label"] for pair in pairs
                                          if pair["exploratory_online_ratio"] < 2],
                       "exploratory_geometric_mean_online_ratio":
                           math.exp(statistics.mean(logs)),
                       "minimum_exploratory_online_ratio": min(ratios),
                       "log_ratio_standard_deviation": statistics.stdev(logs),
                       "all_cold_s3_root_calls_equal_reference": all(
                           pair["reference_cold_s3_root_calls"] ==
                           pair["partial_cold_s3_root_calls"] for pair in pairs),
                       "all_relation_trace_prefixes_equal": all(
                           pair["relation_trace_common_prefix_equal"] for pair in pairs)}
        by_curve[n]["log_ratios"] = logs
    a, b = by_curve[41]["log_ratios"], by_curve[53]["log_ratios"]
    difference = statistics.mean(a) - statistics.mean(b)
    randomizer = random.Random(20261009)
    bootstraps = sorted(math.exp(
        statistics.mean(randomizer.choices(a, k=100)) -
        statistics.mean(randomizer.choices(b, k=100))) for _ in range(20000))
    # Exploratory normal-approximation calculation for equal target counts.
    # Re-estimate variance with repeated paired runs on an isolated host.
    n_80_percent = math.ceil(
        (1.959963984540054 + 0.8416212335729143)**2 *
        (statistics.variance(a) + statistics.variance(b)) / difference**2)
    for curve in by_curve.values():
        del curve["log_ratios"]
    result = {"kind": "s3_partial_rank237_100_target_curve_comparison_v1",
              "first_analysis_sha256": sha(first),
              "extension_analysis_sha256": sha(extension),
              "strict_host_isolation_receipt": None,
              "controlled_speedup": None,
              "curves": {str(n): value for n, value in by_curve.items()},
              "exploratory_ratio_of_curve_geometric_means": math.exp(difference),
              "exploratory_target_bootstrap_95": [bootstraps[499], bootstraps[19499]],
              "observed_log_mean_difference": difference,
              "posthoc_equal_targets_per_curve_for_80_percent_power_at_alpha_0_05": n_80_percent,
              "sample_size_assumptions": (
                  "two_sided_normal_approximation_on_independent_target_log_ratios; "
                  "effect_and_variance_from_contended_mac_panel; recalibrate_with_isolated_"
                  "host_AA_control_and_repeated_one_target_pairs")}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"output": str(output), "ratio_of_curve_geomeans":
                      result["exploratory_ratio_of_curve_geometric_means"],
                      "bootstrap_95": result["exploratory_target_bootstrap_95"],
                      "posthoc_n_per_curve": n_80_percent}, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--first", required=True, type=Path)
    parser.add_argument("--extension", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    main(args.first, args.extension, args.output)
