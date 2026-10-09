"""Analyze all frozen one-target pairs, including failures and fallbacks."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import random
import statistics


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def bootstrap_interval(log_ratios: list[float]) -> list[float] | None:
    if len(log_ratios) < 2:
        return None
    randomizer = random.Random(20261009)
    estimates = sorted(
        math.exp(sum(randomizer.choice(log_ratios) for _ in log_ratios) / len(log_ratios))
        for _ in range(20000)
    )
    return [estimates[499], estimates[19499]]


def wilson(successes: int, count: int) -> list[float]:
    z = 1.959963984540054
    p = successes / count
    denominator = 1 + z * z / count
    center = (p + z * z / (2 * count)) / denominator
    half = z * math.sqrt(p * (1 - p) / count + z * z / (4 * count * count)) / denominator
    return [max(0, center - half), min(1, center + half)]


def main(ledger: Path, archive: Path, output: Path,
         first_target: int, target_count_per_curve: int, plan: Path) -> None:
    if output.exists():
        raise FileExistsError(output)
    rows = [json.loads(line) for line in ledger.read_text().splitlines()]
    freeze = json.loads(plan.read_text())
    assert len(rows) == 4 * target_count_per_curve
    expected_labels = {f"T{index:03}" for index in range(
        first_target, first_target + target_count_per_curve)}
    assert {(row["n"], row["label"], row["arm"]) for row in rows} == {
        (n, label, arm) for n in (41, 53)
        for label in expected_labels for arm in ("reference", "partial")}
    curves = []
    for n in (41, 53):
        by_target: dict[str, dict] = {}
        for row in rows:
            if row["n"] == n:
                by_target.setdefault(row["label"], {})[row["arm"]] = row
        assert set(by_target) == expected_labels
        paired = []
        for label in sorted(by_target):
            arms = by_target[label]
            assert set(arms) == {"reference", "partial"}
            reference, partial = arms["reference"], arms["partial"]
            pair_verified = bool(reference["verified_native_and_fixture"]
                                 and partial["verified_native_and_fixture"])
            ref = reference["report"] if pair_verified else None
            part = partial["report"] if pair_verified else None
            if pair_verified:
                assert reference["target_sha256"] == partial["target_sha256"]
                assert ref["target"] == part["target"]
                assert ref["recovered_scalar"] == part["recovered_scalar"]
                assert ref["factor_base_digest"] == part["factor_base_digest"]
                assert part["prepared_rank_reached"] is True
                assert part["prepare_rank_requested"] == freeze["rank_thresholds"][str(n)]
                assert part["target_s3_calls"] == ref["target_s3_calls"]
                ref_online = ref["timing_ms"]["target_online_after_reusable_setup"]
                part_online = part["timing_ms"]["target_online_after_reusable_setup"]
                ref_calls = ref["rank_target_s3_calls"] + ref["target_s3_calls"]
                part_online_calls = part["fallback_s3_calls"] + part["target_s3_calls"]
                part_cold_calls = part["rank_target_s3_calls"] + part_online_calls
                ref_trace = ref["rank_relation_witnesses"]
                part_trace = part["rank_relation_witnesses"] + part["fallback_relation_witnesses"]
                prefix_equal = ref_trace[:min(len(ref_trace), len(part_trace))] == part_trace[:min(len(ref_trace), len(part_trace))]
                ref_rank_yield = ref["rank_new_rows"] / ref["rank_attempts"]
                part_prep_rank_yield = part["rank_new_rows"] / part["rank_attempts"]
                part_fallback_rank_yield = (part["fallback_new_rows"] / part["fallback_attempts"]
                                            if part["fallback_attempts"] else None)
                ref_cost_per_new_row = ref["timing_ms"]["rank_collection_total"] / ref["rank_new_rows"]
                part_prep_cost_per_new_row = part["timing_ms"]["rank_collection_total"] / part["rank_new_rows"]
                part_fallback_cost_per_new_row = (
                    part["timing_ms"]["fallback_collection_total"] / part["fallback_new_rows"]
                    if part["fallback_new_rows"] else None)
            else:
                ref_online = part_online = ref_calls = part_online_calls = part_cold_calls = None
                prefix_equal = None
                ref_rank_yield = part_prep_rank_yield = part_fallback_rank_yield = None
                ref_cost_per_new_row = part_prep_cost_per_new_row = part_fallback_cost_per_new_row = None
            paired.append({"label": label, "target_sha256": reference["target_sha256"],
                           "pair_verified": pair_verified,
                           "reference_raw_sha256": reference["raw_sha256"],
                           "partial_raw_sha256": partial["raw_sha256"],
                           "reference_online_ms": ref_online,
                           "partial_online_ms": part_online,
                           "exploratory_online_ratio": ref_online / part_online if pair_verified else None,
                           "reference_cold_ms": ref["timing_ms"]["cold_after_launch"] if pair_verified else None,
                           "partial_cold_ms": part["timing_ms"]["cold_after_launch"] if pair_verified else None,
                           "reference_online_s3_root_calls": ref_calls,
                           "partial_online_s3_root_calls": part_online_calls,
                           "partial_cold_s3_root_calls": part_cold_calls,
                           "reference_cold_s3_root_calls": ref_calls,
                           "relation_trace_common_prefix_equal": prefix_equal,
                           "partial_prepared_rank": part["rank"] if pair_verified else None,
                           "partial_fallback_attempts": part["fallback_attempts"] if pair_verified else None,
                           "partial_final_rank": part["final_rank"] if pair_verified else None,
                           "reference_rank_yield_per_attempt": ref_rank_yield,
                           "partial_preparation_rank_yield_per_attempt": part_prep_rank_yield,
                           "partial_fallback_rank_yield_per_attempt": part_fallback_rank_yield,
                           "reference_rank_collection_ms_per_new_row": ref_cost_per_new_row,
                           "partial_preparation_ms_per_new_row": part_prep_cost_per_new_row,
                           "partial_fallback_ms_per_new_row": part_fallback_cost_per_new_row,
                           "reference_factor_base_ms": ref["timing_ms"]["factor_base_construct"] if pair_verified else None,
                           "partial_factor_base_ms": part["timing_ms"]["factor_base_construct"] if pair_verified else None,
                           "reference_index_build_ms": ref["timing_ms"]["index_build"] if pair_verified else None,
                           "partial_index_build_ms": part["timing_ms"]["index_build"] if pair_verified else None,
                           "reference_peak_rss_bytes": ref["peak_rss_bytes"] if pair_verified else None,
                           "partial_peak_rss_bytes": part["peak_rss_bytes"] if pair_verified else None})
        all_verified = all(pair["pair_verified"] for pair in paired)
        complete = [pair for pair in paired if pair["pair_verified"]]
        ratios = [pair["exploratory_online_ratio"] for pair in complete]
        fallback_free = sum(pair["partial_fallback_attempts"] == 0 for pair in complete)
        at_least_2 = sum(ratio >= 2 for ratio in ratios)
        curves.append({"n": n, "target_count": target_count_per_curve,
                       "run_count": 2 * target_count_per_curve,
                       "all_verified": all_verified,
                       "fallback_free_count": fallback_free,
                       "fallback_count": len(complete) - fallback_free,
                       "fallback_free_wilson_95": wilson(fallback_free, len(complete)) if complete else None,
                       "paired_ratio_at_least_2_count": at_least_2,
                       "paired_ratio_at_least_2_wilson_95": wilson(at_least_2, len(complete)) if complete else None,
                       "paired_ratio_below_2_labels": [pair["label"] for pair in complete
                                                        if pair["exploratory_online_ratio"] < 2],
                       "all_common_relation_prefixes_equal": all(pair["relation_trace_common_prefix_equal"] for pair in complete),
                       "exploratory_geometric_mean_online_ratio":
                           math.exp(statistics.mean(map(math.log, ratios))) if all_verified else None,
                       "exploratory_target_bootstrap_95":
                           bootstrap_interval(list(map(math.log, ratios))) if all_verified else None,
                       "exploratory_min_online_ratio": min(ratios) if all_verified else None,
                       "all_pairs_at_least_2": all(ratio >= 2 for ratio in ratios) if all_verified else None,
                       "median_reference_online_ms": statistics.median(pair["reference_online_ms"] for pair in complete) if all_verified else None,
                       "median_partial_online_ms": statistics.median(pair["partial_online_ms"] for pair in complete) if all_verified else None,
                       "median_reference_cold_ms": statistics.median(pair["reference_cold_ms"] for pair in complete) if all_verified else None,
                       "median_partial_cold_ms": statistics.median(pair["partial_cold_ms"] for pair in complete) if all_verified else None,
                       "median_reference_rank_yield_per_attempt": statistics.median(pair["reference_rank_yield_per_attempt"] for pair in complete) if all_verified else None,
                       "median_partial_preparation_rank_yield_per_attempt": statistics.median(pair["partial_preparation_rank_yield_per_attempt"] for pair in complete) if all_verified else None,
                       "median_reference_rank_collection_ms_per_new_row": statistics.median(pair["reference_rank_collection_ms_per_new_row"] for pair in complete) if all_verified else None,
                       "median_partial_preparation_ms_per_new_row": statistics.median(pair["partial_preparation_ms_per_new_row"] for pair in complete) if all_verified else None,
                       "median_reference_factor_base_ms": statistics.median(pair["reference_factor_base_ms"] for pair in complete) if all_verified else None,
                       "median_partial_factor_base_ms": statistics.median(pair["partial_factor_base_ms"] for pair in complete) if all_verified else None,
                       "median_reference_index_build_ms": statistics.median(pair["reference_index_build_ms"] for pair in complete) if all_verified else None,
                       "median_partial_index_build_ms": statistics.median(pair["partial_index_build_ms"] for pair in complete) if all_verified else None,
                       "min_online_s3_root_call_ratio":
                           min(pair["reference_online_s3_root_calls"] / pair["partial_online_s3_root_calls"]
                               for pair in complete) if all_verified else None,
                       "cold_s3_root_calls_equal_on_all_pairs":
                           all(pair["reference_cold_s3_root_calls"] == pair["partial_cold_s3_root_calls"]
                               for pair in complete) if all_verified else None,
                       "pairs": paired})
    result = {"kind": "s3_partial_prepared_fresh_holdout_analysis_v1",
              "scope": (f"one target per process, {target_count_per_curve} points on each curve, "
                        "reference versus fixed partial rank with target-charged fallback"),
              "first_target_index": first_target,
              "target_count_per_curve": target_count_per_curve,
              "frozen_plan_sha256": sha(plan),
              "ledger_sha256": sha(ledger), "raw_archive_sha256": sha(archive),
              "run_count": len(rows), "failure_count": sum(not row["verified_native_and_fixture"] for row in rows),
              "strict_host_isolation_receipt": None, "controlled_speedup": None,
              "a_a_noise_control": None, "curves": curves}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"output": str(output),
                      "curves": [{"n": curve["n"], "all_verified": curve["all_verified"],
                                  "fallback_count": curve["fallback_count"],
                                  "geometric_mean_online_ratio": curve["exploratory_geometric_mean_online_ratio"],
                                  "minimum_online_ratio": curve["exploratory_min_online_ratio"]}
                                 for curve in curves]}, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", required=True, type=Path)
    parser.add_argument("--archive", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--first-target", type=int, default=1)
    parser.add_argument("--target-count-per-curve", type=int, default=24)
    parser.add_argument("--plan", type=Path,
                        default=Path(__file__).resolve().parent / "FREEZE.json")
    args = parser.parse_args()
    main(args.ledger, args.archive, args.output,
         args.first_target, args.target_count_per_curve, args.plan)
