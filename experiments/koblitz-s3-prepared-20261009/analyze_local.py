"""Summarize the frozen paired local panel while retaining every raw run."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import statistics


PHASES = ("target_query", "target_pdp_charged", "target_relation_check",
          "target_descent", "target_recovery_check")


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(ledger: Path, archive: Path, output: Path) -> None:
    if output.exists():
        raise FileExistsError(output)
    rows = [json.loads(line) for line in ledger.read_text().splitlines()]
    assert len(rows) == 48
    curves = []
    for n in (41, 53):
        points: dict[str, dict] = {}
        for row in rows:
            if row["n"] == n:
                points.setdefault(row["label"], {})[row["arm"]] = row
        assert len(points) == 12
        paired = []
        for label, arms in points.items():
            assert set(arms) == {"candidate", "prepared"}
            baseline, prepared = arms["candidate"], arms["prepared"]
            for row in (baseline, prepared):
                report = row["report"]
                timing = report["timing_ms"]
                assert row["verified_native_and_fixture"] is True
                assert row["error"] is None and row["exit_code"] == 0
                assert abs(sum(timing[key] for key in PHASES) -
                           timing["target_online_after_reusable_setup"]) < 0.02
                assert report["group_verified"] is True
                assert report["target_relation_found"] is True
            assert baseline["target_sha256"] == prepared["target_sha256"]
            assert baseline["report"]["target"] == prepared["report"]["target"]
            assert baseline["report"]["recovered_scalar"] == prepared["report"]["recovered_scalar"]
            assert baseline["report"]["factor_base_digest"] == prepared["report"]["factor_base_digest"]
            assert baseline["report"]["target_s3_calls"] == prepared["report"]["target_s3_calls"]
            assert prepared["report"]["rank_complete"] is True
            base_time = baseline["report"]["timing_ms"]["target_online_after_reusable_setup"]
            prep_time = prepared["report"]["timing_ms"]["target_online_after_reusable_setup"]
            base_s3 = (baseline["report"]["rank_target_s3_calls"] +
                       baseline["report"]["target_s3_calls"])
            prep_s3 = prepared["report"]["target_s3_calls"]
            paired.append({"target_label": label, "target_sha256": baseline["target_sha256"],
                           "baseline_online_ms": base_time, "prepared_online_ms": prep_time,
                           "exploratory_online_ratio": base_time / prep_time,
                           "baseline_online_s3_root_calls": base_s3,
                           "prepared_online_s3_root_calls": prep_s3,
                           "online_s3_root_call_ratio": base_s3 / prep_s3,
                           "prepared_relation_setup_s3_root_calls": prepared["report"]["rank_target_s3_calls"],
                           "baseline_reusable_setup_ms": baseline["report"]["timing_ms"]["reusable_setup_total"],
                           "prepared_reusable_setup_ms": prepared["report"]["timing_ms"]["reusable_setup_total"],
                           "baseline_cold_ms": baseline["report"]["timing_ms"]["cold_after_launch"],
                           "prepared_cold_ms": prepared["report"]["timing_ms"]["cold_after_launch"],
                           "baseline_raw_sha256": baseline["raw_sha256"],
                           "prepared_raw_sha256": prepared["raw_sha256"]})
        ratios = [pair["exploratory_online_ratio"] for pair in paired]
        s3_ratios = [pair["online_s3_root_call_ratio"] for pair in paired]
        baseline_runs = [points[label]["candidate"]["report"] for label in points]
        prepared_runs = [points[label]["prepared"]["report"] for label in points]
        curves.append({
            "n": n, "target_count": 12, "run_count": 24,
            "all_verified": True, "all_prepared_full_rank": True,
            "actual_factor_base_points": prepared_runs[0]["factor_base_points"],
            "folded_columns": prepared_runs[0]["orbit_columns"],
            "factor_base_digest": prepared_runs[0]["factor_base_digest"],
            "baseline_rank_attempts": sorted({run["rank_attempts"] for run in baseline_runs}),
            "prepared_rank_attempts": sorted({run["rank_attempts"] for run in prepared_runs}),
            "exploratory_geometric_mean_online_ratio": math.exp(statistics.mean(map(math.log, ratios))),
            "exploratory_median_online_ratio": statistics.median(ratios),
            "exploratory_min_online_ratio": min(ratios),
            "all_target_ratios_at_least_2": all(ratio >= 2 for ratio in ratios),
            "min_online_s3_root_call_ratio": min(s3_ratios),
            "median_online_s3_root_call_ratio": statistics.median(s3_ratios),
            "median_baseline_online_ms": statistics.median(pair["baseline_online_ms"] for pair in paired),
            "median_prepared_online_ms": statistics.median(pair["prepared_online_ms"] for pair in paired),
            "median_baseline_reusable_setup_ms": statistics.median(pair["baseline_reusable_setup_ms"] for pair in paired),
            "median_prepared_reusable_setup_ms": statistics.median(pair["prepared_reusable_setup_ms"] for pair in paired),
            "median_baseline_cold_ms": statistics.median(pair["baseline_cold_ms"] for pair in paired),
            "median_prepared_cold_ms": statistics.median(pair["prepared_cold_ms"] for pair in paired),
            "pairs": paired,
        })
    result = {"kind": "s3_prepared_local_exploratory_pair_summary_v1",
              "scope": "frozen one-target N41/N53 implementation comparison",
              "threshold": "at_least_2x_on_each_curve",
              "accounting": "target_point_input_outside_online; full_rank_base_logs_prepared_before_target_input",
              "strict_host_isolation_receipt": None, "controlled_speedup": None,
              "a_a_noise_control": None, "run_count": len(rows),
              "raw_ledger_sha256": sha(ledger), "raw_archive_sha256": sha(archive),
              "curves": curves}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"output": str(output), "raw_archive_sha256": result["raw_archive_sha256"],
                      "curves": [{"n": curve["n"], "min_ratio": curve["exploratory_min_online_ratio"],
                                  "geomean_ratio": curve["exploratory_geometric_mean_online_ratio"]}
                                 for curve in curves]}, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", required=True, type=Path)
    parser.add_argument("--archive", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    main(args.ledger, args.archive, args.output)
