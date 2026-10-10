#!/usr/bin/env python3
"""Rebuild the frozen C1200/C2800 direct-fill one-repetition screen."""

import gzip
import hashlib
import json
import math
import random
import statistics
import tarfile
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FREEZE = HERE / "SCREEN_FREEZE.json"
LEDGER = HERE / "local/raw_paired.jsonl.gz"
REPORTS = HERE / "local/raw_reports.tar.gz"
HOST = HERE / "local/host.json"


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def geomean(values):
    return math.exp(statistics.mean(math.log(v) for v in values))


def percentile(values, fraction):
    ordered = sorted(values)
    position = (len(ordered) - 1) * fraction
    lower = int(position)
    return ordered[lower] + (ordered[min(lower + 1, len(ordered) - 1)] - ordered[lower]) * (position - lower)


def main():
    freeze = json.loads(FREEZE.read_text())
    host = json.loads(HOST.read_text())
    with gzip.open(LEDGER, "rt") as source:
        rows = [json.loads(line) for line in source]
    assert 1 <= len(rows) <= 16
    assert host["host_exclusive_partition_receipt"] is None
    expected_order = []
    for index, label in enumerate(freeze["target_labels"]):
        arms = (("baseline", "expanded") if index % 2 == 0
                else ("expanded", "baseline"))
        expected_order.extend((label, arm, 1) for arm in arms)
    assert [(row["label"], row["arm"], row["repetition"]) for row in rows] == expected_order[:len(rows)]
    by_key = {}
    failed_runs = []
    with tarfile.open(REPORTS, "r:gz") as archive:
        archive_names = set(archive.getnames())
        for row in rows:
            key = (row["label"], row["arm"], row["repetition"])
            assert key not in by_key
            by_key[key] = row
            assert row["label"] in freeze["target_labels"]
            assert row["arm"] in freeze["arms"]
            assert row["repetition"] == 1
            assert row["candidate_id"] == freeze["arms"][row["arm"]]["candidate_id"]
            assert row["binary_sha256"] == freeze["binary_sha256"]
            assert row["source_sha256"] == freeze["source_sha256"]
            name = f"{row['label']}_{row['arm']}_R{row['repetition']}.json"
            if name in archive_names:
                with archive.extractfile(name) as member:
                    raw = member.read()
                assert (sha256(raw) if raw else None) == row["raw_sha256"]
                if row["report"] is not None:
                    assert json.loads(raw) == row["report"]
            else:
                assert row["raw_sha256"] is None and row["report"] is None
            valid = (row["exit_code"] == 0 and row["error"] is None
                     and row["verified_native_and_fixture"] and row["phase_ok"]
                     and row["resource_within_envelope"] is True)
            if not valid:
                failed_runs.append({"label": row["label"], "arm": row["arm"],
                                    "repetition": row["repetition"],
                                    "exit_code": row["exit_code"], "error": row["error"],
                                    "verified_native_and_fixture": row["verified_native_and_fixture"],
                                    "resource_within_envelope": row["resource_within_envelope"]})
                continue
            report = row["report"]
            assert report["group_verified"] and report["target_relation_found"]
            assert report["final_rank_complete"]
            assert report["prepared_rank"] == freeze["arms"][row["arm"]]["prepare_rank"]
            assert report["factor_base_digest"] == freeze["arms"][row["arm"]]["factor_base_digest"]
            assert report["factor_base_points"] == freeze["arms"][row["arm"]]["factor_base_points"]
            assert report["target_search_schedule"] == freeze["arms"][row["arm"]]["target_search_schedule"]
            assert report["target_search_prefix_columns"] == freeze["arms"][row["arm"]]["target_search_prefix_columns"]
            assert report["target_count"] == 1 and report["target_query_workers"] == 1
            assert report["timing_ms"]["target_online_after_reusable_setup"] > 0
            assert report["target_s3_calls"] > 0
            counters = row["whole_process_resource_counters"]
            assert counters["peak_rss_bytes"] > 0
            assert counters["peak_rss_bytes"] <= freeze["memory_envelope_bytes"]
    assert len(by_key) == len(rows)
    runner_exit = (HERE / "local/runner_exit.txt").read_text().strip()
    stop_file = HERE / "local/stop.json"
    stop = json.loads(stop_file.read_bytes()) if stop_file.exists() else None
    if stop is not None:
        assert stop["rows_started"] == len(rows)
        assert stop["planned_rows"] == 16
        assert failed_runs and stop["label"] == rows[-1]["label"]
        assert stop["arm"] == rows[-1]["arm"]
    complete = len(rows) == 16 and not failed_runs and stop is None and runner_exit == "0"
    if not complete:
        result = {
            "kind": "s3_n53_c2800direct_paired_call_screen_analysis_v1",
            "freeze_sha256": sha256(FREEZE.read_bytes()),
            "ledger_decompressed_sha256": sha256(gzip.decompress(LEDGER.read_bytes())),
            "raw_reports_archive_sha256": sha256(REPORTS.read_bytes()),
            "planned_target_count": len(freeze["target_labels"]),
            "planned_run_count": 16,
            "run_count": len(rows),
            "unstarted_run_count": 16 - len(rows),
            "native_and_fixture_verified_count": sum(bool(row["verified_native_and_fixture"]) for row in rows),
            "resource_within_envelope_count": sum(row["resource_within_envelope"] is True for row in rows),
            "failed_runs": failed_runs,
            "runner_exit": runner_exit,
            "stop": stop,
            "root_call_geomean_baseline_over_expanded": None,
            "online_geomean_baseline_over_expanded": None,
            "screen_advanced": False,
            "host_exclusive_partition_receipt": None,
            "controlled_speedup": None,
        }
        (HERE / "analysis.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
        print(json.dumps(result, indent=2, sort_keys=True))
        return

    target_rows = []
    call_ratios = []
    wall_ratios = []
    for label in freeze["target_labels"]:
        baseline = [by_key[label, "baseline", 1]]
        expanded = [by_key[label, "expanded", 1]]
        point_hashes = {row["point_sha256"] for row in baseline + expanded}
        workload_ids = {row["workload_id"] for row in baseline + expanded}
        assert len(point_hashes) == len(workload_ids) == 1
        calls = {}
        walls = {}
        peaks = {}
        for arm, pair in (("baseline", baseline), ("expanded", expanded)):
            arm_calls = {row["report"]["target_s3_calls"] for row in pair}
            assert len(arm_calls) == 1
            calls[arm] = arm_calls.pop()
            wall_samples = [row["report"]["timing_ms"]["target_online_after_reusable_setup"] for row in pair]
            walls[arm] = geomean(wall_samples)
            peaks[arm] = max(row["whole_process_resource_counters"]["peak_rss_bytes"] for row in pair)
        call_ratio = calls["baseline"] / calls["expanded"]
        wall_ratio = walls["baseline"] / walls["expanded"]
        call_ratios.append(call_ratio)
        wall_ratios.append(wall_ratio)
        target_rows.append({
            "label": label,
            "point_sha256": point_hashes.pop(),
            "workload_id": workload_ids.pop(),
            "baseline_target_s3_calls": calls["baseline"],
            "expanded_target_s3_calls": calls["expanded"],
            "baseline_online_ms": walls["baseline"],
            "expanded_online_ms": walls["expanded"],
            "baseline_peak_rss_bytes": peaks["baseline"],
            "expanded_peak_rss_bytes": peaks["expanded"],
            "baseline_over_expanded_call_ratio": call_ratio,
            "baseline_over_expanded_online_ratio": wall_ratio,
        })

    rng = random.Random(20261009)
    bootstrap = {"call": [], "online": []}
    for _ in range(20000):
        indices = [rng.randrange(len(target_rows)) for _ in target_rows]
        bootstrap["call"].append(geomean([call_ratios[i] for i in indices]))
        bootstrap["online"].append(geomean([wall_ratios[i] for i in indices]))
    call_geomean = geomean(call_ratios)
    online_geomean = geomean(wall_ratios)
    result = {
        "kind": "s3_n53_c2800direct_paired_call_screen_analysis_v1",
        "freeze_sha256": sha256(FREEZE.read_bytes()),
        "ledger_decompressed_sha256": sha256(gzip.decompress(LEDGER.read_bytes())),
        "raw_reports_archive_sha256": sha256(REPORTS.read_bytes()),
        "target_count": len(target_rows),
        "run_count": len(rows),
        "unstarted_run_count": 0,
        "runner_exit": runner_exit,
        "native_and_fixture_verified_count": len(rows),
        "resource_within_envelope_count": len(rows),
        "peak_rss_bytes_max": max(row["whole_process_resource_counters"]["peak_rss_bytes"] for row in rows),
        "failed_runs": [],
        "root_call_geomean_baseline_over_expanded": call_geomean,
        "root_call_bootstrap95_target_resample": [percentile(bootstrap["call"], 0.025), percentile(bootstrap["call"], 0.975)],
        "root_call_targets_at_least_2x": sum(r >= 2 for r in call_ratios),
        "online_geomean_baseline_over_expanded": online_geomean,
        "online_bootstrap95_target_resample": [percentile(bootstrap["online"], 0.025), percentile(bootstrap["online"], 0.975)],
        "online_targets_at_least_2x": sum(r >= 2 for r in wall_ratios),
        "repetitions_per_arm_per_target": 1,
        "screen_advance_root_call_geomean_ratio_min": freeze["screen_advance_root_call_geomean_ratio_min"],
        "screen_advanced": call_geomean >= freeze["screen_advance_root_call_geomean_ratio_min"],
        "host_exclusive_partition_receipt": None,
        "controlled_speedup": None,
        "targets": target_rows,
    }
    (HERE / "analysis.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: value for key, value in result.items() if key != "targets"}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
