#!/usr/bin/env python3
"""Verify the preserved C2800 screen stopped at its frozen RSS gate."""

import gzip
import hashlib
import json
import tarfile
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FREEZE = HERE / "SCREEN_FREEZE.json"
LEDGER = HERE / "local/raw_partial3.jsonl.gz"
REPORTS = HERE / "local/raw_reports.tar.gz"
HOST = HERE / "local/host.json"


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def main():
    freeze = json.loads(FREEZE.read_bytes())
    host = json.loads(HOST.read_bytes())
    rows = [json.loads(line) for line in gzip.decompress(LEDGER.read_bytes()).splitlines()]
    assert freeze["target_labels"] == [f"T{i:03}" for i in range(149, 157)]
    assert len(rows) == 3 and host["host_exclusive_partition_receipt"] is None
    assert [(row["label"], row["arm"], row["repetition"]) for row in rows] == [
        ("T149", "baseline", 1), ("T149", "expanded", 1),
        ("T150", "expanded", 1),
    ]
    parent = json.loads((ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json").read_bytes())
    targets = {item["label"]: item for item in parent["targets"] if item["n"] == 53}
    with tarfile.open(REPORTS, "r:gz") as archive:
        assert set(archive.getnames()) == {
            "T149_baseline_R1.json", "T149_expanded_R1.json", "T150_expanded_R1.json"
        }
        for row in rows:
            frozen = freeze["arms"][row["arm"]]
            target = targets[row["label"]]
            assert row["candidate_id"] == frozen["candidate_id"]
            assert row["workload_id"] == target["workload_id"]
            assert row["point_sha256"] == target["point_sha256"]
            assert row["run_id"] == f"{frozen['candidate_id']}W{target['workload_id']}R1"
            assert row["source_sha256"] == freeze["source_sha256"]
            assert row["binary_sha256"] == freeze["binary_sha256"]
            name = f"{row['label']}_{row['arm']}_R1.json"
            with archive.extractfile(name) as member:
                raw = member.read()
            assert (sha(raw) if raw else None) == row["raw_sha256"]
            if row["label"] == "T150":
                assert row["row_origin"] == "operator_stop_record"
                assert row["error"] == "operator_stop_after_prior_memory_envelope_breach"
                assert row["report"] is None and not raw
                assert row["exit_code"] is None and not row["verified_native_and_fixture"]
                continue
            assert row["exit_code"] == 0 and row["error"] is None
            assert row["verified_native_and_fixture"] and row["phase_ok"]
            report = json.loads(raw)
            assert report == row["report"]
            assert report["group_verified"] and report["recovered_scalar"] == target["fixture_scalar"]
            assert report["target"] == target["public_point"]
            assert report["factor_base_digest"] == frozen["factor_base_digest"]
            assert report["factor_base_points"] == frozen["factor_base_points"]
            assert report["prepared_rank"] == frozen["prepare_rank"]
            assert report["root_table_slots"] == frozen["root_table_slots"]
            assert report["target_search_schedule"] == frozen["target_search_schedule"]
            assert report["fallback_attempts"] == 0

    baseline, expanded, interrupted = rows
    envelope = freeze["memory_envelope_bytes"]
    baseline_peak = baseline["whole_process_resource_counters"]["peak_rss_bytes"]
    expanded_peak = expanded["whole_process_resource_counters"]["peak_rss_bytes"]
    assert baseline_peak <= envelope < expanded_peak
    assert baseline["resource_within_envelope"] is True
    assert expanded["resource_within_envelope"] is False
    assert interrupted["resource_within_envelope"] is None
    assert (HERE / "local/runner_exit.txt").read_text().strip() == "143"
    result = {
        "kind": "s3_n53_c2800prefix_paired_screen_resource_stop_analysis_v1",
        "freeze_sha256": sha(FREEZE.read_bytes()),
        "ledger_decompressed_sha256": sha(gzip.decompress(LEDGER.read_bytes())),
        "raw_reports_archive_sha256": sha(REPORTS.read_bytes()),
        "planned_target_count": 8,
        "planned_run_count": 16,
        "completed_verified_run_count": 2,
        "interrupted_run_count": 1,
        "unstarted_run_count": 13,
        "memory_envelope_bytes": envelope,
        "baseline_peak_rss_bytes": baseline_peak,
        "expanded_peak_rss_bytes": expanded_peak,
        "expanded_rss_excess_bytes": expanded_peak - envelope,
        "t149_baseline_target_s3_calls": baseline["report"]["target_s3_calls"],
        "t149_expanded_target_s3_calls": expanded["report"]["target_s3_calls"],
        "t149_paired_call_ratio": (baseline["report"]["target_s3_calls"]
                                   / expanded["report"]["target_s3_calls"]),
        "t149_paired_online_ratio_exploratory": (
            baseline["report"]["timing_ms"]["target_online_after_reusable_setup"]
            / expanded["report"]["timing_ms"]["target_online_after_reusable_setup"]),
        "screen_advanced": False,
        "screen_aggregate_call_ratio": None,
        "screen_aggregate_online_ratio": None,
        "host_exclusive_partition_receipt": None,
        "controlled_speedup": None,
    }
    (HERE / "analysis.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
