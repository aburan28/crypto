#!/usr/bin/env python3
"""Bind the independently replayed K4096 diagnostic to an IC1 identity."""

import argparse
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ARCHIVE_SHA256 = "41b7be597432218a38dc74081f09c38bfa0d26f4d10570f213ed6cbcde5bd289"
ANALYSIS_SHA256 = "8da0af76705063e04ef4dd82309fd007ba13df0d9dd1068490028be8fe5b08ca"

spec = importlib.util.spec_from_file_location(
    "n83_id30_k4096_materializer", HERE / "materialize.py")
common = importlib.util.module_from_spec(spec)
spec.loader.exec_module(common)
panel = common.panel


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("results", type=Path)
    parser.add_argument("--archive", type=Path, required=True)
    args = parser.parse_args()
    output = args.results.resolve()
    if (panel.file_sha(args.archive) != ARCHIVE_SHA256
            or panel.file_sha(output / "analysis.json") != ANALYSIS_SHA256
            or panel.file_sha(HERE / "analysis_k4096_20261016.json") != ANALYSIS_SHA256):
        raise ValueError("exact K4096 archive or replay digest differs")
    analysis = panel.read(output / "analysis.json")
    source = panel.read(output / "source.json")
    host = panel.read(output / "host.json")
    item = analysis["arm"]
    if (analysis["status"] != "PASS_EXACT_K4096"
            or item["status"] != "PASS_INDEPENDENT_REPLAY"
            or analysis["performance_claim"] is not None
            or source["source_sha256"] != common.SOURCE_SHA256
            or source["runner_sha256"] != panel.file_sha(
                HERE / "evidence_v1/run_k4096.py.txt")
            or host["physical_memory_bytes"] < 40 * 1024**3
            or source["k"] != 4096
            or item["actual_factor_base_points"] != 679936):
        raise ValueError("K4096 source, host or replay status differs")
    folder = output / "runs/opt_in_id30"
    summary = panel.read(folder / "stdout.jsonl")
    field, curve = panel.curve_record()
    record = common.candidate("opt_in_id30", summary,
                              folder / "base_dump.jsonl", field, curve, 4096)
    if (record["factor_base"]["source_base_hash_blake3"] != item["base_hash"]
            or item["root_table_entries"] != summary["root_table_entries"]):
        raise ValueError("K4096 base or exact index identity differs")
    workload = panel.read(HERE / "ic1_rows_control_20261016.json")["workload"]
    workload_id = workload["workload_id"]
    if workload_id != "eebf98ee1031":
        raise ValueError("one-point workload identity differs")
    candidate_id = record["candidate_id"]
    common.write_once(HERE / f"candidate_{candidate_id}.json", record)
    common.write_once(HERE / "ic1_row_k4096_20261016.json", {
        "schema": "n83.id30-k4096-ic1/v1",
        "workload": workload,
        "selected_candidate_id": None,
        "row": {
            "candidate_id": candidate_id,
            "workload_id": workload_id,
            "run_id": candidate_id + "W" + workload_id + "R1",
            "status": item["status"],
            "hardware_tier": "local_correctness_diagnostic",
            "source_sha256": common.SOURCE_SHA256,
            "executed_binary_sha256": source["binary_sha256"],
            "raw_archive_sha256": ARCHIVE_SHA256,
            "analysis_sha256": ANALYSIS_SHA256,
            "actual_usable_points": item["actual_factor_base_points"],
            "effective_columns": 4096,
            "base_hash": item["base_hash"],
            "root_table_entries": item["root_table_entries"],
            "one_target_online_ms_exploratory": item["online_ms_exploratory"],
            "online_phases_ms": item["online_phases_ms"],
            "outer_cold_wall_seconds_exploratory": item["cold_wall_seconds_exploratory"],
            "peak_sampled_rss_bytes": analysis["peak_sampled_rss_bytes"],
            "min_sampled_free_percent": analysis["min_sampled_free_percent"],
            "recovered_scalar": item["recovered_scalar"],
            "point_indices": item["point_indices"],
            "performance_claim": None,
        },
    })
    print(json.dumps({"candidate_id": candidate_id,
                      "workload_id": workload_id,
                      "actual_usable_points": item["actual_factor_base_points"]},
                     sort_keys=True))


if __name__ == "__main__":
    main()
