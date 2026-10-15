#!/usr/bin/env python3
"""Seal an independently replayed fresh N83 confirmation artifact."""

import argparse
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
RUNNER = REPO / "research/n83_edge_confirmation_20261015/run.py"
TARGETS = REPO / "research/n83_edge_confirmation_20261015/targets.jsonl"
SOURCE_COMMIT = "a369212178ffd28ec30e0e20eaa4729ec7b81e8d"
PRODUCER_SHA256 = "decd62935e518ff769cb1198b9da4d368f7b98ec8169d4f9e73c15d0e902ece9"
RUN_ID = 38809654754

spec = importlib.util.spec_from_file_location("n83_confirmation_runner", RUNNER)
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text())


def write_once(path, value):
    blob = json.dumps(value, sort_keys=True, indent=2).encode() + b"\n"
    if path.exists() and path.read_bytes() != blob:
        raise ValueError(f"sealed file differs: {path}")
    if not path.exists():
        path.write_bytes(blob)
    return digest(path)


def gzip_once(path, value):
    blob = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    packed = gzip.compress(blob, compresslevel=9, mtime=0)
    if path.exists() and path.read_bytes() != packed:
        raise ValueError(f"sealed file differs: {path}")
    if not path.exists():
        path.write_bytes(packed)
    return {"path": path.name, "sha256": digest(path), "bytes": path.stat().st_size}


def seal(root, artifact_zip, expected_zip, expected_analysis, artifact_id):
    admission, hashes = runner.frozen()
    if artifact_id <= 0 or digest(artifact_zip) != expected_zip:
        raise ValueError("artifact ID or archived ZIP digest differs")
    analysis_path = root / "analysis.json"
    if digest(analysis_path) != expected_analysis:
        raise ValueError("independent replay digest differs")
    analysis, source, host = read(analysis_path), read(root / "source.json"), read(root / "host.json")
    candidates = admission["candidate_order"]
    order = [arm for arm, _ in runner.arms(candidates)]
    if (analysis["schema"] != "n83.edge-confirmation-analysis/v1"
            or analysis["candidate_order"] != [row["pipeline"] for row in candidates]
            or analysis["arm_order"] != order
            or len(analysis["arms"]) != len(order)
            or {row["arm"] for row in analysis["arms"]} != set(order)
            or analysis["target_count"] != 4000
            or source["schema"] != "n83.edge-confirmation-source/v1"
            or source["hashes"] != hashes
            or source["binary_sha256"] != PRODUCER_SHA256
            or source["panel_runner_sha256"] != digest(RUNNER)
            or source["admission_sha256"] != digest(runner.ADMISSION)
            or source["arm_order"] != order
            or source["candidate_ids"] != [row["candidate_id"] for row in candidates]
            or host["git_sha"] != SOURCE_COMMIT):
        raise ValueError("source, admission, host or replay identity differs")
    balanced = analysis["balanced_comparison"]
    if (analysis["status"] == "COMPLETE") != (balanced is not None):
        raise ValueError("comparison completion status differs")
    selected = analysis["selected_pipeline"]
    if (selected != (balanced["selected_pipeline"] if balanced else None)
            or (selected is not None and selected not in
                {row["pipeline"] for row in candidates})
            or analysis["decision_status"] !=
                ("SELECTED" if selected else "BOUNDED_UNKNOWN")):
        raise ValueError("selection status differs from frozen comparison")

    points = [json.loads(line) for line in TARGETS.read_text().splitlines() if line]
    if len(points) != 4000 or len({tuple(point) for point in points}) != 4000:
        raise ValueError("frozen public-point panel differs")
    workload = {
        "schema_version": 1, "curve_id": "EC1N83Ckb1a1haacb69ee9f37",
        "target_count": 4000, "public_targets": points,
        "target_sha256": admission["target_sha256"],
        "answer_sidecar_sha256": admission["answer_sidecar_sha256"],
        "input_law": {
            "domain": "n83-edge-confirm-fresh-4000-20261015-v1",
            "answer_sidecar_sha256": admission["answer_sidecar_sha256"]},
        "cold_or_warm": "four_mirrored_fresh_processes_with_empty_index_and_logs; reusable_setup_outside_per_target_online_intervals",
    }
    workload_id = hashlib.sha256(json.dumps(workload, sort_keys=True,
        separators=(",", ":")).encode()).hexdigest()[:12]
    workload_record = gzip_once(HERE / "workload.json.gz", {"workload_id": workload_id, **workload})
    by_arm = {row["arm"]: row for row in analysis["arms"]}
    measurements = []
    for index, (arm, candidate) in enumerate(runner.arms(candidates), 1):
        row = by_arm[arm]
        folder = root / "runs" / arm
        usage = read(folder / "usage.json") if (folder / "usage.json").exists() else None
        phase_rows = None
        final_rank = None
        final_la_ms = None
        if row["status"] == "PASS_INDEPENDENT_REPLAY":
            targets = [json.loads(line) for line in
                       (folder / "targets_result.jsonl").read_text().splitlines() if line]
            if len(targets) != 4000 or len(row["targets"]) != 4000:
                raise ValueError(f"{arm}: verified target count differs")
            keys = ("target_query_ms", "target_pdp_ms", "target_relation_check_ms",
                    "target_descent_ms", "target_recovery_check_ms")
            phase_rows = [{key: target[key] for key in keys} for target in targets]
            summary = read(folder / "stdout.jsonl")
            final_rank = summary["rank"]
            final_la_ms = summary["timing_ms"]["linear_algebra"]
        measurements.append({
            "arm": arm, "candidate_id": candidate["candidate_id"],
            "workload_id": workload_id,
            "run_id": candidate["candidate_id"] + "W" + workload_id + "R" + str(index),
            "status": row["status"], "hardware_tier": "cloud_instance_exploratory",
            "workload_role": "secondary_fresh_4000_point_factor_base_confirmation",
            "controlled_speedup": None,
            "candidate_manifest_sha256": candidate["candidate_manifest_sha256"],
            "actual_factor_base_points": row.get("actual_factor_base_points"),
            "effective_columns": row.get("effective_columns"),
            "base_hash": row.get("base_hash"),
            "root_table_entries": row.get("root_table_entries"),
            "peak_rss_bytes": row.get("peak_rss_bytes"),
            "resource_envelope_status": row.get("resource_envelope_status"),
            "cpu_throttled_periods_delta": row.get("cpu_throttled_periods_delta"),
            "oom_kills_delta": row.get("oom_kills_delta"),
            "system_steal_ticks_delta": row.get("system_steal_ticks_delta"),
            "setup_ms": row.get("setup_ms"),
            "final_rank": final_rank, "final_la_ms": final_la_ms,
            "mean_online_ms_exploratory": row.get("mean_online_ms_exploratory"),
            "median_online_ms_exploratory": row.get("median_online_ms_exploratory"),
            "target_records": row.get("targets"),
            "target_online_phases_ms": phase_rows,
            "full_process_wall_seconds": usage.get("wall_seconds") if usage else None,
            "stop_reason": usage.get("stop_reason") if usage else None,
            "exit_status": usage.get("exit_status") if usage else None,
        })
    measurement_record = gzip_once(HERE / "measurements.json.gz", measurements)
    write_once(HERE / "selection.json", {
        "schema": "n83.edge-confirmation-selection/v1",
        "status": analysis["decision_status"] if analysis["status"] == "COMPLETE" else "INCOMPLETE",
        "selected_pipeline": selected,
        "selected_candidate_id": balanced["selected_candidate_id"] if balanced else None,
        "balanced_comparison": balanced,
        "admission_sha256": digest(runner.ADMISSION),
        "workload_id": workload_id,
        "artifact_run_id": RUN_ID, "artifact_id": artifact_id,
        "artifact_sha256": expected_zip, "analysis_sha256": expected_analysis,
        "source_commit": SOURCE_COMMIT,
    })
    write_once(HERE / "records_index.json", {
        "schema": "n83.edge-confirmation-compressed-records/v1",
        "artifact_run_id": RUN_ID, "artifact_id": artifact_id,
        "artifact_sha256": expected_zip, "analysis_sha256": expected_analysis,
        "host_sha256": digest(root / "host.json"),
        "source_receipt_sha256": digest(root / "source.json"),
        "workload_id": workload_id, "workload": workload_record,
        "measurement_rows": len(measurements), "measurements": measurement_record,
    })
    print(json.dumps({"status": analysis["status"],
                      "decision_status": analysis["decision_status"],
                      "selected_pipeline": selected,
                      "workload_id": workload_id}, sort_keys=True))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifact-zip", required=True, type=Path)
    parser.add_argument("--artifact-sha256", required=True)
    parser.add_argument("--analysis-sha256", required=True)
    parser.add_argument("--artifact-id", required=True, type=int)
    parser.add_argument("artifact_results", type=Path)
    args = parser.parse_args()
    seal(args.artifact_results.resolve(), args.artifact_zip.resolve(),
         args.artifact_sha256, args.analysis_sha256, args.artifact_id)


if __name__ == "__main__":
    main()
