#!/usr/bin/env python3
"""Seal independently replayed same-point Tier 1 IC/rho pairs."""

import argparse
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
RUNNER = REPO / "research/n83_unresolved_tier1_20261015/run.py"
SOURCE_COMMIT = "e1a77618587533364cf2236f09cc41b5bc816429"
RUN_ID = 38818664674

spec = importlib.util.spec_from_file_location("n83_same_point_runner", RUNNER)
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
    if path.exists():
        try:
            if gzip.decompress(path.read_bytes()) != blob:
                raise ValueError(f"sealed content differs: {path}")
        except (EOFError, OSError) as error:
            raise ValueError(f"invalid sealed gzip: {path}") from error
    else:
        path.write_bytes(gzip.compress(blob, compresslevel=9, mtime=0))
    return {"path": path.name, "sha256": digest(path), "bytes": path.stat().st_size}


def seal(root, artifact_zip, expected_zip, expected_analysis, artifact_id):
    hashes, point, _, candidates = runner.frozen()
    if artifact_id <= 0 or digest(artifact_zip) != expected_zip:
        raise ValueError("artifact ID or archived ZIP digest differs")
    analysis_path = root / "analysis.json"
    if digest(analysis_path) != expected_analysis:
        raise ValueError("independent replay digest differs")
    analysis, source, host = read(analysis_path), read(root / "source.json"), read(root / "host.json")
    names = [row["pipeline"] for row in candidates]
    rows = analysis["pairs"]
    if (analysis["schema"] != "n83.unresolved-tier1-analysis/v1"
            or analysis["factor_base_selection_status"] != "BOUNDED_UNKNOWN"
            or analysis["target_count"] != 1
            or analysis["candidate_order"] != names
            or len(rows) != 2 or [row["pipeline"] for row in rows] != names
            or source["schema"] != "n83.unresolved-tier1-source/v1"
            or source["hashes"] != hashes
            or source["candidate_order"] != candidates
            or source["runner_sha256"] != digest(RUNNER)
            or source["confirmation_selection_sha256"] != digest(runner.RESULT)
            or host["schema"] != "n83.unresolved-tier1-host/v1"
            or host["git_sha"] != SOURCE_COMMIT):
        raise ValueError("same-point source, host, workload or replay identity differs")
    complete = analysis["status"] == "PASS_TIER1_SAME_POINT_PAIR_SET"
    if complete != (analysis["same_allocation"] is True and
                   all(row["status"] == "PASS_TIER1_PAIR" for row in rows)):
        raise ValueError("pair completion differs from replay status")
    if analysis["status"] not in ("PASS_TIER1_SAME_POINT_PAIR_SET",
                                  "TIER1_GATE_INCOMPLETE"):
        raise ValueError("unknown Tier 1 replay status")
    workload = {
        "schema_version": 1, "curve_id": "EC1N83Ckb1a1haacb69ee9f37",
        "target_count": 1, "public_targets": [point],
        "target_sha256": runner.TARGET_SHA256,
        "answer_sidecar_sha256": runner.ANSWER_SHA256,
        "input_law": {"domain": "n83-unresolved-tier1-one-fresh-target-20261015-v1",
                      "answer_sidecar_sha256": runner.ANSWER_SHA256},
        "cold_or_warm": "two_fresh_IC_processes_with_empty_indexes_and_logs; one_independent_rho_solve_per_candidate_on_same_public_point",
    }
    workload_id = hashlib.sha256(json.dumps(workload, sort_keys=True,
        separators=(",", ":")).encode()).hexdigest()[:12]
    workload_record = gzip_once(HERE / "workload.json.gz", {"workload_id": workload_id, **workload})
    measurements = []
    for index, (candidate, pair) in enumerate(zip(candidates, rows), 1):
        passed = pair["status"] == "PASS_TIER1_PAIR"
        ic = pair.get("ic") or {}
        rho = pair.get("rho") or {}
        boundary_overhead_ms = None
        if passed:
            boundary_overhead_ms = ic["online_ms"] - sum(ic["phases_ms"].values())
            if not 0 <= boundary_overhead_ms <= 1:
                raise ValueError("IC online phase boundary overhead differs")
        measurements.append({
            "candidate_id": candidate["candidate_id"],
            "candidate_manifest_sha256": candidate["candidate_manifest_sha256"],
            "pipeline": candidate["pipeline"],
            "workload_id": workload_id,
            "run_id": candidate["candidate_id"] + "W" + workload_id + "R" + str(index),
            "status": pair["status"], "hardware_tier": "tier1_cloud_instance",
            "target_count": 1, "target_point": point,
            "pair_order": pair.get("order"),
            "candidate_factor_base_points": candidate["actual_usable_points"],
            "effective_columns": candidate["k"],
            "source_base_hash_blake3": candidate["source_base_hash_blake3"],
            "root_table_entries": candidate["root_table_entries"],
            "ic_online_ms": ic.get("online_ms") if passed else None,
            "rho_online_ms": rho.get("online_ms") if passed else None,
            "online_speedup_rho_over_ic": pair.get("ratio_rho_over_ic") if passed else None,
            "ic_online_phases_ms": ic.get("phases_ms") if passed else None,
            "ic_online_boundary_overhead_ms": boundary_overhead_ms,
            "ic_setup_ms": ic.get("setup_ms"),
            "ic_peak_rss_bytes": ic.get("peak_rss_bytes"),
            "rho_peak_rss_bytes": rho.get("peak_rss_bytes"),
            "noise_receipt": pair.get("noise"),
            "correctness": "PASS_INDEPENDENT_REPLAY" if passed else None,
            "ic_usage": ic.get("usage"), "rho_usage": rho.get("usage"),
        })
    measurement_record = gzip_once(HERE / "measurements.json.gz", measurements)
    write_once(HERE / "result.json", {
        "schema": "n83.unresolved-tier1-result/v1",
        "status": analysis["status"],
        "factor_base_selection_status": "BOUNDED_UNKNOWN",
        "candidate_order": names,
        "online_ratio_rho_over_ic": analysis["online_ratio_rho_over_ic"],
        "descriptive_ic_online_ms": analysis["descriptive_ic_online_ms"],
        "workload_id": workload_id,
        "artifact_run_id": RUN_ID, "artifact_id": artifact_id,
        "artifact_sha256": expected_zip, "analysis_sha256": expected_analysis,
        "source_commit": SOURCE_COMMIT,
    })
    write_once(HERE / "records_index.json", {
        "schema": "n83.unresolved-tier1-compressed-records/v1",
        "artifact_run_id": RUN_ID, "artifact_id": artifact_id,
        "artifact_sha256": expected_zip, "analysis_sha256": expected_analysis,
        "host_sha256": digest(root / "host.json"),
        "source_receipt_sha256": digest(root / "source.json"),
        "workload_id": workload_id, "workload": workload_record,
        "measurement_rows": len(measurements), "measurements": measurement_record,
    })
    print(json.dumps({"status": analysis["status"],
                      "online_ratio_rho_over_ic": analysis["online_ratio_rho_over_ic"],
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
