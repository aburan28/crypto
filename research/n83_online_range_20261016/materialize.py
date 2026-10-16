#!/usr/bin/env python3
"""Bind the replayed expanded N83 online screen to exact IC1 candidates."""

import argparse
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RUN_ID = 38896107047
SOURCE_COMMIT = "c001093fb47304814bdee2a6b3d31ea8389d9213"
SOURCE_SHA256 = "3ee6bc11399b382e72789b152043e5a67abc5b783c4e998ad8c534d60baebc03"
ARTIFACT_ID = 12047173997
ARCHIVE_SHA256 = "53cdc2fa1ccc46a4617cf9ccc80609d0bc19a5e6d90b2b6a5d849d33ed8c0606"
ANALYSIS_SHA256 = "0e33e7ec12ad35f19ec15314455155a86f5137a998535d9144aa77d5c5859717"


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


dynamic = module("n83_range_candidate_helper",
                 ROOT / "research/n83_root11_dynamic_20261016/materialize.py")
panel = dynamic.panel
panel.SOURCE_COMMIT = SOURCE_COMMIT
panel.SOURCE_SHA256 = SOURCE_SHA256
runner = module("n83_expanded_range_runner", HERE / "run.py")


def write_once(path, value):
    blob = json.dumps(value, sort_keys=True, indent=2) + "\n"
    if path.exists() and path.read_text() != blob:
        raise ValueError(f"sealed record differs: {path}")
    path.write_text(blob)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("results", type=Path)
    parser.add_argument("--archive", required=True, type=Path)
    args = parser.parse_args()
    root = args.results.resolve()
    versioned = HERE / f"analysis_{RUN_ID}.json"
    if (panel.file_sha(args.archive) != ARCHIVE_SHA256
            or panel.file_sha(root / "analysis.json") != ANALYSIS_SHA256
            or panel.file_sha(versioned) != ANALYSIS_SHA256):
        raise ValueError("raw ZIP or independent analysis hash differs")
    analysis = panel.read(versioned)
    source = panel.read(root / "source.json")
    host = panel.read(root / "host.json")
    if (analysis["status"] != "BOUNDED_UNKNOWN"
            or analysis["evidence_status"] != "PASS_INDEPENDENT_REPLAY"
            or analysis["execution_order_status"] != "PASS"
            or analysis["selected_factor_base"] is not None
            or host["git_sha"] != SOURCE_COMMIT
            or source["source_sha256"] != SOURCE_SHA256
            or source["hashes"]["examples/koblitz_orbit_dlp_fast.rs"] != SOURCE_SHA256
            or source["hashes"]["Cargo.lock"] != panel.LOCK_SHA256
            or source["runner_sha256"] != panel.file_sha(HERE / "run.py")
            or source["arm_order"] != [name for name, _, _, _ in runner.ARMS]):
        raise ValueError("executed source, order or replay differs")
    workload = panel.read(ROOT / "research/n83_root11_dynamic_20261016"
                          / "ic1_rows_38876385086.json")["workload"]
    workload_id = workload["workload_id"]
    if workload_id != "eebf98ee1031":
        raise ValueError("exposed one-point workload differs")
    field, curve = panel.curve_record()
    by_arm = {item["arm"]: item for item in analysis["arms"]}
    if set(by_arm) != {name for name, _, _, _ in runner.ARMS}:
        raise ValueError("replayed arm set differs")
    manifests = {}
    identity_by_layout = {}
    repetitions = {}
    rows = []
    for name, k, layout, _ in runner.ARMS:
        item = by_arm[name]
        if (item["status"] != "PASS_INDEPENDENT_REPLAY"
                or item["resources"]["status"] != "PASS"):
            rows.append({"arm": name, "candidate_id": None,
                         "workload_id": workload_id, "run_id": None,
                         "status": item["status"], "raw_evidence": item})
            continue
        folder = root / "runs" / name
        summary = panel.read(folder / "stdout.jsonl")
        record = dynamic.candidate(k, layout, summary,
                                   folder / "base_dump.jsonl", field, curve)
        if (item["base_hash"] != summary["base_hash"]
                or item["actual_factor_base_points"] != 166 * k):
            raise ValueError(f"{name}: replayed base differs")
        candidate_id = record["candidate_id"]
        key = (k, layout)
        if key in identity_by_layout and identity_by_layout[key] != candidate_id:
            raise ValueError(f"K{k} {layout}: A/A candidate identities differ")
        identity_by_layout[key] = candidate_id
        manifests[candidate_id] = record
        repetitions[candidate_id] = repetitions.get(candidate_id, 0) + 1
        rows.append({
            "arm": name, "candidate_id": candidate_id,
            "workload_id": workload_id,
            "run_id": candidate_id + "W" + workload_id + "R"
                      + str(repetitions[candidate_id]),
            "status": item["status"], "hardware_tier": analysis["hardware_tier"],
            "source_commit": SOURCE_COMMIT,
            "executed_binary_sha256": source["binary_sha256"],
            "rho_binary_sha256": source["rho_binary_sha256"],
            "host_receipt_sha256": panel.file_sha(root / "host.json"),
            "source_receipt_sha256": panel.file_sha(root / "source.json"),
            "correctness_receipt_sha256": ANALYSIS_SHA256,
            "raw_artifact_run_id": RUN_ID,
            "raw_artifact_id": ARTIFACT_ID,
            "one_target_online_ms": item["one_target_online_ms"],
            "online_phases_ms": item["online_phases_ms"],
            "cold_process_ms_exploratory": item["cold_process_ms_exploratory"],
            "actual_usable_points": item["actual_factor_base_points"],
            "effective_columns": item["effective_columns"],
            "peak_rss_bytes": item["resources"]["peak_rss_bytes"],
        })
    for candidate_id, record in manifests.items():
        write_once(HERE / f"candidate_{candidate_id}.json", record)
    write_once(HERE / f"ic1_rows_{RUN_ID}.json", {
        "schema": "n83.online-range-ic1/v1",
        "workload": workload, "source_commit": SOURCE_COMMIT,
        "archive_sha256": ARCHIVE_SHA256,
        "analysis_sha256": ANALYSIS_SHA256,
        "candidates": sorted(manifests), "rows": rows,
        "selection_status": analysis["status"],
        "selected_candidate_id": None,
        "online_aa_control": analysis["online_aa_control"],
        "same_point_rho": analysis["rho"],
    })
    print(json.dumps({"candidate_count": len(manifests),
                      "row_count": len(rows), "selected_candidate_id": None},
                     sort_keys=True))


if __name__ == "__main__":
    main()
