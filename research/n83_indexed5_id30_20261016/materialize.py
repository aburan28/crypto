#!/usr/bin/env python3
"""Materialize the replayed N83 ID29/ID30 control as exact IC1 records."""

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE_SHA256 = "f14d9747cd24627d39d7cae4c9bef69ed2819544125ef15f31306ce5d75f3fc7"
ARCHIVE_SHA256 = "f831bdb50bff9dd27e47c5e21554aed6358d42b5d62859e10cc4ba147c85ddf7"
ANALYSIS_SHA256 = "69f77eaaae6e08185046c1701450e03e032c5892e6644dc370e2d15c07e58552"


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


dynamic = module("n83_id30_candidate_helper", ROOT / "research/n83_root11_dynamic_20261016/materialize.py")
panel = dynamic.panel
panel.SOURCE_SHA256 = SOURCE_SHA256


def write_once(path, value):
    blob = json.dumps(value, sort_keys=True, indent=2) + "\n"
    if path.exists() and path.read_text() != blob:
        raise ValueError(f"sealed identity differs: {path}")
    path.write_text(blob)


def candidate(name, summary, dump_path, field, curve, k=1024):
    if (summary["root_table_indexed5_id30_enabled"] is not
            (name == "opt_in_id30")):
        raise ValueError(f"{name}: executed ID mode differs")
    template = dict(summary)
    template["root_table_layout"] = "indexed83_5byte_exact"
    record = dynamic.candidate(k, "root11_dynamic", template,
                               dump_path, field, curve)
    record["implementation"]["source_commit"] = None
    record["implementation"]["source_sha256"]["ic"] = SOURCE_SHA256
    record["implementation"]["algorithm_flags"]["indexed5_id30"] = (
        name == "opt_in_id30")
    slot = record["point_decomposition"]["root_slot"]
    slot["state_id_bits"] = 30 if name == "opt_in_id30" else 29
    slot["fingerprint_bits"] = 9 if name == "opt_in_id30" else 10
    record["point_decomposition"]["variant"] += (
        "_state_id30_fingerprint9" if name == "opt_in_id30"
        else "_state_id29_fingerprint10")
    for section in ("point_decomposition", "relation_collection",
                    "relation_linear_algebra", "target_descent"):
        key = "solver_source_sha256" if section == "point_decomposition" else "source_digest"
        record[section][key] = SOURCE_SHA256
    body = {key: value for key, value in record.items() if key != "candidate_id"}
    record["candidate_id"] = (
        f"IC1N83Ckb1a1fb{166*k}PDP4rootRCguidedLAgaussTDdirectISO0h"
        + panel.digest(body)[:12])
    return record


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("results", type=Path)
    parser.add_argument("--archive", type=Path, required=True)
    args = parser.parse_args()
    output = args.results.resolve()
    analysis_path = output / "analysis.json"
    historical_source = subprocess.check_output(
        ["git", "show", "5a88a28299de2b1ae82a3a1a6a6a1cd6ed2dc076:examples/koblitz_orbit_dlp_fast.rs"],
        cwd=ROOT)
    if (panel.file_sha(args.archive) != ARCHIVE_SHA256
            or panel.file_sha(analysis_path) != ANALYSIS_SHA256
            or panel.file_sha(HERE / "analysis_control_20261016.json") != ANALYSIS_SHA256
            or hashlib.sha256(historical_source).hexdigest() != SOURCE_SHA256):
        raise ValueError("raw archive, replay or source digest differs")
    analysis = panel.read(analysis_path)
    source = panel.read(output / "source.json")
    if (analysis["status"] != "PASS_EXACT_PARITY"
            or analysis["same_base_index_and_scalar"] is not True
            or source["source_sha256"] != SOURCE_SHA256
            or source["runner_sha256"] != panel.file_sha(
                HERE / "evidence_v1/run_control.py.txt")
            or source["arm_order"] != ["legacy_id29", "opt_in_id30"]):
        raise ValueError("control identity or independent replay differs")
    targets = ROOT / "research/n83_cold_factor_base_20261016/targets.jsonl"
    answers = ROOT / "research/n83_cold_factor_base_20261016/answers.jsonl"
    point = panel.read(targets)
    answer = panel.read(answers)
    if answer["point"] != point or answer["fixture_index"] != 0:
        raise ValueError("frozen public point differs")
    workload = {
        "schema_version": 1, "curve_id": panel.CURVE_ID,
        "target_count": 1, "public_targets": [point],
        "input_law": {
            "domain": "n83-cold-factor-base-one-plus-four-20261016-v1",
            "fixture_index": 0,
            "scalar_rule": "1+(little_endian(SHA256(domain||u32_le(index))) mod (r-1))",
            "seed_sha256": answer["seed_sha256"],
            "generator_sha256": source["hashes"][
                "research/n83_cold_factor_base_20261016/generate.py"],
            "answer_sidecar_sha256": panel.file_sha(answers)},
        "target_sha256": panel.file_sha(targets),
        "cold_or_warm": "one_public_point_per_fresh_process_empty_base_index_and_logs",
    }
    workload_id = panel.digest(workload)[:12]
    if workload_id != "eebf98ee1031":
        raise ValueError("one-point workload identity differs")
    field, curve = panel.curve_record()
    rows = []
    candidates = []
    for item in analysis["arms"]:
        name = item["arm"]
        if item["status"] != "PASS_INDEPENDENT_REPLAY":
            raise ValueError(f"{name}: independent replay required")
        folder = output / "runs" / name
        summary = panel.read(folder / "stdout.jsonl")
        record = candidate(name, summary, folder / "base_dump.jsonl", field, curve)
        if (item["base_hash"] != summary["base_hash"]
                or item["root_table_entries"] != summary["root_table_entries"]):
            raise ValueError(f"{name}: base or index receipt differs")
        candidate_id = record["candidate_id"]
        candidates.append(candidate_id)
        write_once(HERE / f"candidate_{candidate_id}.json", record)
        rows.append({
            "arm": name, "candidate_id": candidate_id,
            "workload_id": workload_id,
            "run_id": candidate_id + "W" + workload_id + "R1",
            "status": item["status"], "hardware_tier": "local_correctness_control",
            "source_sha256": SOURCE_SHA256,
            "executed_binary_sha256": source["binary_sha256"],
            "raw_archive_sha256": ARCHIVE_SHA256,
            "analysis_sha256": ANALYSIS_SHA256,
            "actual_usable_points": item["actual_factor_base_points"],
            "effective_columns": 1024,
            "base_hash": item["base_hash"],
            "root_table_entries": item["root_table_entries"],
            "one_target_online_ms_exploratory": item["online_ms_exploratory"],
            "online_phases_ms": item["online_phases_ms"],
            "outer_cold_wall_seconds_exploratory": item["cold_wall_seconds_exploratory"],
            "recovered_scalar": item["recovered_scalar"],
            "point_indices": item["point_indices"],
            "performance_claim": None,
        })
    write_once(HERE / "ic1_rows_control_20261016.json", {
        "schema": "n83.id30-control-ic1/v1",
        "workload": {"workload_id": workload_id, **workload},
        "candidates": candidates, "rows": rows,
        "selection_status": analysis["status"],
        "selected_candidate_id": None,
    })
    print(json.dumps({"candidate_count": len(candidates),
                      "workload_id": workload_id, "row_count": len(rows)},
                     sort_keys=True))


if __name__ == "__main__":
    main()
